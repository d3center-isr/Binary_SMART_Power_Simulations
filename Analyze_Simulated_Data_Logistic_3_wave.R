analyze_simulated_data_logistic_3_wave <- function(sim_data_wide, 
                                          regimen1,
                                          regimen2) { 
  # Define regimens for estimands:
  regimens_1wave <- rbind(mm=c(1,-1,-1,1),
                                 mp=c(1,-1,1,-1),
                                 pm=c(1,1,-1,-1),
                                 pp=c(1,1,1,1));
  colnames(regimens_1wave) <- c("int","a1","a2","a1a2");
  contrast_coefficients <- regimens_1wave[regimen1,]-
    regimens_1wave[regimen2,];
  dummy_contrast_coefficients <- c(mm=-1,mp=0,pm=1,pp=0);
  sim_data_wide$known_weight <- 4*(sim_data_wide$R==0)+
    2*(sim_data_wide$R==1);
  sim0 <- sim_data_wide; sim0$time <- 0; sim0$Y <- sim0$Y0;
  sim1 <- sim_data_wide; sim1$time <- 1; sim1$Y <- sim1$Y1;
  sim2 <- sim_data_wide; sim2$time <- 2; sim2$Y <- sim2$Y2;

  sim_long <- rbind(sim0[,c("id","time","A1","R","A2","Y","known_weight")],
                    sim1[,c("id","time","A1","R","A2","Y","known_weight")],
                    sim2[,c("id","time","A1","R","A2","Y","known_weight")]);
  sim_long <- sim_long[order(sim_long$id,
                             sim_long$time),];
  sim_long$post <- 1 * (sim_long$time > 0);
  sim_long$last <- 1 * (sim_long$time > 1);
  ################################################################
  # Apply pseudo replication to long data:
  sim_long$wave <- sim_long$time;  # Wave is the same thing as time, except that
  # the pseudo-replicates will be treated as additional waves.;
  rows_to_replicate_long <- sim_long[which(sim_long$R==1),];
  rows_not_to_replicate_long <- sim_long[which(sim_long$R==0),];
  rows_not_to_replicate_long$total_wave <- rows_not_to_replicate_long$time + 1;
  positive_pseudodata_long <- rows_to_replicate_long;
  positive_pseudodata_long$A2 <- 1;
  positive_pseudodata_long$total_wave <- positive_pseudodata_long$time + 1;
  negative_pseudodata_long <- rows_to_replicate_long;
  negative_pseudodata_long$A2 <- -1;
  negative_pseudodata_long$total_wave <- negative_pseudodata_long$time + 4;
  # We keep the same subject ID to show that we don't really have all those
  # new participants.  So we have to distinguish the new observations somehow,
  # and so we treat them as new waves of data on the same person.   
  data_for_analysis_long <- rbind(positive_pseudodata_long,
                                  negative_pseudodata_long,
                                  rows_not_to_replicate_long);
  data_for_analysis_long <- data_for_analysis_long[order(data_for_analysis_long$id,
                                                         data_for_analysis_long$total_wave),];
  ################################################################
  # Start analysis:
  # Fit one-wave models:
  
  gee_mid_only <- geeglm(formula = Y ~ A1 + A2 + A1:A2,
                          id=id, 
                          weights = known_weight, 
                          data=data_for_analysis_long[which(data_for_analysis_long$time==1),],
                          corstr = "independence"); 
  final_mean_estimates_mid_only <- drop(regimens_1wave %*% gee_mid_only$geese$beta);
  target_contrast_mid_only <-
    final_mean_estimates_mid_only[regimen1] -
    final_mean_estimates_mid_only[regimen2] ;
  variance_of_target_contrast_mid_only <-
    dummy_contrast_coefficients %*% regimens_1wave %*%
    gee_mid_only$geese$vbeta %*%
    t(regimens_1wave) %*% dummy_contrast_coefficients;
  std_err_target_contrast_mid_only <-
    drop(sqrt(variance_of_target_contrast_mid_only));
  
  gee_last_only <- geeglm(formula = Y ~ A1 + A2 + A1:A2,
                       id=id, 
                       weights = known_weight, 
                       data=data_for_analysis_long[which(data_for_analysis_long$time==2),],
                       corstr = "independence"); 
  final_mean_estimates_last_only <- drop(regimens_1wave %*% gee_last_only$geese$beta);
  target_contrast_last_only <-
    final_mean_estimates_last_only[regimen1] -
    final_mean_estimates_last_only[regimen2] ;
  variance_of_target_contrast_last_only <-
    dummy_contrast_coefficients %*% regimens_1wave %*%
    gee_last_only$geese$vbeta %*%
    t(regimens_1wave) %*% dummy_contrast_coefficients;
  std_err_target_contrast_last_only <-
    drop(sqrt(variance_of_target_contrast_last_only));

    
  # Fit three-wave model
  # Test contrast under working independence two-wave GEE;
  regimens_after_pre <- cbind(regimens_1wave[,"int"],
                             1, 
                             1,
                             regimens_1wave[,2:4]);
  gee_formula_pre_post <- Y ~  post + last + post:A1 +post:A2 + post:A1:A2;
  # assumes piecewise constant for pretest (time 0) and posttest (times 1-2);
  # or equivalently piecewise linear for time 0 to 1 and time 1 to 2, with the
  # latter slope constrained to zero;
   gee_ind <- geeglm(formula = gee_formula_pre_post,
                             id=id,
                             weights = known_weight, 
                             data=data_for_analysis_long,
                             corstr = "independence" );
  final_mean_estimates_ind <-
     drop(regimens_after_pre %*% gee_ind$geese$beta);
   target_contrast_ind <-
     final_mean_estimates_ind[regimen1] -
     final_mean_estimates_ind[regimen2] ;
   variance_of_target_contrast_ind <-
     dummy_contrast_coefficients %*% regimens_after_pre %*%
     gee_ind$geese$vbeta %*%
     t(regimens_after_pre) %*% dummy_contrast_coefficients;
   std_err_target_contrast_ind <- 
     drop(sqrt(variance_of_target_contrast_ind));

   
   all_crossprods01 <- NULL;
   all_crossprods12 <- NULL;
   all_crossprods02 <- NULL;
   all_weights <- NULL;
   resids0 <- gee_ind$residuals[which(data_for_analysis_long$time==0)];
   resids1 <- gee_ind$residuals[which(data_for_analysis_long$time==1)];
   resids2 <- gee_ind$residuals[which(data_for_analysis_long$time==2)];
   weights2 <- data_for_analysis_long$known_weight[which(data_for_analysis_long$time==2)];
   std_resids0 <- resids0 / sd(resids0);
   std_resids1 <- resids1 / sd(resids1);
   std_resids2 <- resids2 / sd(resids2);
   n_subjects <- length(unique(gee_ind$id));
   for (this_subject in 1:n_subjects) { 
     all_crossprods01 <- c(all_crossprods01,
                           std_resids0[this_subject]*std_resids1[this_subject]);
     all_crossprods12 <- c(all_crossprods12,
                           std_resids1[this_subject]*std_resids2[this_subject]);
     all_crossprods02 <- c(all_crossprods02,
                           std_resids0[this_subject]*std_resids2[this_subject]);
     all_weights <- c(all_weights, gee_ind$weights[min(this_subject)]);
   }
   
   #  Fit working AR1 GEE:
   #...  first Estimate AR1 correlation parameter  
   rho_hat_ar1 <- pmax(0,mean(c(all_crossprods01,all_crossprods12),
                              w=rep(all_weights,each=2))); 
   block_work_corr_ar1 <- rbind(c(1,rho_hat_ar1,rho_hat_ar1^2),
                                c(rho_hat_ar1,1,rho_hat_ar1),
                                c(rho_hat_ar1^2,rho_hat_ar1,1));
   work_corr_ar1 <- rbind(cbind(block_work_corr_ar1,0*block_work_corr_ar1),
                          cbind(0*block_work_corr_ar1,block_work_corr_ar1));
   work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr_ar1,
                                   id=data_for_analysis_long$id,
                                   waves=data_for_analysis_long$total_wave);
   gee_ar1 <- geeglm(formula = gee_formula_pre_post,
                     id=id,
                     weights = known_weight,
                     data=data_for_analysis_long,
                     corstr = "fixed",
                     zcor=work_corr_as_zcor);
   final_mean_estimates_ar1 <-
     drop(regimens_after_pre %*% gee_ar1$geese$beta);
   target_contrast_ar1 <-
     final_mean_estimates_ar1[regimen1] -
     final_mean_estimates_ar1[regimen2] ;
   variance_of_target_contrast_ar1 <-
     dummy_contrast_coefficients %*% regimens_after_pre %*%
     gee_ar1$geese$vbeta %*%
     t(regimens_after_pre) %*% dummy_contrast_coefficients;
   std_err_target_contrast_ar1 <-
     drop(sqrt(variance_of_target_contrast_ar1));
   
 
   #  Fit working exchangeable GEE:
     #...  first Estimate exchangeable correlation parameter 
   rho_hat_ar1 <- pmax(0,mean(c(all_crossprods01,all_crossprods12),
                              w=rep(all_weights,each=2)));
   rho_hat_exch <- pmax(0,mean(c(all_crossprods01,all_crossprods12,all_crossprods12),
                               w=rep(all_weights,each=3)));
   block_work_corr_exch <- rbind(c(1,rho_hat_exch,rho_hat_exch),
                            c(rho_hat_exch,1,rho_hat_exch),
                            c(rho_hat_exch,rho_hat_exch,1));
   work_corr_exch <- rbind(cbind(block_work_corr_exch,0*block_work_corr_exch),
                      cbind(0*block_work_corr_exch,block_work_corr_exch));
   work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr_exch,
                                   id=data_for_analysis_long$id,
                                   waves=data_for_analysis_long$total_wave);
   gee_exch <- geeglm(formula = gee_formula_pre_post,
                              id=id,
                              weights = known_weight,
                              data=data_for_analysis_long,
                              corstr = "fixed",
                              zcor=work_corr_as_zcor);
   final_mean_estimates_exch <-
     drop(regimens_after_pre %*% gee_exch$geese$beta);
   target_contrast_exch <-
     final_mean_estimates_exch[regimen1] -
     final_mean_estimates_exch[regimen2] ;
   variance_of_target_contrast_exch <-
     dummy_contrast_coefficients %*% regimens_after_pre %*%
     gee_exch$geese$vbeta %*%
     t(regimens_after_pre) %*% dummy_contrast_coefficients;
   std_err_target_contrast_exch <-
     drop(sqrt(variance_of_target_contrast_exch));
   
  answer <- list(
    final_means_mid_only = final_mean_estimates_mid_only,
    final_means_last_only = final_mean_estimates_last_only,
    final_means_ind = final_mean_estimates_ind,
    final_means_ar1 = final_mean_estimates_ar1,
    final_means_exch = final_mean_estimates_exch,
    target_contrast_mid_only = target_contrast_mid_only,
    target_contrast_last_only = target_contrast_last_only,
    target_contrast_ind = target_contrast_ind,
    target_contrast_ar1 = target_contrast_ar1,
    target_contrast_exch = target_contrast_exch,
    std_err_target_contrast_mid_only = std_err_target_contrast_mid_only,
    std_err_target_contrast_last_only = std_err_target_contrast_last_only,
    std_err_target_contrast_ind = std_err_target_contrast_ind,
    std_err_target_contrast_ar1 = std_err_target_contrast_ar1,
    std_err_target_contrast_exch = std_err_target_contrast_exch,
    rho_hat_ar1=rho_hat_ar1,
    rho_hat_exch=rho_hat_exch
  );
  return(answer);
}