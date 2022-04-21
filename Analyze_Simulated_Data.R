analyze_simulated_data <- function(sim_data_wide, 
                                   regimen1,
                                   regimen2,
                                   which_link) { 
  # Define regimens for estimands:
  regimens_1_wave_model <- rbind(pp=c(1,1,1,1),
                                 pm=c(1,1,-1,-1),
                                 mp=c(1,-1,1,-1),
                                 mm=c(1,-1,-1,1));
  colnames(regimens_1_wave_model) <- c("int","a1","a2","a1a2");
  regimens_2_wave_model <- rbind(pp=c(1,1,1,1,1),
                                 pm=c(1,1,1,-1,-1),
                                 mp=c(1,1,-1,1,-1),
                                 mm=c(1,1,-1,-1,1));
  colnames(regimens_2_wave_model) <- c("int","t","ta1","ta2","ta1a2"); 
  # Specify which contrast to focus on; 
  contrast_coefficients_1_wave <- regimens_1_wave_model[regimen1,]-
    regimens_1_wave_model[regimen2,];
  contrast_coefficients_2_wave <- regimens_2_wave_model[regimen1,]-
    regimens_2_wave_model[regimen2,];
  dummy_contrast_coefficients <- c(pp=0,pm=1,mp=0,mm=-1);
  sim_data_wide$known_weight <- 4*(sim_data_wide$R==0)+2*(sim_data_wide$R==1);
  sim_pre <- sim_data_wide; sim_pre$time <- 0; sim_pre$Y <- sim_pre$Ypre;
  sim_post <- sim_data_wide; sim_post$time <- 1; sim_post$Y <- sim_post$Ypost;
  sim_long <- rbind(sim_pre[,c("id","time","A1","R","A2","Y","known_weight")],
                    sim_post[,c("id","time","A1","R","A2","Y","known_weight")]);
  sim_long <- sim_long[order(sim_long$id,sim_long$time),];
  ################################################################
  # Apply pseudo replication to wide data:
  rows_to_replicate_wide <- sim_data_wide[which(sim_data_wide$R==1),];
  rows_not_to_replicate_wide <- sim_data_wide[which(sim_data_wide$R==0),];
  rows_not_to_replicate_wide$pseudo_wave <- 1;
  positive_pseudodata_wide <- rows_to_replicate_wide;
  positive_pseudodata_wide$A2 <- 1;
  positive_pseudodata_wide$pseudo_wave <- 1;
  negative_pseudodata_wide <- rows_to_replicate_wide;
  negative_pseudodata_wide$A2 <- -1;
  negative_pseudodata_wide$pseudo_wave <- 2;
  # We keep the same subject ID to show that we don't really have all those
  # new participants.  So we have to distinguish the new observations somehow,
  # and so we treat them as new waves of data on the same person.   
  data_for_analysis_wide <- rbind(positive_pseudodata_wide,
                                  negative_pseudodata_wide,
                                  rows_not_to_replicate_wide);
  data_for_analysis_wide <- data_for_analysis_wide[order(data_for_analysis_wide$id,
                                                         data_for_analysis_wide$pseudo_wave),];
  mean_ypre <- weighted.mean(data_for_analysis_wide$Ypre,
                              w=data_for_analysis_wide$known_weight);
  mean_ypost <- weighted.mean(data_for_analysis_wide$Ypost,
                              w=data_for_analysis_wide$known_weight);
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
  negative_pseudodata_long$total_wave <- negative_pseudodata_long$time + 3;
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
  # Fit one-wave model:
  gee_formula_1 <- Ypost ~ A1 + A2 + A1:A2; 
  mean_ypre <- weighted.mean(data_for_analysis_wide$Ypre,
                              w=data_for_analysis_wide$known_weight);
  mean_ypost <- weighted.mean(data_for_analysis_wide$Ypost,
                              w=data_for_analysis_wide$known_weight);
  mustart_wide <- rep(mean_ypost,nrow(data_for_analysis_wide));
  gee_1_wave <- geeglm(formula = gee_formula_1,
                       id=id,
                       family=binomial(link=which_link),
                       weights = known_weight, 
                       data=data_for_analysis_wide,
                      mustart=mustart_wide,
                       corstr = "independence");
  estimated_log_odds_by_time_and_regimen_1_wave <-
    drop(regimens_1_wave_model %*% gee_1_wave$geese$beta);
  log_odds_estimates_1_wave <- estimated_log_odds_by_time_and_regimen_1_wave;
  target_contrast_1_wave <-
    estimated_log_odds_by_time_and_regimen_1_wave[regimen1] -
    estimated_log_odds_by_time_and_regimen_1_wave[regimen2] ;
  variance_of_target_contrast_1_wave <-
    dummy_contrast_coefficients %*% regimens_1_wave_model %*%
    gee_1_wave$geese$vbeta %*%
    t(regimens_1_wave_model) %*% dummy_contrast_coefficients;
  std_err_target_contrast_1_wave <-
    drop(sqrt(variance_of_target_contrast_1_wave));
  #   # Test contrast under working independence two-wave GEE;
  gee_formula_2 <- Y ~  time + time:A1 +time:A2 + time:A1:A2;
  mustart_long <- rep(NA,nrow(data_for_analysis_long));  
  mustart_long[which(data_for_analysis_long$time==0)] <- mean_ypre;
  mustart_long[which(data_for_analysis_long$time==1)] <- mean_ypost;
  gee_2_waves_ind <- geeglm(formula = gee_formula_2,
                            id=id,
                            weights = known_weight,
                            family=binomial(link=which_link),
                            mustart=mustart_long,
                            data=data_for_analysis_long,
                            corstr = "independence" );
  estimated_log_odds_by_time_and_regimen_2_waves_ind <-
    drop(regimens_2_wave_model %*% gee_2_waves_ind$geese$beta);
  log_odds_estimates_2_waves_ind <-
    estimated_log_odds_by_time_and_regimen_2_waves_ind;
  target_contrast_2_waves_ind <-
    estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen1] -
    estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen2] ;
  variance_of_target_contrast_2_waves_ind <-
    dummy_contrast_coefficients %*% regimens_2_wave_model %*%
    gee_2_waves_ind$geese$vbeta %*%
    t(regimens_2_wave_model) %*% dummy_contrast_coefficients;
  # Estimate exchangeable correlation parameter
  raw_residuals <- drop(gee_2_waves_ind$y - gee_2_waves_ind$fitted.values);
  residual_variances <- gee_2_waves_ind$fitted.values *
    (1-gee_2_waves_ind$fitted.values);
  scaled_residuals <- raw_residuals / sqrt(residual_variances);
  all_crossprods <- NULL;
  all_weights <- NULL;
  n_subjects <- length(unique(gee_2_waves_ind$id));
  for (this_subject in 1:n_subjects) { # subject index;
    these_rows <- which(gee_2_waves_ind$id==this_subject);
    resids <- scaled_residuals[these_rows];
    all_crossprods <- c(all_crossprods,
                        resids[1]*resids[2]);
    all_weights <- c(all_weights, gee_2_waves_ind$weights[these_rows[1]]);
  }
  rho_hat <- pmax(0,mean(all_crossprods,w=all_weights));
  #  Fit working exchangeable GEE:
  block_work_corr <- rbind(c(1,rho_hat),c(rho_hat,1));
  work_corr <- rbind(cbind(block_work_corr,0*block_work_corr),
                     cbind(0*block_work_corr,block_work_corr));
  work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr,
                                  id=data_for_analysis_long$id,
                                  waves=data_for_analysis_long$total_wave);
  gee_2_waves_exch <- geeglm(formula = gee_formula_2,
                             id=id,
                             weights = known_weight,
                             data=data_for_analysis_long,
                             family=binomial(link=which_link),
                             mustart=mustart_long,
                             corstr = "fixed",
                             zcor=work_corr_as_zcor);
  estimated_log_odds_by_time_and_regimen_2_waves_exch <-
    drop(regimens_2_wave_model %*% gee_2_waves_exch$geese$beta);
  variance_of_target_contrast_2_waves_exch <-
    dummy_contrast_coefficients %*% regimens_2_wave_model %*%
    gee_2_waves_exch$geese$vbeta %*%
    t(regimens_2_wave_model) %*% dummy_contrast_coefficients;
  
  answer <- list(
    coef_estimates_1_wave = gee_1_wave$geese$beta,
    coef_estimates_2_waves_ind = gee_2_waves_ind$geese$beta,
    coef_estimates_2_waves_exch = gee_2_waves_exch$geese$beta,
    log_odds_estimates_1_wave = estimated_log_odds_by_time_and_regimen_1_wave,
    log_odds_estimates_2_waves_ind =
      estimated_log_odds_by_time_and_regimen_2_waves_ind,
    log_odds_estimates_2_waves_exch =
      estimated_log_odds_by_time_and_regimen_2_waves_exch,
    marginal_prob_estimates_1_wave = plogis(estimated_log_odds_by_time_and_regimen_1_wave),
    marginal_prob_estimates_2_waves_ind =
      plogis(estimated_log_odds_by_time_and_regimen_2_waves_ind),
    marginal_prob_estimates_2_waves_exch =
      plogis(estimated_log_odds_by_time_and_regimen_2_waves_exch),
    target_contrast_1_wave =
      estimated_log_odds_by_time_and_regimen_1_wave[regimen1] -
      estimated_log_odds_by_time_and_regimen_1_wave[regimen2] ,
    target_contrast_2_waves_ind =
      estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen1] -
      estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen2] ,
    target_contrast_2_waves_exch =
      estimated_log_odds_by_time_and_regimen_2_waves_exch[regimen1] -
      estimated_log_odds_by_time_and_regimen_2_waves_exch[regimen2] ,
    std_err_target_contrast_1_wave =
      drop(sqrt(variance_of_target_contrast_1_wave)),
    std_err_target_contrast_2_waves_ind =
      drop(sqrt(variance_of_target_contrast_2_waves_ind)),
    std_err_target_contrast_2_waves_exch =
      drop(sqrt(variance_of_target_contrast_2_waves_exch)),
    rho_hat = rho_hat);
  return(answer);
}
