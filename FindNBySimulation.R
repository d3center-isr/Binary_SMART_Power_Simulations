find_n_by_simulation <- function(n_to_try,
                                 n_sims_each,
                                 true_cell_conditional_means,
                                 prob_R_given_A1_plus,
                                 prob_R_given_A1_minus,
                                 true_corr_matrix) { 
###################################################
# Define model information;
  n_sample_size_scenarios <- length(n_to_try);
  gee_formula <- y ~  s1 +
  s2 +
  s1:a1 +  
  s2:a1 +
  s2:a2 + 
  s2:a1:a2;
n_params <- 7;
n_regimens <- 4;
regimens_design_matrix_time3 <- rbind( c( 1, 1, 1, +1, +1, +1, +1),
                                       c( 1, 1, 1, +1, +1, -1, -1),
                                       c( 1, 1, 1, -1, -1, +1, -1),
                                       c( 1, 1, 1, -1, -1, -1, +1));
rownames(regimens_design_matrix_time3) <- c("pp", "pm", "mp", "mm");
design_matrices_by_regimen <- array(NA,c(3,n_params,n_regimens));
for (j in 1:n_regimens) {
  design_matrices_by_regimen[1,,j] <- c(1,0,0,0,0,0,0);
  design_matrices_by_regimen[2,,j] <- c(1,1,0,ifelse(j<3,+1,-1),0,0,0);
  design_matrices_by_regimen[3,,j] <- regimens_design_matrix_time3[j,];
}
dimnames(design_matrices_by_regimen) <- list(paste("Time",1:3,sep=""),
                                             colnames(regimens_design_matrix_time3),
                                             rownames(regimens_design_matrix_time3));
linear_weights_for_last_wave <- c(0,0,1);
###################################################
# Get ready for first loop;
power_ind <- rep(NA, n_sample_size_scenarios);
power_ar1 <- rep(NA, n_sample_size_scenarios);
power_exch <- rep(NA, n_sample_size_scenarios);
mean_corr_conditional <- matrix(NA,n_sample_size_scenarios,n_sims_each);
mean_rho_estimated_ar1 <- matrix(NA,n_sample_size_scenarios,n_sims_each);
mean_rho_estimated_exch <- matrix(NA,n_sample_size_scenarios,n_sims_each);
for (this_sample_size_scenario in 1:n_sample_size_scenarios) {
  n_subjects <- n_to_try[this_sample_size_scenario];
  log_odds_contrast_ind <- rep(NA,n_sims_each);
  std_err_log_odds_contrast_ind <- rep(NA,n_sims_each);
  log_odds_contrast_exch <- rep(NA,n_sims_each);
  std_err_log_odds_contrast_exch <- rep(NA,n_sims_each);
  log_odds_contrast_ar1 <- rep(NA,n_sims_each);
  std_err_log_odds_contrast_ar1 <- rep(NA,n_sims_each);
  for (this_sim in 1:n_sims_each) {
    simulation <- simulate_binary_smart_3_waves(n_subjects,
                                                true_cell_conditional_means,
                                                prob_R_given_A1_plus,
                                                prob_R_given_A1_minus,
                                                true_corr_matrix);
    sim_data <- simulation$long_data;
    mean_corr_conditional[this_sample_size_scenario,this_sim] <- simulation$mean_corr_conditional;
    #######################################################################################
    # Start analysis:
    #######################################################################################
    # Apply pseudo replication 
    sim_data$known_weight <- 4*(sim_data$R==0)+2*(sim_data$R==1);   
    sim_data$wave <- sim_data$time;  # Wave is the same thing as time, except that 
    # the pseudo-replicates will be treated as additional waves.;
    rows_to_replicate <- sim_data[which(sim_data$R==1),];
    rows_not_to_replicate <- sim_data[which(sim_data$R==0),];
    positive_pseudodata <- rows_to_replicate;
    positive_pseudodata$a2 <- 1;
    negative_pseudodata <- rows_to_replicate;
    negative_pseudodata$a2 <- -1;
    negative_pseudodata$wave <- negative_pseudodata$wave + 3;  
    # We keep the same subject ID to show that we don't really have all those
    # new participants.  So we have to distinguish the new observations somehow,
    # and so we treat them as new waves of data on the same person.  Although
    # it seems very ad-hoc, this method has been shown to be valid.
    # Create the final analysis dataset including replicates.
    data_for_analysis <- rbind(positive_pseudodata,
                               negative_pseudodata, 
                               rows_not_to_replicate);
    data_for_analysis <- data_for_analysis[order(data_for_analysis$subject_id,
                                                 data_for_analysis$wave),]; 
    
    
    # Fit working independence GEE:
    
    gee_ind <- geeglm(formula = gee_formula,  
                      id=subject_id,   
                      weights = known_weight,
                      data=data_for_analysis,
                      corstr = "independence",
                      family=binomial()); 
    gee_cov_beta_ind <- gee_ind$geese$vbeta;
    # Test contrast under working independence GEE;
    estimated_log_odds_by_time_and_regimen_ind <- cbind(
      design_matrices_by_regimen[,,"pp"]%*%gee_ind$geese$beta,
      design_matrices_by_regimen[,,"pm"]%*%gee_ind$geese$beta,
      design_matrices_by_regimen[,,"mp"]%*%gee_ind$geese$beta,
      design_matrices_by_regimen[,,"mm"]%*%gee_ind$geese$beta);
    colnames(estimated_log_odds_by_time_and_regimen_ind) <- c("pp","pm","mp","mm");
    log_odds_contrast_ind[this_sim] <-
      estimated_log_odds_by_time_and_regimen_ind["Time3",regimen1] - 
      estimated_log_odds_by_time_and_regimen_ind["Time3",regimen2] ;
    log_odds_ratio_var_ind <- 
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave) %*%
      rbind(design_matrices_by_regimen[,,regimen1],
            design_matrices_by_regimen[,,regimen2]) %*%
      gee_cov_beta_ind %*%
      t(rbind(design_matrices_by_regimen[,,regimen1],
              design_matrices_by_regimen[,,regimen2])) %*%
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave); 
    std_err_log_odds_contrast_ind[this_sim] <- sqrt(log_odds_ratio_var_ind);
    
    
    
    # Estimate exchangeable correlation parameter
    raw_residuals <- drop(gee_ind$y - gee_ind$fitted.values);
    residual_variances <- gee_ind$fitted.values * (1-gee_ind$fitted.values);
    scaled_residuals <- raw_residuals / sqrt(residual_variances);
    all_crossprods <- NULL;
    for (this_subject in 1:n_subjects) { # subject index;
      these_rows <- which(data_for_analysis$subject_id==this_subject);
      resids <- scaled_residuals[these_rows];
      all_crossprods <- c(all_crossprods,
                          resids[1]*resids[2],
                          resids[2]*resids[3],
                          resids[1]*resids[3]);
    }
    rho_hat_exch <- pmax(0,mean(all_crossprods));
    mean_rho_estimated_exch[this_sample_size_scenario,this_sim] <- rho_hat_exch;
    # Fit working exchangeable GEE:
    block_work_corr <- diag(as.vector(rep(1-rho_hat_exch,3)))+matrix(rho_hat_exch,3,3);
    work_corr <- rbind(cbind(block_work_corr,0*block_work_corr),
                       cbind(0*block_work_corr,block_work_corr)); 
    work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr, 
                                    id=data_for_analysis$subject_id,  
                                    waves=data_for_analysis$wave); 
    gee_exch <- geeglm(formula = gee_formula,  
                       id=subject_id,  
                       family=binomial,
                       weights = known_weight,
                       data=data_for_analysis,
                       scale.fix=TRUE,
                       corstr = "fixed",
                       zcor=work_corr_as_zcor);
    gee_cov_beta_exch <- gee_exch$geese$vbeta;
    # Test contrast under working exchangeable GEE;
    estimated_log_odds_by_time_and_regimen_exch <- cbind(
      design_matrices_by_regimen[,,"pp"]%*%gee_exch$geese$beta,
      design_matrices_by_regimen[,,"pm"]%*%gee_exch$geese$beta,
      design_matrices_by_regimen[,,"mp"]%*%gee_exch$geese$beta,
      design_matrices_by_regimen[,,"mm"]%*%gee_exch$geese$beta);
    colnames(estimated_log_odds_by_time_and_regimen_exch) <- c("pp","pm","mp","mm");
    log_odds_contrast_exch[this_sim] <-
      estimated_log_odds_by_time_and_regimen_exch["Time3",regimen1] - 
      estimated_log_odds_by_time_and_regimen_exch["Time3",regimen2] ;
    log_odds_ratio_var_exch <- 
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave) %*%
      rbind(design_matrices_by_regimen[,,regimen1],
            design_matrices_by_regimen[,,regimen2]) %*%
      gee_cov_beta_exch %*%
      t(rbind(design_matrices_by_regimen[,,regimen1],
              design_matrices_by_regimen[,,regimen2])) %*%
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave); 
    std_err_log_odds_contrast_exch[this_sim] <- sqrt(log_odds_ratio_var_exch);
    
    
    # Estimate AR-1 correlation parameter;
    raw_residuals <- drop(gee_ind$y - gee_ind$fitted.values);
    residual_variances <- gee_ind$fitted.values * (1-gee_ind$fitted.values);
    scaled_residuals <- raw_residuals / sqrt(residual_variances);
    all_crossprods <- NULL;
    for (this_subject in 1:n_subjects) { # subject index;
      these_rows <- which(data_for_analysis$subject_id==this_subject);
      resids <- scaled_residuals[these_rows];
      all_crossprods <- c(all_crossprods,
                          resids[1]*resids[2],
                          resids[2]*resids[3]);
    }
    rho_hat_ar1 <- mean(all_crossprods); 
    mean_rho_estimated_ar1[this_sample_size_scenario,this_sim] <- rho_hat_ar1;
    # Fit working AR-1 GEE
    block_work_corr <- matrix(0,3,3);
    for (this_row in 1:3) {
      for (this_col in 1:3) {
        block_work_corr[this_row,this_col] <- rho_hat_ar1^abs(this_row-this_col);
      }
    }
    work_corr <- rbind(cbind(block_work_corr,0*block_work_corr),
                       cbind(0*block_work_corr,block_work_corr)); 
    work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr, 
                                    id=data_for_analysis$subject_id,  
                                    waves=data_for_analysis$wave); 
    gee_ar1 <- geeglm(formula = gee_formula,  
                      id=subject_id,  
                      family=binomial,
                      weights = known_weight,
                      data=data_for_analysis,
                      scale.fix=TRUE,
                      corstr = "fixed",
                      zcor=work_corr_as_zcor);
    gee_cov_beta_ar1 <- gee_ar1$geese$vbeta;
    estimated_log_odds_by_time_and_regimen_ar1 <- cbind(
      design_matrices_by_regimen[,,"pp"]%*%gee_ar1$geese$beta,
      design_matrices_by_regimen[,,"pm"]%*%gee_ar1$geese$beta,
      design_matrices_by_regimen[,,"mp"]%*%gee_ar1$geese$beta,
      design_matrices_by_regimen[,,"mm"]%*%gee_ar1$geese$beta);
    colnames(estimated_log_odds_by_time_and_regimen_ar1) <- c("pp","pm","mp","mm");
    log_odds_contrast_ar1[this_sim] <-
      estimated_log_odds_by_time_and_regimen_ar1["Time3",regimen1] - 
      estimated_log_odds_by_time_and_regimen_ar1["Time3",regimen2] ;
    log_odds_ratio_var_ar1 <- 
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave) %*%
      rbind(design_matrices_by_regimen[,,regimen1],
            design_matrices_by_regimen[,,regimen2]) %*%
      gee_cov_beta_ar1 %*%
      t(rbind(design_matrices_by_regimen[,,regimen1],
              design_matrices_by_regimen[,,regimen2])) %*%
      c(linear_weights_for_last_wave,
        -linear_weights_for_last_wave); 
    std_err_log_odds_contrast_ar1[this_sim] <- sqrt(log_odds_ratio_var_ar1);
    
  }
  
  power_ind[this_sample_size_scenario] <- mean(abs(log_odds_contrast_ind)>
                                                 qnorm(1-target_alpha/2)*std_err_log_odds_contrast_ind);
  power_ar1[this_sample_size_scenario] <- mean(abs(log_odds_contrast_ar1)>
                                                 qnorm(1-target_alpha/2)*std_err_log_odds_contrast_ar1);
  power_exch[this_sample_size_scenario] <- mean(abs(log_odds_contrast_exch)>
                                                  qnorm(1-target_alpha/2)*std_err_log_odds_contrast_exch);
}


model_ind <- glm(power_ind~n_to_try,family=quasibinomial(link="probit"));
n_needed_ind <- ceiling(as.numeric((qnorm(target_power)-model_ind$coef[1]) / model_ind$coef[2]));

model_ar1 <- glm(power_ar1~n_to_try,family=quasibinomial(link="probit"));
n_needed_ar1 <- ceiling(as.numeric((qnorm(target_power)-model_ar1$coef[1]) / model_ar1$coef[2]));

model_exch <- glm(power_exch~n_to_try,family=quasibinomial(link="probit"));
n_needed_exch <- ceiling(as.numeric((qnorm(target_power)-model_exch$coef[1]) / model_exch$coef[2]));
answer <- list(n_ind=n_needed_ind,
            n_ar1=n_needed_ar1,
            n_exch=n_needed_exch,
            corr_conditional=mean(mean_corr_conditional,na.rm=TRUE),
            corr_conditional_na=mean(is.na(mean_corr_conditional)),
            rho_ar1=mean(mean_rho_estimated_ar1),
            rho_exch=mean(mean_rho_estimated_exch),
            n_tried = n_to_try,
            power_ind =power_ind,
            power_ar1 = power_ar1,
            power_exch = power_exch);
return(answer);
}
