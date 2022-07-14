##################################################################
# Do initial setup;
##################################################################
rm(list = ls());
start_time <- Sys.time(); 
set.seed(16801); 
source("Simulate_Data_Logistic_3_wave.R");
source("Analyze_Simulated_Data_Logistic.R");
source("Analyze_Simulated_Data_Logistic_3_wave.R");
library(geepack);   
##################################################################
# Define basic settings and true parameters:
##################################################################
regimen1 <- "pm";
regimen2 <- "mm";
# Target contrast is assumed to be marginal contrast of 
# pm versus mm, i.e., (+1,-1) versus (-1,-1)
n_sims <- 5000;
target_alpha <- .05; 
overall_results <- NULL; 
target_power <- .8;
zQ <- qnorm(target_power);
mean_by_method_and_regimen <- list();
power_by_method <- list();
mean_target_contrasts <- list();
average_rho_estimates <- list();
this_scenario <- 0;

for (effect_size in c("low","medium","high")) {
  for (n_subjects in c(300,500)) {
    this_scenario <- this_scenario +1;
    beta_Y2_y0 <- 3;
    beta_Y2_a1 <- switch(effect_size,
                         "low"=.22, 
                         "medium"=.45,
                         "high"=.78);  
    beta_Y2_intercept <- -1.55;
    E_Y0 <- .4;
    theta_R <- c(intercept = -.62, 
                 y0 = 1,
                 a1 = .5);
    beta_Y1 <- c(intercept = -1.55,
                 y0 = 3,
                 r = 1,
                 a1 = beta_Y2_a1,
                 a2NR = 0,
                 a1a2NR = 0);
    #####################################
    n_regimens <- 4;
    n_methods <- 7;
    means <- matrix(NA, n_sims,6);
    cors <- matrix(NA, n_sims,6);
    final_mean_estimates <- array(NA,c(n_sims,
                                       n_regimens,
                                       n_methods));
    target_contrasts <- matrix(NA,n_sims, n_methods);
    std_err_target_contrasts <- matrix(NA,n_sims, n_methods);
    rho_estimates <- matrix(NA,n_sims, 2); 
    dimnames(final_mean_estimates)[[2]] <- c("mm","mp","pm","pp");
    dimnames(final_mean_estimates)[[3]] <- 
      colnames(target_contrasts) <-
      colnames(std_err_target_contrasts) <- c("mid",
                                              "last",
                                              "alt_mid",
                                              "alt_last",
                                              "ind",
                                              "ar1",
                                              "exch");
    ##################################################################
    # Begin loop over simulations
    ##################################################################
    for (this_sim in 1:n_sims) { 
      ##################################################################
      # Simulate data;
      ################################################################## 
      sim_data_wide <- simulate_data_logistic_3_wave(n_subjects, 
                                                     E_Y0,
                                                     theta_R,
                                                     beta_Y1);
      means[this_sim,] <- apply(sim_data_wide[,-1],2,mean);
      temp <- cor(sim_data_wide[,c("Y0","R","Y1","Y2")])
      cors[this_sim,] <- as.vector(temp[upper.tri(temp)]);
      colnames(cors) <- c("Y0_R",
                          "Y0_Y1",
                          "Y0_Y2",
                          "R_Y1",
                          "R_Y2",
                          "Y1_Y2");
      analysis_results <- analyze_simulated_data_logistic_3_wave(sim_data_wide, 
                                                                 regimen1, 
                                                                 regimen2);
      two_wave_data_1 <- sim_data_wide;  
      two_wave_data_1$Y2 <- two_wave_data_1$Y1; two_wave_data_1$Y1 <- NA;
      # Note that in the original version of the two-wave code, the outcomes
      # were called Y0 and Y2 in order to emphasize that Y2 occurs after A2, not
      # before. We later changed the notation on the two-wave case to Y0 and Y1
      # anyway. So in order to reuse old code,whichever wave is being treated
      # as posttest in this simulation experiment is temporarily renamed as Y2.
      two_wave_data_2 <- sim_data_wide;  
      two_wave_data_2$Y1 <- NA;
      analysis_results_old_way_1 <- analyze_simulated_data_logistic(two_wave_data_1, 
                                                                    regimen1, 
                                                                    regimen2);
      # pretest and middle only
      analysis_results_old_way_2 <- analyze_simulated_data_logistic(two_wave_data_2, 
                                                                    regimen1, 
                                                                    regimen2);
      # pretest and end only
      final_mean_estimates[this_sim,,] <- cbind(analysis_results$final_means_mid_only,
                                                analysis_results$final_means_last_only,
                                                analysis_results_old_way_1$final_means_2_waves_exch,
                                                analysis_results_old_way_2$final_means_2_waves_exch,
                                                analysis_results$final_means_ind,
                                                analysis_results$final_means_ar1,
                                                analysis_results$final_means_exch);
      target_contrasts[this_sim,] <- c(analysis_results$target_contrast_mid_only,
                                       analysis_results$target_contrast_last_only,
                                       analysis_results_old_way_1$target_contrast_2_waves_exch,
                                       analysis_results_old_way_2$target_contrast_2_waves_exch,
                                       analysis_results$target_contrast_ind,
                                       analysis_results$target_contrast_ar1,
                                       analysis_results$target_contrast_exch);
      std_err_target_contrasts[this_sim,] <- c(analysis_results$std_err_target_contrast_mid_only,
                                               analysis_results$std_err_target_contrast_last_only,
                                               analysis_results_old_way_1$std_err_target_contrast_2_waves_exch,
                                               analysis_results_old_way_2$std_err_target_contrast_2_waves_exch,
                                               analysis_results$std_err_target_contrast_ind,
                                               analysis_results$std_err_target_contrast_ar1,
                                               analysis_results$std_err_target_contrast_exch);
      rho_estimates[this_sim,] <-   c(ar1=analysis_results$rho_hat_ar1,
                                      exch=analysis_results$rho_hat_exch);
    }
    
    mean_by_method_and_regimen[[this_scenario]] <- apply(final_mean_estimates,c(2,3),mean)
    print(mean_by_method_and_regimen[[this_scenario]]);
    
    power_by_method[[this_scenario]] <- apply(abs(target_contrasts/std_err_target_contrasts)>1.96,2,mean);
    print(power_by_method[[this_scenario]]); 
    
    mean_target_contrasts[[this_scenario]] <- apply(target_contrasts,2,mean);
    print(mean_target_contrasts[[this_scenario]]);
    
    average_rho_estimates[[this_scenario]] <- apply(rho_estimates,2,mean);
    print(average_rho_estimates[[this_scenario]]);
  }
  save.image(paste("Ran-",
                   as.integer(Sys.time()),".rdata",sep=""));
}
finish_time <- Sys.time();
print(finish_time-start_time); 
print("Power by method and scenario:")
print(rbind(power_by_method[[1]],
            power_by_method[[2]],
            power_by_method[[3]],
            power_by_method[[4]],
            power_by_method[[5]]))
save.image("Did_3_Wave_Sim.rdata");