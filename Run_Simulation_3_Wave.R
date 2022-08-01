##################################################################
# Do initial setup;
##################################################################
rm(list = ls());
start_time <- Sys.time(); 
set.seed(18503); 
source("Simulate_Data_Logistic_3_wave.R");
source("Analyze_Simulated_Data_Logistic.R");
source("Analyze_Simulated_Data_Logistic_3_wave.R");
#debug(analyze_simulated_data_logistic_3_wave);
library(geepack);   
##################################################################
# Define basic settings and true parameters:
##################################################################
regimen1 <- "pm";
regimen2 <- "mm";
# Target contrast is assumed to be marginal contrast of 
# pm versus mm, i.e., (+1,-1) versus (-1,-1)
n_sims <- 10000;
target_alpha <- .05; 
overall_results <- NULL; 
target_power <- .8;
zQ <- qnorm(target_power);
mean_fitted_value_by_method_and_regimen <- list();
power_by_method <- list();
mean_target_contrasts <- list();
average_rho_estimates <- list();
this_scenario <- 0;
for (n_subjects in c(300,500)) { 
for (delayed_effect in c(FALSE,TRUE)) { 
  this_scenario <- this_scenario +1; 
  E_Y0 <- .4;
  theta_R <- c(intercept = -.62, 
               y0 = 1,
               a1 = .5);  
  if (delayed_effect==FALSE) {  
    beta_Y2_a1 <- 0;
    beta_Y2_y1 <- 3;
  }
  if (delayed_effect==TRUE) { 
    beta_Y2_a1 <- .275;   # 2.0831 from .3, .5, 1.89 from .25
    beta_Y2_y1 <- .5;  
  }
  beta_Y1 <- c(intercept = -1.55,
               y0 = 3,
               r = 1,
               a1 = .78, # "high" effect size; 
               a2NR = 0,
               a1a2NR = 0);
  beta_Y2 <- c(intercept=-1.4, 
               a1=beta_Y2_a1, 
               y1=beta_Y2_y1);
  #####################################
  n_regimens <- 4;
  n_methods <- 7;
  means <- matrix(NA, n_sims,6);
  cors <- matrix(NA, n_sims,6);
  final_logodds_estimates <- array(NA,c(n_sims,
                                        n_regimens,
                                        n_methods));
  target_contrasts <- matrix(NA,n_sims, n_methods);
  std_err_target_contrasts <- matrix(NA,n_sims, n_methods);
  rho_estimates <- matrix(NA,n_sims, 2); 
  dimnames(final_logodds_estimates)[[2]] <- c("mm","mp","pm","pp");
  dimnames(final_logodds_estimates)[[3]] <- 
    colnames(target_contrasts) <-
    colnames(std_err_target_contrasts) <- c("simple_1",
                                            "exch_01",
                                            "simple_2",
                                            "exch_02",
                                            "ind_012",
                                            "ar1_012",
                                            "exch_012");
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
                                                             beta_Y1,
                                                             beta_Y2);
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
    analysis_results_old_way_1 <- analyze_simulated_data_logistic(two_wave_data_1, 
                                                                  regimen1, 
                                                                  regimen2);
    # pretest and middle only
    two_wave_data_2 <- sim_data_wide;   
    two_wave_data_2$Y1 <- two_wave_data_2$Y2;   
    analysis_results_old_way_2 <- analyze_simulated_data_logistic(two_wave_data_2, 
                                                                  regimen1, 
                                                                  regimen2);
    # pretest and end only
    final_logodds_estimates[this_sim,,] <- cbind(analysis_results$final_means_mid_only,
                                                 analysis_results_old_way_1$final_means_2_waves_exch,
                                                 analysis_results$final_means_last_only,
                                                 analysis_results_old_way_2$final_means_2_waves_exch,
                                                 analysis_results$final_means_ind,
                                                 analysis_results$final_means_ar1,
                                                 analysis_results$final_means_exch);
    target_contrasts[this_sim,] <- c(analysis_results$target_contrast_mid_only,
                                     analysis_results_old_way_1$target_contrast_2_waves_exch,
                                     analysis_results$target_contrast_last_only,
                                     analysis_results_old_way_2$target_contrast_2_waves_exch,
                                     analysis_results$target_contrast_ind,
                                     analysis_results$target_contrast_ar1,
                                     analysis_results$target_contrast_exch);
    std_err_target_contrasts[this_sim,] <- c(analysis_results$std_err_target_contrast_mid_only,
                                             analysis_results_old_way_1$std_err_target_contrast_2_waves_exch,
                                             analysis_results$std_err_target_contrast_last_only,
                                             analysis_results_old_way_2$std_err_target_contrast_2_waves_exch,
                                             analysis_results$std_err_target_contrast_ind,
                                             analysis_results$std_err_target_contrast_ar1,
                                             analysis_results$std_err_target_contrast_exch);
    rho_estimates[this_sim,] <-   c(ar1=analysis_results$rho_hat_ar1,
                                    exch=analysis_results$rho_hat_exch);
  }
  
  mean_fitted_value_by_method_and_regimen[[this_scenario]] <- apply(final_logodds_estimates,c(2,3),mean)
  print(mean_fitted_value_by_method_and_regimen[[this_scenario]]);
  
  power_by_method[[this_scenario]] <- apply(abs(target_contrasts/std_err_target_contrasts)>1.96,2,mean);
  print(power_by_method[[this_scenario]]); 
  
  mean_target_contrasts[[this_scenario]] <- apply(target_contrasts,2,mean);
  print(mean_target_contrasts[[this_scenario]]);
  
  average_rho_estimates[[this_scenario]] <- apply(rho_estimates,2,mean);
  print(average_rho_estimates[[this_scenario]]);
}
}
save.image(paste("Ran-",
                 as.integer(Sys.time()),".rdata",sep="")); 

finish_time <- Sys.time();

print("Final log odds estimates:");
print(apply(final_logodds_estimates,c(2,3),mean));

print("Mean fitted values:");
print(mean_fitted_value_by_method_and_regimen);

print("Time taken:")
print(finish_time-start_time); 

print("Power by method and scenario:")
print(rbind(power_by_method[[1]],
            power_by_method[[2]],
            power_by_method[[3]],
            power_by_method[[4]]));

print("Mean log odds ratios:");
print(round(rbind(mean_target_contrasts[[1]],
                  mean_target_contrasts[[2]],
                  mean_target_contrasts[[3]],
                  mean_target_contrasts[[4]]),4));


print("Mean odds ratios:");
print(round(exp(rbind(mean_target_contrasts[[1]],
                      mean_target_contrasts[[2]],
                      mean_target_contrasts[[3]],
                      mean_target_contrasts[[4]])),4));

save.image("Did_3_Wave_Sim.rdata")