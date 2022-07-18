##################################################################
# Do initial setup;
##################################################################
rm(list = ls());
start_time <- Sys.time(); 
set.seed(1855); 
source("Simulate_Data_Logistic.R");
source("Analyze_Simulated_Data_Logistic.R");
library(geepack);   
##################################################################
# Define basic settings and true parameters:
##################################################################
regimen1 <- "pm";
regimen2 <- "mm";
# Target contrast is assumed to be marginal contrast of 
# pm versus mm, i.e., (+1,-1) versus (-1,-1)
n_sims <- 1000;
target_alpha <- .05; 
target_power <- .8;
previous_results <- read.csv(file="results-step-1.csv");
results_all_conditions_and_sizes <- NULL; 
results_calculated_targets <- NULL;
theta_R_a1  <- .5;
beta_Y2_r <- 1;
theta_R_y0 <- 1;
for (condition in 1:nrow(previous_results)) {
  if (previous_results$n_subjects[condition]==500) {  # We are recalculating n_subjects
    # based on effect size so we do not
    # need every n_subjects value from
    # the previous analysis.
    min_n_to_try <- floor(.8 * previous_results$n_from_cond_2[condition]);
    max_n_to_try <- ceiling(1.2 * previous_results$n_from_marg_1[condition]);
    pre_post_corr <- previous_results$pre_post_corr[condition];
    effect_size <- previous_results$effect_size[condition];
    results_this_condition <- NULL;
    for (n_subjects  in round(seq(min_n_to_try, max_n_to_try,length=10))) {
      beta_Y2_y0 <- switch(pre_post_corr,
                           "no"=0,
                           "low"=1.2,
                           "high"=3);
      if (pre_post_corr=="no") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.10,
                             "medium"=.25,
                             "high"=.47);
        beta_Y2_intercept <- -0.44;  
      }
      if (pre_post_corr=="low") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.115, 
                             "medium"=.29,  
                             "high"=.52);  
        beta_Y2_intercept <- -0.9;   
      }
      if (pre_post_corr=="high") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.22, 
                             "medium"=.45,
                             "high"=.78);  
        beta_Y2_intercept <- -1.55;
      }
      E_Y0 <- .4;
      theta_R <- c(intercept = -.62, 
                   y0 = theta_R_y0,
                   a1 = theta_R_a1);
      beta_Y2 <- c(intercept = beta_Y2_intercept,
                   y0 = beta_Y2_y0,
                   r = beta_Y2_r,
                   a1 = beta_Y2_a1,
                   a2NR = 0,
                   a1a2NR = 0);
      ##################################################################
      # Set up data structures to record simulation results;
      ##################################################################
      n_regimens <- 4;
      n_methods <- 3;
      target_contrasts <- matrix(NA,n_sims, n_methods);
      std_err_target_contrasts <- matrix(NA,n_sims, n_methods);
      colnames(target_contrasts) <-
        colnames(std_err_target_contrasts) <- c("1wave","2_ind","2_exch");
      ######################################################
      print(paste(theta_R_y0,beta_Y2_y0,beta_Y2_r));
      ##################################################################
      # Begin loop over simulations
      ##################################################################
      for (this_sim in 1:n_sims) { 
        ##################################################################
        # Simulate data;
        ################################################################## 
        sim_wide <- simulate_data_logistic(n_subjects, 
                                           E_Y0,
                                           theta_R,
                                           beta_Y2); 
        ##################################################################
        # Perform analyses;
        ##################################################################
        analysis_results <- analyze_simulated_data_logistic(sim_wide, 
                                                            regimen1, 
                                                            regimen2);
        ##################################################################
        # Record results;
        ##################################################################
        target_contrasts[this_sim,] <- 
          c(analysis_results$target_contrast_1_wave,
            analysis_results$target_contrast_2_waves_ind,
            analysis_results$target_contrast_2_waves_exch);
        std_err_target_contrasts[this_sim,] <- 
          c(analysis_results$std_err_target_contrast_1_wave,
            analysis_results$std_err_target_contrast_2_waves_ind,
            analysis_results$std_err_target_contrast_2_waves_exch);
      } 
      ##################################################################
      # Report summary of results;
      ##################################################################
      reject_H0 <- abs(target_contrasts / std_err_target_contrasts) > 
        qnorm(1-target_alpha/2); 
      power_estimates <- apply(reject_H0,2,mean);
      print(power_estimates);
      results_this_condition <- rbind(results_this_condition,
                                      c(n_subjects=n_subjects,
                                        power_estimates)); 
      results_all_conditions_and_sizes <- rbind(results_all_conditions_and_sizes,
                                                results_this_condition);
    }
    # Calculate needed sample size according to simulations
    model_1wave <- glm(results_this_condition[,"1wave"]~results_this_condition[,"n_subjects"],
                       family=quasibinomial(link="probit"));
    n_needed_1wave <- ceiling(as.numeric((qnorm(target_power)-model_1wave$coef[1])/
                                           model_1wave$coef[2]));
    model_2wave <- glm(results_this_condition[,"2_exch"]~results_this_condition[,"n_subjects"],
                       family=quasibinomial(link="probit"));
    n_needed_2wave <- ceiling(as.numeric((qnorm(target_power)-model_2wave$coef[1])/
                                           model_2wave$coef[2]));
    # Append to data frame of needed sample sizes  
    results_calculated_targets <- rbind(results_calculated_targets,
                                        data.frame(pre_post_corr=pre_post_corr,
                                                   effect_size=effect_size,
                                                   n_from_marg_1=previous_results$n_from_marg_1[condition],
                                                   n_from_cond_1=previous_results$n_from_cond_1[condition],
                                                   n_needed_1wave=n_needed_1wave,
                                                   n_from_marg_2=previous_results$n_from_marg_2[condition],
                                                   n_from_cond_2=previous_results$n_from_cond_2[condition],
                                                   n_needed_2wave=n_needed_2wave));
  } 
}
print(results_calculated_targets);
print("Simulation duration:");
finish_time <- Sys.time();
print(finish_time - start_time);
save.image(paste("Ran-",
                 as.integer(Sys.time()),".rdata",sep=""));
write.csv(x=results_calculated_targets,file="results-step-2.csv")