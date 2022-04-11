##################################################################
# Do initial setup;
##################################################################
rm(list = ls());
start_time <- Sys.time(); 
set.seed(16801);
library(bindata);
library(geepack);  
library(dplyr);
source("Simulate_Data.R");  
source("Analyze_Simulated_Data.R");
##################################################################
# Define basic settings:
##################################################################
# Target contrast is assumed to be marginal contrast of 
# pm versus mm, i.e., (+1,-1) versus (-1,-1)
regimen1 <- "pm";
regimen2 <- "mm";
n_sims <- 1000;
target_alpha <- .05; 
link_names <- c("identity","logit");
######################################################
# Loop over sample sizes
######################################################
mean_coef_estimates <- NULL;
power_estimates <- NULL;
rho_hat_summary <- NULL;
scenario <- 0;
for (n_subjects in c(100,186,200)) {  
  scenario <- scenario+1;
  print(paste("Beginning scenario",scenario));
  ######################################################
  # Define true parameters
  ######################################################
  mu_pre <- .50;
  beta_R_0 <- .25;
  beta_R_Ypre <- .40;
  beta_R_A1 <- .10;
  beta_Ypost_0 <- .30;
  beta_Ypost_Ypre <- .20;
  beta_Ypost_R <- 0.25;
  beta_Ypost_A1 <- .10;
  beta_Ypost_A2 <- .05;
  # Intuition for parameters: Y is something bad, R is something good, 
  # A1=+1 is better, A2 does not matter; 
  rho <- .3;  # this is from numerical results; I want a better way
  ######################################################
  ##################################################################
  # Set up data structures to record results;
  ##################################################################
  rho_hats <- rep(NA,n_sims);
  n_regimens <- 4;
  n_vars_wide <- 5;
  n_methods <- 3;
  n_links <- length(link_names);
  coef_estimates <- array(NA,c(n_sims,n_vars_wide,n_methods,n_links));
  target_contrasts <- array(NA,c(n_sims,n_methods,n_links));
  std_err_target_contrasts <- array(NA,c(n_sims,n_methods,n_links));
  dimnames(coef_estimates)[[2]] <- 
    c("intercept","model-specific","time_A1","time_A2","time_A1_A2"); 
  dimnames(coef_estimates)[[3]] <- 
    colnames(target_contrasts) <-
    colnames(std_err_target_contrasts) <- c("1wave","2_ind","2_exch");
  dimnames(coef_estimates)[[4]] <- 
    dimnames(target_contrasts)[[3]] <- 
    dimnames(std_err_target_contrasts)[[3]] <- link_names;
  ##################################################################
  # Estimate power from one-wave formulas based on Kidwell et al (2018);
  ##################################################################
  r1 <- beta_R_0 + beta_R_Ypre*mu_pre + beta_R_A1;
  #P(R=1) for regimen 1 of pairwise comparison, involving A1=+1;
  r2 <- beta_R_0 + beta_R_Ypre*mu_pre - beta_R_A1; 
  #P(R=1) for regimen 2 of pairwise comparison, involving A1=-1;
  pPOST <- c(pp = beta_Ypost_0 + 
               beta_Ypost_Ypre * mu_pre + 
               beta_Ypost_R * r1 + 
               beta_Ypost_A1 +
               beta_Ypost_A2*(1-r1),
             pm = beta_Ypost_0 + 
               beta_Ypost_Ypre * mu_pre + 
               beta_Ypost_R * r1 + 
               beta_Ypost_A1 +
               -beta_Ypost_A2*(1-r1),
             mp = beta_Ypost_0 + 
               beta_Ypost_Ypre * mu_pre + 
               beta_Ypost_R * r2 +  
               -beta_Ypost_A1 +
               beta_Ypost_A2*(1-r2),
             mm = beta_Ypost_0 + 
               beta_Ypost_Ypre * mu_pre + 
               beta_Ypost_R * r2 + 
               -beta_Ypost_A1 +
               -beta_Ypost_A2*(1-r2));
  # pPOST is the expected value (=probability=mean) of the final
  # outcome under each regimen
  vPOST <- pPOST * (1-pPOST); # marginal variance from formula
  # for variance of a random binary variable (p(1-p));
  p1 <- pPOST[regimen1]; #P(Y_post=1) for regimen 1 of pairwise comparison, 
  # which is assumed to be (+1,-1);
  p2 <- pPOST[regimen2]; #P(Y_post=1) for regimen 2 of pairwise comparison,
  # which is assumed to be (-1,-1);
  OR <- abs(as.numeric((p1/(1-p1)) / (p2/(1-p2))));
  v1 <- p1 * (1-p1); # Var(Y_post) for everybody in regimen 1;
  v2 <- p2 * (1-p2); # Var(Y_post) for everybody in regimen 2;
  z_one_minus_half_alpha <- qnorm(1-target_alpha/2);
  var1_margprob_1_wave <- as.numeric( 2*((2-r1)/v1)+2*(2-r2)/(v2));
  # the approximate variance of the log odds ratio, divided by the sample size,
  # from the marginal probabilities formula;
  predicted_power_margprob <- pnorm(sqrt(n_subjects*(log(OR)^2)/
                                           var1_margprob_1_wave)
                                    - z_one_minus_half_alpha ); 
  # this is the version of formula using marginal regimen means as inputs,
  # as Kidwell et al (2018) did; 
  ###################################################
  # Now we will use the version of the formula using conditional means;
  # First calculate Final outcome probabilities for responders, given each regimen:
  pPOSTr <- c(pp = beta_Ypost_0 + 
                beta_Ypost_Ypre * mu_pre + 
                beta_Ypost_R + 
                beta_Ypost_A1,
              pm = beta_Ypost_0 + 
                beta_Ypost_Ypre * mu_pre + 
                beta_Ypost_R + 
                beta_Ypost_A1, # same as pp because A2 doesn't matter for responders;
              mp = beta_Ypost_0 + 
                beta_Ypost_Ypre * mu_pre + 
                beta_Ypost_R -
                beta_Ypost_A1,
              mm = beta_Ypost_0 + 
                beta_Ypost_Ypre * mu_pre + 
                beta_Ypost_R -
                beta_Ypost_A1 ); # same as pp because A2 doesn't matter for responders;
  # pPOSTr is the marginal final outcome probability 
  # for responders in each regimen.;
  vPOSTr <- pPOSTr * (1-pPOSTr);  
  
  
  p1r <- pPOSTr[regimen1];   #P(Y_post=1) for responders in regimen 1, i.e., (+1,-1);
  p2r <- pPOSTr[regimen2]; #P(Y_post=1) for responders in regimen 2, i.e., (-1,-1);
  v1r <- p1r * (1-p1r); # Var(Y_post) for responders in regimen 1;
  v2r <- p2r * (1-p2r); # Var(Y_post) for responders in regimen 2;
  p1nr <- beta_Ypost_0 + 
    beta_Ypost_Ypre * mu_pre + 
    beta_Ypost_A1 - 
    beta_Ypost_A2 ;  #P(Y_post=1) for nonresponders in regimen 1, i.e., (+1,-1);
  p2nr <- beta_Ypost_0 + 
    beta_Ypost_Ypre * mu_pre - 
    beta_Ypost_A1 - 
    beta_Ypost_A2 ;  #P(Y_post=1) for nonresponders in regimen 2, i.e., (-1,-1);
  v1nr <- p1nr * (1-p1nr); # Var(Y_post) for nonresponders in regimen 1;
  v2nr <- p2nr * (1-p2nr); # Var(Y_post) for nonresponders in regimen 2;
  # Compute conditional probabilities version of Kidwell formula;
  v10 <- v1nr + (p1nr-p1)^2;
  v20 <- v2nr + (p2nr-p2)^2;
  v11 <- v1r + (p1r-p1)^2;
  v21 <- v2r + (p2r-p2)^2; 
  var1_condprob_1wave <- (4*(1-r1)*v10 + 2*r1*v11)/(v1^2) + 
    (4*(1-r2)*v20 + 2*r2*v21)/(v2^2); 
  # the approximate variance of the log odds ratio, divided by the sample size,
  # from the conditional probabilities formula;
  predicted_power_condprob <- pnorm(sqrt(n_subjects*(log(OR)^2)/
                                           var1_condprob_1wave)
                                    - z_one_minus_half_alpha ); 
  # this is the modified version of the formula using
  # conditional cell means as inputs; 
  ##################################################################
  # Estimate power from new two-wave formulas ;
  ##################################################################
  # Marginal variances for pretest and for posttest under each regimen:
  vPRE <- mu_pre * (1-mu_pre); # Marginal pretest variance
  temp <- rbind(c(0,-rho*sqrt(vPRE*vPOST)),
                0,
                0,
                0,
                0);
  arrowhead_bread_matrix <- temp + t(temp) + 
    diag(c(4*vPRE,vPOST));
  inv_bread_matrix_two_waves <- (1-rho^2)*solve(arrowhead_bread_matrix);
  prob_R <- (r1+r2)/2;
  dummy_contrast_coefficients <- c(pp=0,pm=1,mp=0,mm=-1);
  # assuming regimen 1 is (+,-) and regimen 2 is (-,-)
  var1_contrast_approx_2_wave <- (4-2*prob_R) * c(0,dummy_contrast_coefficients) %*%
    inv_bread_matrix_two_waves %*% c(0,dummy_contrast_coefficients) ;
  predicted_power_margprob_2_wave <- 
    pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_approx_2_wave) -
            z_one_minus_half_alpha);
  ###################################################
  # Two-wave power formula using conditional probabilities: 
  pPOSTnr <- c(pp = beta_Ypost_0 + 
                 beta_Ypost_Ypre * mu_pre + 
                 beta_Ypost_A1 +
                 beta_Ypost_A2,
               pm = beta_Ypost_0 + 
                 beta_Ypost_Ypre * mu_pre + 
                 beta_Ypost_A1 +
                 -beta_Ypost_A2,
               mp = beta_Ypost_0 + 
                 beta_Ypost_Ypre * mu_pre + 
                 beta_Ypost_A1 +
                 -beta_Ypost_A2,
               mm = beta_Ypost_0 + 
                 beta_Ypost_Ypre * mu_pre + 
                 -beta_Ypost_A1 +
                 -beta_Ypost_A2);
  vPOSTnr <- pPOSTnr * (1-pPOSTnr);  
  Rmat <- rbind(c(1,rho),c(rho,1));
  Rinv <- solve(Rmat); 
  arrowhead_meat_matrix <- matrix(0,5,5);
  regimens_response_rates <- c(r1,r1,r2,r2); 
        # this is P(R=1) given each regimen, noting A2 does not matter;
  for (d in 1:4) {
    Gd <- diag(c(vPRE,vPOST[d]))
    Gd0 <- diag(c(vPRE,vPOSTnr[d]));
    Gd1 <- diag(c(vPRE,vPOSTr[d]));
    rd <- regimens_response_rates[d];
    Qd <- sqrt(Gd) %*% Rinv %*% solve(sqrt(Gd)) %*% 
      (4*(1-rd)*sqrt(Gd0)%*%Rmat%*%sqrt(Gd0) +
         2*rd*sqrt(Gd1)%*%Rmat%*%sqrt(Gd1)) %*% solve(sqrt(Gd)) %*% Rinv %*%sqrt(Gd) ;
    arrowhead_meat_matrix[1,1] <-  arrowhead_meat_matrix[1,1] + Qd[1,1];
    arrowhead_meat_matrix[d+1,1] <- Qd[1,2];
    arrowhead_meat_matrix[1,d+1] <- Qd[1,2];
    arrowhead_meat_matrix[d+1,d+1] <- Qd[2,2];
  }
  var1_contrast_sharp_2_wave <-   c(0,dummy_contrast_coefficients) %*%
    inv_bread_matrix_two_waves %*% arrowhead_meat_matrix  %*%
    inv_bread_matrix_two_waves %*% c(0,dummy_contrast_coefficients) ;
  predicted_power_condprob_2_wave <-  
    pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_sharp_2_wave) -
            z_one_minus_half_alpha);  
  
  
  
  
  
  
  
  ##################################################################
  # Begin loop over simulations
  ##################################################################
  for (this_sim in 1:n_sims) {
    ##################################################################
    # Simulate data;
    ################################################################## 
    sim_wide <- simulate_data(n_subjects,
                              mu_pre,
                              beta_R_0,
                              beta_R_Ypre,
                              beta_R_A1,
                              beta_Ypost_0,
                              beta_Ypost_A1,
                              beta_Ypost_Ypre,
                              beta_Ypost_R,
                              beta_Ypost_A2); 
    ##################################################################
    # Analyze;
    ##################################################################
    for (link_index in 1:n_links) {
      analysis_results <- analyze_simulated_data(sim_wide, 
                                                 regimen1, 
                                                 regimen2,
                                                 which_link=link_names[link_index]);
      ##################################################################
      # Record results;
      ##################################################################
      coef_estimate_1_wave_padded <- c(analysis_results$coef_estimates_1_wave[1],
                                       0,
                                       analysis_results$coef_estimates_1_wave[2:4]);
      coef_estimates[this_sim,,,link_index] <- cbind(coef_estimate_1_wave_padded,
                                                     analysis_results$coef_estimates_2_waves_ind,
                                                     analysis_results$coef_estimates_2_waves_exch);
      target_contrasts[this_sim,,link_index] <- c(analysis_results$target_contrast_1_wave,
                                                  analysis_results$target_contrast_2_waves_ind,
                                                  analysis_results$target_contrast_2_waves_exch);
      std_err_target_contrasts[this_sim,,link_index] <- c(analysis_results$std_err_target_contrast_1_wave,
                                                          analysis_results$std_err_target_contrast_2_waves_ind,
                                                          analysis_results$std_err_target_contrast_2_waves_exch);
      rho_hats[this_sim] <- analysis_results$rho_hat;
    }
  }
  ##################################################################
  # Report summary of results;
  ##################################################################
  finish_time <- Sys.time();
  mean_coef_estimates <- rbind(mean_coef_estimates,
                               apply(coef_estimates,c(2,3,4),mean));
  power_estimates <- rbind(power_estimates,
                           c(n=n_subjects,
                                   pred_marg1 = as.numeric(predicted_power_margprob), 
                                   pred_cond1 = as.numeric(predicted_power_condprob),
                                   pred_marg2 = as.numeric(predicted_power_margprob_2_wave),
                                   pred_cond2 = as.numeric(predicted_power_condprob_2_wave),
                                   line = apply(abs(target_contrasts[,,1])>1.96*std_err_target_contrasts[,,1],2,mean),
                                   lgst = apply(abs(target_contrasts[,,2])>1.96*std_err_target_contrasts[,,2],2,mean)));
  rho_hat_summary <- rbind(rho_hat_summary, summary(rho_hats));
}
save.image(paste("Ran-method-with-probs-",as.integer(Sys.time()),".rdata",sep=""));
print("Simulation duration:");
print(finish_time - start_time);
print("Compare methods:");
  print(round(power_estimates[,c("n",
                                 "line.1wave",
                                 "line.2_ind",
                                 "lgst.1wave",
                                 "lgst.2_ind",
                                 "line.2_exch",
                                 "lgst.2_exch")],4));
print("Compare simulated to expected power:");
  print(round(power_estimates[,c("n",
                               "pred_marg1",
                               "pred_cond1",
                               "lgst.1wave",
                               "pred_marg2",
                               "pred_cond2",
                               "lgst.2_exch")],4));
print("Correlation coefficient estimates:")
print(rho_hat_summary);
