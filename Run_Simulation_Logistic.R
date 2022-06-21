##################################################################
# Do initial setup;
##################################################################
rm(list = ls());
start_time <- Sys.time(); 
set.seed(16801); 
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
n_sims <- 200;
target_alpha <- .05; 
overall_results <- NULL; 
theta_R_a1  <- .5;
beta_Y2_r <- 1;
theta_R_y0 <- 1;
beta_Y2_intercept <- -1;
for (pre_post_corr in c("no",
                        "low",
                        "high")) {
  for (effect_size in c("low",
                        "medium",
                        "high")) {
    for (n_subjects  in c(300,500)) {
      beta_Y2_y0 <- switch(pre_post_corr,
                           "no"=0,
                           "low"=1,
                           "high"=2.5);
      if (pre_post_corr=="no") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.10,
                             "medium"=.25,
                             "high"=.46);
      }
      if (pre_post_corr=="low") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.14,
                             "medium"=.33,
                             "high"=.58);
      }
      if (pre_post_corr=="high") {
        beta_Y2_a1 <- switch(effect_size,
                             "low"=.18,
                             "medium"=.40,
                             "high"=.7);
      }
      E_Y0 <- .4;
      theta_R <- c(intercept = -.5,
                   y0 = theta_R_y0,
                   a1 = theta_R_a1);
      beta_Y2 <- c(intercept = beta_Y2_intercept,
                   y0 = beta_Y2_y0,
                   r = beta_Y2_r,
                   a1 = beta_Y2_a1,
                   a2NR = 0,
                   a1a2NR = 0);
      ##################################################################
      # Calculate marginal probabilities for R and Y2:
      ##################################################################
      # Model for generating response status
      # (from pretest and initial treatment option);
      design_R <- expand.grid(a1=c(-1,1),
                              y0=c(0,1))[,2:1];
      cond_E_R <- plogis(cbind(1,as.matrix(design_R)) %*% theta_R);
      rownames(cond_E_R) <- paste(design_R$y0,
                                  design_R$a1);
      marg_E_R <- c(m=NA,
                    p=NA);
      marg_E_R["m"] <- (1-E_Y0)*cond_E_R[which(design_R$y0==0 & design_R$a1==-1)] +
        E_Y0*cond_E_R[which(design_R$y0==1 & design_R$a1==-1)];
      marg_E_R["p"] <- (1-E_Y0)*cond_E_R[which(design_R$y0==0 & design_R$a1==+1)] +
        E_Y0*cond_E_R[which(design_R$y0==1 & design_R$a1==+1)];
      print("Calculated regimen marginal expected values for R:")
      print(marg_E_R); 
      design_Y2 <- expand.grid(a2=c(-1,1),
                               a1=c(-1,1),
                               r=c(0,1),
                               y0=c(0,1),
                               intercept=1)[,5:1];
      design_Y2$a2NR <- (1-design_Y2$r)*design_Y2$a2; # code responders as zero here;
      design_Y2$a1a2NR <- (1-design_Y2$r)*design_Y2$a1*design_Y2$a2; 
      these_columns <- c("intercept","y0","r","a1","a2NR","a1a2NR");
      cond_E_Y2 <- plogis(as.matrix(design_Y2[,these_columns]) %*%
                            beta_Y2[these_columns]);
      print("Calculated conditional expected values for Y2 by path:")
      print(cbind(design_Y2,cond_E_Y2));
      design_Y2_marg <- expand.grid(a1=c(-1,1),
                                    a2NR=c(-1,1));
      regimen_names <- paste( ifelse(design_Y2_marg$a1>0,"p","m"), 
                              ifelse(design_Y2_marg$a2NR>0,"p","m"), 
                              sep="");
      marg_E_Y2 <- rep(0,4);
      names(marg_E_Y2) <- regimen_names;
      for (row in 1:4) {
        a1 <- design_Y2_marg$a1[row];
        a2 <- design_Y2_marg$a2NR[row];
        for (y0 in c(0,1)) {
          for (r in c(0,1)) {
            p_this_y0_value <- ifelse(y0==1, 
                                      E_Y0, 
                                      1-E_Y0);
            this_cond_E_R <- cond_E_R[which(design_R$y0==y0 &
                                              design_R$a1==a1)]; 
            p_this_r_value <- ifelse(r==1, 
                                     this_cond_E_R, 
                                     1-this_cond_E_R);
            this_cond_E_Y2 <- cond_E_Y2[which(design_Y2$y0==y0 &
                                                design_Y2$r==r &  
                                                design_Y2$a1==a1 &
                                                design_Y2$a2==a2)];
            marg_E_Y2[row] <- marg_E_Y2[row] + 
              p_this_y0_value * p_this_r_value * this_cond_E_Y2;
          }
        }
      }
      print("Calculated regimen marginal expected values for Y2:");
      print(marg_E_Y2);
      ##################################################################
      # Set up data structures to record simulation results;
      ##################################################################
      n_regimens <- 4;
      n_methods <- 3;
      final_mean_estimates <- array(NA,c(n_sims,
                                         n_regimens,
                                         n_methods));
      target_contrasts <- matrix(NA,n_sims, n_methods);
      std_err_target_contrasts <- matrix(NA,n_sims, n_methods);
      rho_hat_estimates <- rep(NA, n_sims);
      cor_Y0_R <- rep(NA, n_sims);
      cor_Y0_Y2 <- rep(NA, n_sims);
      cor_R_Y2 <- rep(NA, n_sims);
      mean_R <- rep(NA, n_sims);
      mean_Y2 <- rep(NA, n_sims);
      dimnames(final_mean_estimates)[[2]] <- c("mm","mp","pm","pp");
      dimnames(final_mean_estimates)[[3]] <- 
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
        final_mean_estimates[this_sim,,] <- 
          cbind(analysis_results$final_means_1_wave,
                analysis_results$final_means_2_waves_ind,
                analysis_results$final_means_2_waves_exch); 
        target_contrasts[this_sim,] <- 
          c(analysis_results$target_contrast_1_wave,
            analysis_results$target_contrast_2_waves_ind,
            analysis_results$target_contrast_2_waves_exch);
        std_err_target_contrasts[this_sim,] <- 
          c(analysis_results$std_err_target_contrast_1_wave,
            analysis_results$std_err_target_contrast_2_waves_ind,
            analysis_results$std_err_target_contrast_2_waves_exch);
        rho_hat_estimates[this_sim] <- analysis_results$rho_hat;
        mean_R[this_sim] <- mean(sim_wide$R); 
        mean_Y2[this_sim] <- mean(sim_wide$Y2);
        cor_Y0_R[this_sim] <- cor(sim_wide)["Y0","R"]; 
        cor_Y0_Y2[this_sim] <- cor(sim_wide)["Y0","Y2"];
        cor_R_Y2[this_sim] <- cor(sim_wide)["R","Y2"];
      }
      ##################################################################
      # Use simulation estimate of correlation parameter;
      ##################################################################
      print("Correlation parameter rho:")
      rho <- mean(rho_hat_estimates);  
      ##################################################################
      # Estimate power from one-wave formulas based on Kidwell et al (2018)
      # for each sample size scenario;
      ################################################################## 
      r1 <- marg_E_R[substr(regimen1,1,1)];
      #P(R=1 | A1) given regimen 1 of pairwise comparison;
      r2 <- marg_E_R[substr(regimen2,1,1)]; 
      #P(R=1 | A1) given regimen 2 of pairwise comparison;
      p1 <- marg_E_Y2[regimen1]; #P(Y_post=1) for regimen 1 of pairwise comparison;
      p2 <- marg_E_Y2[regimen2]; #P(Y_post=1) for regimen 2 of pairwise comparison;
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
      names(predicted_power_margprob) <- paste("n=",n_subjects);
      # this is the version of formula using marginal regimen means as inputs,
      # as Kidwell et al (2018) did; 
      print("Predicted power, working independence, from marginal probabilities:");
      print(predicted_power_margprob);
      ###################################################
      # Now we try the version of the power formula using conditional means;
      # First calculate final outcome probabilities for responders, 
      # given each regimen:
      design_for_responders <- design_Y2[which(design_Y2$r==1 & design_Y2$a2==-1),
                                         c("y0","r","a1")];
      # a2 is irrelevant here (because a2NR=0*a2) so we just pick one value
      E_Y2_for_responders <- cond_E_Y2[which(design_Y2$r==1 & design_Y2$a2==-1)];
      marg_E_Y2_for_responders <- c(
        m=(1-E_Y0)*E_Y2_for_responders[which(design_for_responders$y0==0 &
                                               design_for_responders$a1==-1)] + 
          E_Y0 * E_Y2_for_responders[which(design_for_responders$y0==1 &
                                             design_for_responders$a1==-1)],
        p=(1-E_Y0)*E_Y2_for_responders[which(design_for_responders$y0==0 &
                                               design_for_responders$a1==+1)] + 
          E_Y0 * E_Y2_for_responders[which(design_for_responders$y0==1 &
                                             design_for_responders$a1==+1)]) 
      # = E(Y2 | R=1, A1=a1) marginal over Y0;
      pPOSTr <- c(  mm=as.numeric(marg_E_Y2_for_responders["m"]),
                    mp=as.numeric(marg_E_Y2_for_responders["m"]),
                    pm=as.numeric(marg_E_Y2_for_responders["p"]),
                    pp=as.numeric(marg_E_Y2_for_responders["p"]))
      # = E(Y2|A1) for each regimen, which involves duplicating rows 
      # because the second-phase treatment does not matter for responders;
      vPOSTr <- pPOSTr * (1-pPOSTr);  
      p1r <- pPOSTr[regimen1];   #P(Y_post=1) for responders in regimen 1, i.e., (+1,-1);
      p2r <- pPOSTr[regimen2]; #P(Y_post=1) for responders in regimen 2, i.e., (-1,-1);
      v1r <- p1r * (1-p1r); # Var(Y_post) for responders in regimen 1;
      v2r <- p2r * (1-p2r); # Var(Y_post) for responders in regimen 2;
      # Now calculate analogous quantities for nonresponders;
      design_for_nonresponders <- design_Y2[which(design_Y2$r==0),];
      E_Y2_for_nonresponders <- cond_E_Y2[which(design_Y2$r==0)];
      marg_E_Y2_for_nonresponders <- c(
        mm=(1-E_Y0)*E_Y2_for_nonresponders[which(design_for_nonresponders$y0==0 &
                                                   design_for_nonresponders$a1==-1 &
                                                   design_for_nonresponders$a2NR==-1)] + 
          E_Y0 * E_Y2_for_nonresponders[which(design_for_nonresponders$y0==1 &
                                                design_for_nonresponders$a1==-1 &
                                                design_for_nonresponders$a2NR==-1)],
        mp=(1-E_Y0)*E_Y2_for_nonresponders[which(design_for_nonresponders$y0==0 &
                                                   design_for_nonresponders$a1==-1 &
                                                   design_for_nonresponders$a2NR==+1)] + 
          E_Y0 * E_Y2_for_nonresponders[which(design_for_nonresponders$y0==1 &
                                                design_for_nonresponders$a1==-1 &
                                                design_for_nonresponders$a2NR==+1)],
        pm=(1-E_Y0)*E_Y2_for_nonresponders[which(design_for_nonresponders$y0==0 &
                                                   design_for_nonresponders$a1==+1 &
                                                   design_for_nonresponders$a2NR==-1)] + 
          E_Y0 * E_Y2_for_nonresponders[which(design_for_nonresponders$y0==1 &
                                                design_for_nonresponders$a1==+1 &
                                                design_for_nonresponders$a2NR==-1)],
        pp=(1-E_Y0)*E_Y2_for_nonresponders[which(design_for_nonresponders$y0==0 &
                                                   design_for_nonresponders$a1==+1 &
                                                   design_for_nonresponders$a2NR==+1)] + 
          E_Y0 * E_Y2_for_nonresponders[which(design_for_nonresponders$y0==1 &
                                                design_for_nonresponders$a1==+1 &
                                                design_for_nonresponders$a2NR==+1)])
      p1nr <- marg_E_Y2_for_nonresponders[regimen1];
      p2nr <- marg_E_Y2_for_nonresponders[regimen2];
      v1nr <- p1nr * (1-p1nr); # Var(Y_post) for nonresponders in regimen 1;
      v2nr <- p2nr * (1-p2nr); # Var(Y_post) for nonresponders in regimen 2;
      v10 <- v1nr + (p1nr-p1)^2;
      v20 <- v2nr + (p2nr-p2)^2;
      v11 <- v1r + (p1r-p1)^2;
      v21 <- v2r + (p2r-p2)^2; 
      var1_condprob_1wave <- as.numeric((4*(1-r1)*v10 + 2*r1*v11)/(v1^2) + 
                                          (4*(1-r2)*v20 + 2*r2*v21)/(v2^2)); 
      # the approximate variance of the log odds ratio, divided by the sample size,
      # from the conditional probabilities formula;
      predicted_power_condprob <- pnorm(sqrt(n_subjects*(log(OR)^2)/
                                               var1_condprob_1wave)
                                        - z_one_minus_half_alpha ); 
      # this is the modified version of the formula using
      # conditional cell means as inputs; 
      print("Predicted power, working independence, from conditional probabilities:") 
      print(predicted_power_condprob);
      ##################################################################
      # Estimate power from new two-wave formulas ;
      ##################################################################
      # Marginal variances for pretest and for posttest under each regimen:
      vPRE <- E_Y0 * (1-E_Y0); # Marginal pretest variance
      vPOST <- marg_E_Y2 * (1-marg_E_Y2); # marginal variance from formula
      # for variance of a random binary variable (p(1-p));
      temp <- rbind(c(0,-rho*sqrt(vPRE*vPOST)),
                    0,
                    0,
                    0,
                    0);
      arrowhead_bread_matrix <- temp + t(temp) + 
        diag(c(4*vPRE,vPOST));
      inv_bread_matrix_two_waves <- (1-rho^2)*solve(arrowhead_bread_matrix);
      prob_R <- (r1+r2)/2;
      dummy_contrast_coefficients <- 0*marg_E_Y2;
      dummy_contrast_coefficients[regimen1] <- +1;
      dummy_contrast_coefficients[regimen2] <- -1;
      # assuming regimen 1 is (+,-) and regimen 2 is (-,-)
      var1_contrast_approx_2_wave <- as.numeric( (4-2*prob_R) * c(0,dummy_contrast_coefficients) %*%
                                                   inv_bread_matrix_two_waves %*% c(0,dummy_contrast_coefficients)) ;
      predicted_power_margprob_2_wave <- 
        pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_approx_2_wave) -
                z_one_minus_half_alpha); 
      print("Predicted power, correlated, from marginal probabilities:");
      print(predicted_power_margprob_2_wave);
      ###################################################
      # Two-wave power formula using conditional probabilities: 
      vPOSTnr <- marg_E_Y2_for_nonresponders*(1-marg_E_Y2_for_nonresponders);  
      Rmat <- rbind(c(1,rho),c(rho,1));
      Rinv <- solve(Rmat); 
      arrowhead_meat_matrix <- matrix(0,5,5);
      regimens_response_rates <- 0*vPOST;
      regimens_response_rates["mm"] <- marg_E_R["m"];
      regimens_response_rates["mp"] <- marg_E_R["m"];
      regimens_response_rates["pm"] <- marg_E_R["p"];
      regimens_response_rates["pp"] <- marg_E_R["p"];
      # this is P(R=1) given each regimen, noting A2 does not matter 
      # because it is randomized and happens after R;
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
      var1_contrast_sharp_2_wave <- as.numeric( c(0,dummy_contrast_coefficients) %*%
                                                  inv_bread_matrix_two_waves %*% arrowhead_meat_matrix  %*%
                                                  inv_bread_matrix_two_waves %*% c(0,dummy_contrast_coefficients)) ;
      predicted_power_condprob_2_wave <-  
        pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_sharp_2_wave) -
                z_one_minus_half_alpha);  
      names(predicted_power_condprob_2_wave) <- paste("n=",n_subjects);
      print("Predicted power, correlated, from conditional probabilities:");
      print(predicted_power_condprob_2_wave);
      ##################################################################
      # Report summary of results;
      ##################################################################
      reject_H0 <- abs(target_contrasts / std_err_target_contrasts) > 
        qnorm(1-target_alpha/2);
      print("Comparison of power across methods:");
      power_estimates <- apply(reject_H0,2,mean);
      print(power_estimates);
      power_gain <- power_estimates["2_exch"]-power_estimates["1wave"];
      print("Comparison of working independence power formula:");
      print(cbind(n=n_subjects,
                  from_marg=predicted_power_margprob,
                  from_cond=predicted_power_condprob,
                  observed=power_estimates["2_ind"]));
      print("Comparison of working non-independence power formula:");
      print(cbind(n=n_subjects,
                  from_marg=predicted_power_margprob_2_wave,
                  from_cond=predicted_power_condprob_2_wave,
                  observed=power_estimates["2_exch"]));
      print("Correlation coefficient estimates:");
      print(summary(rho_hat_estimates));
      print("Simulation duration:");
      finish_time <- Sys.time();
      print(finish_time - start_time);
      overall_results <- rbind(overall_results,
                               data.frame(OR=OR,
                                          n_subjects=n_subjects,
                                          pre_post_corr=pre_post_corr,
                                          effect_size=effect_size,
                                          beta_Y2_a1=beta_Y2_a1,  
                                          beta_Y2_y0=beta_Y2_y0, 
                                          est_mm=mean(final_mean_estimates[,1,1]),
                                          est_mp=mean(final_mean_estimates[,2,1]),
                                          est_pm=mean(final_mean_estimates[,3,1]),
                                          est_pp=mean(final_mean_estimates[,4,1]), 
                                          mean_R=round(mean(mean_R),2),
                                          mean_Y2=round(mean(mean_Y2),2),
                                          cor_Y0_R=round(mean(cor_Y0_R),2),
                                          cor_Y0_Y2=round(mean(cor_Y0_Y2),2),
                                          cor_R_Y2=round(mean(cor_R_Y2),2),
                                          rho_hat=round(mean(rho_hat_estimates),2),
                                          r1=r1,
                                          r2=r2,
                                          p1=p1,
                                          p2=p2,
                                          power1=round(power_estimates["1wave"],2),
                                          power2=round(power_estimates["2_exch"],2),
                                          from_marg_1=as.numeric(predicted_power_margprob),
                                          from_cond_1=as.numeric(predicted_power_condprob),
                                          from_marg_2=as.numeric(predicted_power_margprob_2_wave),
                                          from_cond_2=as.numeric(predicted_power_condprob_2_wave),
                                          marg_power_off_1_wave=as.numeric(power_estimates["1wave"]-predicted_power_margprob),
                                          marg_power_off_2_waves=as.numeric(power_estimates["2_exch"]-predicted_power_margprob_2_wave),
                                          cond_power_off_1_wave=as.numeric(power_estimates["1wave"]-predicted_power_condprob),
                                          cond_power_off_2_waves=as.numeric(power_estimates["2_exch"]-predicted_power_condprob_2_wave),
                                          power_gain_add_wave=round(power_gain,2))); 
    }
  }
}
rownames(overall_results) <- NULL;
write.csv(x=overall_results,file="results.csv");
print(overall_results);
finish_time <- Sys.time();
print(finish_time-start_time);
save.image(paste("Ran-",
                 as.integer(Sys.time()),
                 sep=""))