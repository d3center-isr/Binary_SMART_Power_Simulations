#####################################################
# Set up simulation
#####################################################
rm(list = ls());
set.seed(16801);
library(bindata);
library(geepack); 
source("SimulateBinarySmartTwoWaves.R");
#debug(simulate_binary_smart_2_waves);
n_sims <- 2000;
#########################################################################
# Define parameters;
#########################################################################
target_power <- .80;
target_alpha <- .05;
rho <- 0.5;  
prob_R_given_A1_minus <-   .6;
prob_R_given_A1_plus <-   .7; 
prob_R <- (prob_R_given_A1_minus+prob_R_given_A1_plus)/2;
gee_formula_1 <- y ~ a1 + a2 + a1:a2;
gee_formula_2 <- y ~  time + time:a1 +time:a2 + time:a1:a2;
n_params <- 7;
n_regimens <- 4;
n_contrasts <- 4;
regimens_1_wave_model <- rbind(pp=c(1,1,1,1),
                               pm=c(1,1,-1,-1),
                               mp=c(1,-1,1,-1),
                               mm=c(1,-1,-1,1));
colnames(regimens_1_wave_model) <- c("int","a1","a2","a1a2");
regimens_2_wave_model <- rbind(pp=c(1,1,1,1,1),
                               pm=c(1,1,1,-1,-1),
                               mp=c(1,1,-1,1,-1),
                               mm=c(1,1,-1,-1,1));
colnames(regimens_2_wave_model) <- c("int","t","a1","a2","a1a2")
# Define data structures to hold results;
effect_size_coefficients <- log(c(1.435,1.92,2.89));  
n_effect_size_coefficients <- length(effect_size_coefficients)
n_sample_sizes <- 2;
sample_sizes <- c(300,500);
predicted_power_approx_1_wave <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_sharp_1_wave <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_approx_2_wave <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_sharp_2_wave <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
mean_corr_conditional <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
all_rho_hat <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
contrasts_1_wave <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_1_wave <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
contrasts_2_waves_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_2_waves_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
contrasts_2_waves_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_2_waves_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
##################################################################
# Decide which contrast to focus on;
##################################################################
regimen1 <- "pm";
regimen1_responders <- "post_p1";
regimen1_nonresponders <- "post_p0m";
regimen2 <- "mm";
regimen2_responders <- "post_m1";
regimen2_nonresponders <- "post_m0m";
contrast_coefficients_1_wave <- regimens_1_wave_model[regimen1,]-
  regimens_1_wave_model[regimen2,];
contrast_coefficients_2_wave <- regimens_2_wave_model[regimen1,]-
  regimens_2_wave_model[regimen2,];
dummy_contrast_coefficients <- c(pp=0,pm=1,mp=0,mm=-1);
##################################################################
# Do main loop;
##################################################################
start_time <- Sys.time();
for (effect_size_scenario in 1:n_effect_size_coefficients) {
  print(paste("Starting effect size scenario",effect_size_scenario));
  effect_size_coefficient <- effect_size_coefficients[effect_size_scenario];
  #######################################################################################
  # Define parameters (conditional means and response rates) to begin simulation;
  #######################################################################################
  true_cell_conditional_means <- c(pre_0  = .70,  # baseline, future nonresponders;
                                   pre_1  = .70,  # baseline, future responders;
                                   post_p1=  plogis(1.0-.5*effect_size_coefficient),  # final, given A1=+1 and responder
                                   post_p0p= plogis(1.5-.5*effect_size_coefficient),  # time 1, given A1=+1 and nonresponder, given A2=+1
                                   post_p0m= plogis(1.5-.5*effect_size_coefficient),  # time 1, given A1=+1 and nonresponder, given A2=-1
                                   post_m1=  plogis(1.0+.5*effect_size_coefficient),  # time 1, given A1=-1 and responder
                                   post_m0p= plogis(1.5+.5*effect_size_coefficient),  # time 1, given A1=-1 and nonresponder, given A2=+1
                                   post_m0m= plogis(1.5+.5*effect_size_coefficient));  # time 1, given A1=-1 and nonresponder, given A2=-1 
  print(true_cell_conditional_means);
  #########################################################################
  # Define true values of marginal means and marginal estimands.;
  #########################################################################
  n_regimens <- 4;
  p <- true_cell_conditional_means;  # just to make the formulas below look shorter;
  true_cell_conditional_variances <- true_cell_conditional_means * 
    (1-true_cell_conditional_means);
  marginal_pretest_mean <- as.numeric(prob_R*p["pre_1"]+
                                        (1-prob_R)*p["pre_0"]);
  true_regimen_marginal_means <- c(
    pp=as.numeric(prob_R_given_A1_plus*p["post_p1"]+ 
                    (1-prob_R_given_A1_plus)*p["post_p0p"]),
    pm=as.numeric(prob_R_given_A1_plus*p["post_p1"]+ 
                    (1-prob_R_given_A1_plus)*p["post_p0m"]),
    mp=as.numeric(prob_R_given_A1_minus*p["post_m1"]+ 
                    (1-prob_R_given_A1_minus)*p["post_m0p"]),
    mm=as.numeric(prob_R_given_A1_minus*p["post_m1"]+ 
                    (1-prob_R_given_A1_minus)*p["post_m0m"]));
  print(true_regimen_marginal_means);
  # The suffixes denote (A1, A2NR).
  # For example, t3mp means A1=-1 for everyone and then
  # A2=+1 for nonresponders.
  linear_weights_for_last_wave <- c(0,0,1);
  true_regimen_marginal_means <- c(
    pp=as.numeric(true_regimen_marginal_means["pp"]),
    pm=as.numeric(true_regimen_marginal_means["pm"]),
    mp=as.numeric(true_regimen_marginal_means["mp"]),
    mm=as.numeric(true_regimen_marginal_means["mm"]));
  regimens_true_final_responder_probabilities <- 
    c(pp=as.numeric(p["post_p1"]),
      pm=as.numeric(p["post_p1"]),
      mp=as.numeric(p["post_m1"]),
      mm=as.numeric(p["post_m1"]));
  regimens_true_final_nonresponder_probabilities <- 
    c(pp=as.numeric(p["post_p0p"]),
      pm=as.numeric(p["post_p0m"]),
      mp=as.numeric(p["post_m0p"]),
      mm=as.numeric(p["post_m0m"]));
  regimens_response_rates <- c(
    prob_R_given_A1_plus,
    prob_R_given_A1_plus,
    prob_R_given_A1_minus,
    prob_R_given_A1_minus);
  ##################################################################
  # Begin loop over scenarios
  ##################################################################
  for (sample_size_scenario in 1:n_sample_sizes) {
    print(paste("Starting sample size scenario",sample_size_scenario));
    n_subjects <- sample_sizes[sample_size_scenario]; 
    # Calculate expected power using conservative approximation from Kidwell et al (2018)
    p1 <- true_regimen_marginal_means[regimen1];
    p2 <- true_regimen_marginal_means[regimen2]; 
    OR <- as.numeric((p1/(1-p1)) / (p2/(1-p2))); 
    print(paste("Odds ratio target is",max(OR,1/OR)));
    zQ <- qnorm(target_power);
    z_one_minus_half_alpha <- qnorm(1-target_alpha/2);
    var1_approx_1_wave <- as.numeric(2*((2-prob_R_given_A1_plus)/(p1*(1-p1))+
                                          (2-prob_R_given_A1_minus)/(p2*(1-p2))));
    predicted_power_approx_1_wave[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_approx_1_wave) -
              z_one_minus_half_alpha );
    ##################################################################
    # Try other power formulas
    ################################################################## 
    # Sharp one-wave:
    V1 <- p1*(1-p1);
    V2 <- p2*(1-p2);
    V11 <- true_cell_conditional_variances[regimen1_responders] +
      (regimens_true_final_responder_probabilities[regimen1]-
         true_regimen_marginal_means[regimen1])^2;
    V10 <- true_cell_conditional_variances[regimen1_nonresponders] +
      (regimens_true_final_nonresponder_probabilities[regimen1]-
         true_regimen_marginal_means[regimen1])^2;
    V21 <- true_cell_conditional_variances[regimen2_responders] +
      (regimens_true_final_responder_probabilities[regimen2]-
         true_regimen_marginal_means[regimen2])^2;
    V20 <- true_cell_conditional_variances[regimen2_nonresponders] +
      (regimens_true_final_nonresponder_probabilities[regimen2]-
         true_regimen_marginal_means[regimen2])^2;
    var1_contrast_sharp_1_wave <- (4*(1-prob_R_given_A1_plus)*V10 + 2*prob_R_given_A1_plus*V11)/(V1^2) + 
      (4*(1-prob_R_given_A1_minus)*V20 + 2*prob_R_given_A1_minus*V21)/(V2^2); 
    predicted_power_sharp_1_wave[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_sharp_1_wave) -
              z_one_minus_half_alpha);
    # Approximate two-wave power formula;
    marginal_pretest_variance <- marginal_pretest_mean * (1-marginal_pretest_mean);
    true_regimen_marginal_variance <- true_regimen_marginal_means * (1-true_regimen_marginal_means);
    temp <- rbind(c(0,
                    -rho*sqrt(marginal_pretest_variance*true_regimen_marginal_variance["pp"]),
                    -rho*sqrt(marginal_pretest_variance*true_regimen_marginal_variance["pm"]),
                    -rho*sqrt(marginal_pretest_variance*true_regimen_marginal_variance["mp"]),
                    -rho*sqrt(marginal_pretest_variance*true_regimen_marginal_variance["mm"])),
                  0,0,0,0);
    arrowhead_bread_matrix <- temp + t(temp) + 
      diag(c(4*marginal_pretest_variance,
             true_regimen_marginal_variance["pp"],
             true_regimen_marginal_variance["pm"],
             true_regimen_marginal_variance["mp"],
             true_regimen_marginal_variance["mm"]));
    inv_bread_matrix_two_waves <- (1-rho^2)*solve(arrowhead_bread_matrix);
    prob_R <- .5*prob_R_given_A1_plus + .5*prob_R_given_A1_minus;
    var1_contrast_approx_2_wave <- (4-2*prob_R) * c(0,dummy_contrast_coefficients) %*%
      inv_bread_matrix_two_waves %*% c(0,dummy_contrast_coefficients) ;
    predicted_power_approx_2_wave[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_approx_2_wave) -
              z_one_minus_half_alpha);
    # Sharp two-wave power formula;
    Rmat <- rbind(c(1,rho),c(rho,1));
    Rinv <- solve(Rmat);
    v0 <- marginal_pretest_mean*(1-marginal_pretest_mean);
    arrowhead_meat_matrix <- matrix(0,5,5);
    for (d in 1:4) {
      Gd <- diag(c(v0,
                   true_regimen_marginal_means[d] * 
                     (1-true_regimen_marginal_means[d])));
      Gd0 <- diag(c(v0,regimens_true_final_nonresponder_probabilities[d] * 
                      (1-regimens_true_final_nonresponder_probabilities[d])));
      Gd1 <- diag(c(v0,regimens_true_final_responder_probabilities[d] * 
                      (1-regimens_true_final_responder_probabilities[d])));
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
    predicted_power_sharp_2_wave[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_sharp_2_wave) -
              z_one_minus_half_alpha); 
    ##################################################################
    # Begin loop over simulations
    ##################################################################
    print("Begin simulation");
    for (this_sim in 1:n_sims) {
      simulation <- simulate_binary_smart_2_waves(n_subjects,
                                                  true_cell_conditional_means,
                                                  prob_R_given_A1_plus,
                                                  prob_R_given_A1_minus,
                                                  rho);
      sim_data <- simulation$long_data;
      mean_corr_conditional[effect_size_scenario,sample_size_scenario,this_sim] <- simulation$mean_corr_conditional;
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
      negative_pseudodata$wave <- negative_pseudodata$wave + 2;  
      # We keep the same subject ID to show that we don't really have all those
      # new participants.  So we have to distinguish the new observations somehow,
      # and so we treat them as new waves of data on the same person.  Although
      # it seems very ad-hoc, this method has been shown to be valid.
      # Create the final analysis dataset including replicates.
      data_for_analysis_2_waves <- rbind(positive_pseudodata,
                                         negative_pseudodata, 
                                         rows_not_to_replicate);
      data_for_analysis_2_waves <- data_for_analysis_2_waves[order(data_for_analysis_2_waves$subject_id,
                                                                   data_for_analysis_2_waves$wave),]; 
      data_for_analysis_1_wave <- data_for_analysis_2_waves[which(data_for_analysis_2_waves$time==1),];
      # Fit one-wave model:
      gee_1_wave <- geeglm(formula = gee_formula_1,  
                           id=subject_id,   
                           weights = known_weight,
                           data=data_for_analysis_1_wave,
                           corstr = "independence",
                           family=binomial()); 
      # Test contrast under one-wave GEE;
      estimated_log_odds_by_time_and_regimen_1_wave <- 
        drop(regimens_1_wave_model %*% gee_1_wave$geese$beta);
      contrasts_1_wave[effect_size_scenario, sample_size_scenario,this_sim] <-
        estimated_log_odds_by_time_and_regimen_1_wave[regimen1] - 
        estimated_log_odds_by_time_and_regimen_1_wave[regimen2] ;
      log_odds_ratio_var_1_wave <-
        dummy_contrast_coefficients %*% t(regimens_1_wave_model) %*%
        gee_1_wave$geese$vbeta %*%
        regimens_1_wave_model %*% dummy_contrast_coefficients; 
      std_err_contrasts_1_wave[effect_size_scenario, sample_size_scenario,this_sim] <- 
        sqrt(log_odds_ratio_var_1_wave);
      # Test contrast under working independence two-wave GEE; 
      gee_2_waves_ind <- geeglm(formula = gee_formula_2,  
                                id=subject_id,   
                                weights = known_weight,
                                data=data_for_analysis_2_waves,
                                corstr = "independence",
                                family=binomial()); 
      estimated_log_odds_by_time_and_regimen_2_waves_ind <- 
        drop(regimens_2_wave_model %*% gee_2_waves_ind$geese$beta);
      contrasts_2_waves_ind[effect_size_scenario, sample_size_scenario,this_sim] <-
        estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen1] - 
        estimated_log_odds_by_time_and_regimen_2_waves_ind[regimen2] ;
      log_odds_ratio_var_2_waves_ind <-
        dummy_contrast_coefficients %*% regimens_2_wave_model %*%
        gee_2_waves_ind$geese$vbeta %*%
        t(regimens_2_wave_model) %*% dummy_contrast_coefficients; 
      std_err_contrasts_2_waves_ind[effect_size_scenario, sample_size_scenario,this_sim] <- 
        sqrt(log_odds_ratio_var_2_waves_ind);
      
      # # Estimate exchangeable correlation parameter
      raw_residuals <- drop(gee_2_waves_ind$y - gee_2_waves_ind$fitted.values);
      residual_variances <- gee_2_waves_ind$fitted.values * 
        (1-gee_2_waves_ind$fitted.values);
      scaled_residuals <- raw_residuals / sqrt(residual_variances);
      all_crossprods <- NULL;
      for (this_subject in 1:n_subjects) { # subject index;
        these_rows <- which(data_for_analysis_2_waves$subject_id==this_subject);
        resids <- scaled_residuals[these_rows];
        all_crossprods <- c(all_crossprods,
                            resids[1]*resids[2]);
      }
      rho_hat <- pmax(0,mean(all_crossprods));
      all_rho_hat[effect_size_scenario,sample_size_scenario,this_sim] <- rho_hat;
      # Fit working exchangeable GEE:
      block_work_corr <- rbind(c(1,rho_hat),c(rho_hat,1));
      work_corr <- rbind(cbind(block_work_corr,0*block_work_corr),
                         cbind(0*block_work_corr,block_work_corr)); 
      work_corr_as_zcor <- fixed2Zcor(cor.fixed=work_corr, 
                                      id=data_for_analysis_2_waves$subject_id,  
                                      waves=data_for_analysis_2_waves$wave+1); 
      gee_2_waves_exch <- geeglm(formula = gee_formula_2,  
                                 id=subject_id,  
                                 family=binomial,
                                 weights = known_weight,
                                 data=data_for_analysis_2_waves,
                                 scale.fix=TRUE,
                                 corstr = "fixed",
                                 zcor=work_corr_as_zcor);
      estimated_log_odds_by_time_and_regimen_2_waves_exch <- 
        drop(regimens_2_wave_model %*% gee_2_waves_exch$geese$beta);
      contrasts_2_waves_exch[effect_size_scenario, sample_size_scenario,this_sim] <-
        estimated_log_odds_by_time_and_regimen_2_waves_exch[regimen1] - 
        estimated_log_odds_by_time_and_regimen_2_waves_exch[regimen2] ;
      log_odds_ratio_var_2_waves_exch <-
        dummy_contrast_coefficients %*% regimens_2_wave_model %*%
        gee_2_waves_exch$geese$vbeta %*%
        t(regimens_2_wave_model) %*% dummy_contrast_coefficients; 
      std_err_contrasts_2_waves_exch[effect_size_scenario, sample_size_scenario,this_sim] <- 
        sqrt(log_odds_ratio_var_2_waves_exch);
    }
  }
}
finish_time <- Sys.time();
power_estimates <- cbind(
  approx_OR = rep(exp(effect_size_coefficients), each=2), 
  n = rep(sample_sizes,times=3),
  pred_approx_1 = as.vector(t(predicted_power_approx_1_wave)),
  pred_sharp_1 = as.vector(t(predicted_power_sharp_1_wave)),
  obs_1 = as.vector(t(apply(abs(contrasts_1_wave) > 1.96* std_err_contrasts_1_wave,c(1,2),mean))),
  obs_2_ind = as.vector(t(apply(abs(contrasts_2_waves_ind) > 1.96* std_err_contrasts_2_waves_ind,c(1,2),mean))),
  pred_approx_2 = as.vector(t(predicted_power_approx_2_wave)),
  pred_sharp_2 = as.vector(t(predicted_power_sharp_2_wave)),
  obs_2_exch = as.vector(t(apply(abs(contrasts_2_waves_exch) > 1.96* std_err_contrasts_2_waves_exch,c(1,2),mean)))
);

print(round(power_estimates,4)); 
print(finish_time-start_time);

save.image(paste("Ran-simple-",as.integer(Sys.time()),".rdata",sep=""));
