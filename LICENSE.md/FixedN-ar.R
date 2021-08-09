#####################################################
# Set up simulation
#####################################################
rm(list = ls());
set.seed(16802);
library(bindata);
library(geepack); 
source("SimulateBinarySmartThreeWaves.R");
n_sims <- 2000;
#########################################################################
# Define parameters;
#########################################################################
target_power <- .80;
target_alpha <- .05;
rho <- 0.3;  
true_corr_matrix <- rbind(c(1,rho,rho^2),c(rho,1,rho),c(rho^2,rho,1));
prob_R_given_A1_minus <-   .6;
prob_R_given_A1_plus <-   .7; 
prob_R <- (prob_R_given_A1_minus+prob_R_given_A1_plus)/2;
gee_formula <- y ~  s1 +
  s2 +
  s1:a1 +  
  s2:a1 +
  s2:a2 + 
  s2:a1:a2;
n_params <- 7;
n_regimens <- 4;
n_contrasts <- 4;
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

# Define data structures to hold results;
n_effect_size_coefficients <- 4;
effect_size_coefficients <- c(0,log(c(1.435,1.92,2.89)));  
n_sample_sizes <- 2;
sample_sizes <- c(300,500);
predicted_power_Kidwell <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_sharp <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_two_waves_approx_first_way <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_two_waves_approx_second_way <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_two_waves_pretty <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
predicted_power_two_waves_elaborate <- matrix(NA,n_effect_size_coefficients, n_sample_sizes);
mean_corr_conditional <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
all_rho_hat_ar1 <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
all_rho_hat_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
params_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims,n_params));
std_err_params_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims,n_params));
params_ar1 <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims, n_params));
std_err_params_ar1 <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims,n_params));
params_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims, n_params));
std_err_params_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims, n_params));
contrasts_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_ind <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
contrasts_ar1 <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_ar1 <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
contrasts_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
std_err_contrasts_exch <- array(NA,c(n_effect_size_coefficients, n_sample_sizes, n_sims));
##################################################################
# Decide which contrast to focus on;
##################################################################
regimen1 <- "pm";
regimen1_responders <- "t3p1";
regimen1_nonresponders <- "t3p0m";
regimen2 <- "mm";
regimen2_responders <- "t3m1";
regimen2_nonresponders <- "t3m0m";
contrast_coefficients <- c(0,1,0,-1);
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
  true_cell_conditional_means <- c(t1  = .70,  # time 1;
                                   t2p= .55,  # time 2, given A1=+1 
                                   t2m= .60,  # time 2, given A1=-1 
                                   t3p1=  plogis(1.0-.5*effect_size_coefficient),  # time 3, given A1=+1 and responder
                                   t3p0p= plogis(1.5-.5*effect_size_coefficient),  # time 3, given A1=+1 and nonresponder, given A2=+1
                                   t3p0m= plogis(1.5-.5*effect_size_coefficient),  # time 3, given A1=+1 and nonresponder, given A2=-1
                                   t3m1=  plogis(1.0+.5*effect_size_coefficient),  # time 3, given A1=-1 and responder
                                   t3m0p= plogis(1.5+.5*effect_size_coefficient),  # time 3, given A1=-1 and nonresponder, given A2=+1
                                   t3m0m= plogis(1.5+.5*effect_size_coefficient));  # time 3, given A1=-1 and nonresponder, given A2=-1 
  print(true_cell_conditional_means);
  #########################################################################
  # Define true values of marginal means and marginal estimands.;
  #########################################################################
  n_regimens <- 4;
  p <- true_cell_conditional_means;  # just to make the formulas below look shorter;
  true_cell_conditional_variances <- true_cell_conditional_means * 
    (1-true_cell_conditional_means);
  true_regimen_marginal_means <- c(
    t1=as.numeric(p["t1"]),
    t2p=as.numeric(prob_R_given_A1_plus*p["t2p"]+
                     (1-prob_R_given_A1_plus)*p["t2p"]),
    t2m=as.numeric(prob_R_given_A1_minus*p["t2m"]+
                     (1-prob_R_given_A1_minus)*p["t2m"]),
    t3pp=as.numeric(prob_R_given_A1_plus*p["t3p1"]+ 
                      (1-prob_R_given_A1_plus)*p["t3p0p"]),
    t3pm=as.numeric(prob_R_given_A1_plus*p["t3p1"]+ 
                      (1-prob_R_given_A1_plus)*p["t3p0m"]),
    t3mp=as.numeric(prob_R_given_A1_minus*p["t3m1"]+ 
                      (1-prob_R_given_A1_minus)*p["t3m0p"]),
    t3mm=as.numeric(prob_R_given_A1_minus*p["t3m1"]+ 
                      (1-prob_R_given_A1_minus)*p["t3m0m"]));
  print(true_regimen_marginal_means);
  # The suffixes denote (A1, A2NR).
  # For example, t3mp means A1=-1 for everyone and then
  # A2=+1 for nonresponders.
  q <- true_regimen_marginal_means;
  linear_weights_for_last_wave <- c(0,0,1);
  regimens_true_final_marginal_probabilities <- c(
    pp=as.numeric(q["t3pp"]),
    pm=as.numeric(q["t3pm"]),
    mp=as.numeric(q["t3mp"]),
    mm=as.numeric(q["t3mm"]));
  regimens_true_final_responder_probabilities <- 
    c(
      pp=as.numeric(p["t3p1"]),
      pm=as.numeric(p["t3p1"]),
      mp=as.numeric(p["t3m1"]),
      mm=as.numeric(p["t3m1"]));
    regimens_true_final_nonresponder_probabilities <- 
      c(
        pp=as.numeric(p["t3p0p"]),
        pm=as.numeric(p["t3p0m"]),
        mp=as.numeric(p["t3m0p"]),
        mm=as.numeric(p["t3m0m"]));
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
    p1 <- regimens_true_final_marginal_probabilities[regimen1];
    p2 <- regimens_true_final_marginal_probabilities[regimen2];
    if (effect_size_scenario ==1) {
    OR <- NA;
    } else {
    OR <- as.numeric((p1/(1-p1)) / (p2/(1-p2)));
    }
    print(paste("Odds ratio target is",max(OR,1/OR)));
    zQ <- qnorm(target_power);
    z_one_minus_half_alpha <- qnorm(1-target_alpha/2);
    var1_contrast_Kidwell_formula <- as.numeric(2*((2-prob_R_given_A1_plus)/(p1*(1-p1))+
                                                     (2-prob_R_given_A1_minus)/(p2*(1-p2))));
    predicted_power_Kidwell[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_Kidwell_formula) -
              z_one_minus_half_alpha );
    ##################################################################
    # Try other power formulas
    ################################################################## 
    V1 <- p1*(1-p1);
    V2 <- p2*(1-p2);
    V11 <- true_cell_conditional_variances[regimen1_responders] +
      (regimens_true_final_responder_probabilities[regimen1]-
         regimens_true_final_marginal_probabilities[regimen1])^2;
    V10 <- true_cell_conditional_variances[regimen1_nonresponders] +
      (regimens_true_final_nonresponder_probabilities[regimen1]-
         regimens_true_final_marginal_probabilities[regimen1])^2;
    V21 <- true_cell_conditional_variances[regimen2_responders] +
      (regimens_true_final_responder_probabilities[regimen2]-
         regimens_true_final_marginal_probabilities[regimen2])^2;
    V20 <- true_cell_conditional_variances[regimen2_nonresponders] +
      (regimens_true_final_nonresponder_probabilities[regimen2]-
         regimens_true_final_marginal_probabilities[regimen2])^2;
    var1_contrast_sharp_one_wave <- (4*(1-prob_R_given_A1_plus)*V10 + 2*prob_R_given_A1_plus*V11)/(V1^2) + 
      (4*(1-prob_R_given_A1_minus)*V20 + 2*prob_R_given_A1_minus*V21)/(V2^2); 
    predicted_power_sharp[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_sharp_one_wave) -
              z_one_minus_half_alpha);
    var1_contrast_pretty <- 2*(2-prob_R)*((4-3*rho^2)/(4*p1*(1-p1)) -
                 (rho^2)/(2*sqrt(p1*(1-p1)*p2*(1-p2))) + 
                   (4-3*rho^2)/(4*p2*(1-p2)));
    predicted_power_two_waves_pretty[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_pretty) -
              z_one_minus_half_alpha);
    
    temp <- rbind(c(0,
                    -rho*sqrt(q["t1"]*(1-q["t1"])*q["t3pp"]*(1-q["t3pp"])),
                    -rho*sqrt(q["t1"]*(1-q["t1"])*q["t3pm"]*(1-q["t3pm"])),
                    -rho*sqrt(q["t1"]*(1-q["t1"])*q["t3mp"]*(1-q["t3mp"])),
                    -rho*sqrt(q["t1"]*(1-q["t1"])*q["t3mm"]*(1-q["t3mm"]))),
                  0,0,0,0);
    arrowhead_bread_matrix <- temp + t(temp) + 
      diag(c(4*(q["t1"]*(1-q["t1"])),
             (q["t3pp"]*(1-q["t3pp"])),
             (q["t3pm"]*(1-q["t3pm"])),
             (q["t3mp"]*(1-q["t3mp"])),
             (q["t3mm"]*(1-q["t3mm"]))));
    inv_bread_matrix_two_waves <- (1-rho^2)*solve(arrowhead_bread_matrix);
    prob_R <- .5*prob_R_given_A1_plus + .5*prob_R_given_A1_minus;
    var1_contrast_two_waves_approx_first_way <- (4-2*prob_R) * c(0,contrast_coefficients) %*%
       inv_bread_matrix_two_waves %*% c(0,contrast_coefficients) ;
    predicted_power_two_waves_approx_first_way[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_two_waves_approx_first_way) -
              z_one_minus_half_alpha);
              
    var1_contrast_two_waves_approx_second_way <- c(0,contrast_coefficients) %*%
      (diag(c(4,2*rep(2-prob_R,4))) %*% inv_bread_matrix_two_waves) %*% c(0,contrast_coefficients) ;
    predicted_power_two_waves_approx_second_way[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_two_waves_approx_second_way) -
              z_one_minus_half_alpha);
    Rmat <- rbind(c(1,rho),c(rho,1));
    Rinv <- solve(Rmat);
    v0 <- p["t1"]*(1-p["t1"]);
    arrowhead_meat_matrix <- matrix(0,5,5);
    for (d in 1:4) {
      Gd <- diag(c(v0,
        regimens_true_final_marginal_probabilities[d] * 
        (1-regimens_true_final_marginal_probabilities[d])));
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
    var1_contrast_two_waves_elaborate <-   c(0,contrast_coefficients) %*%
      inv_bread_matrix_two_waves %*% arrowhead_meat_matrix  %*%
      inv_bread_matrix_two_waves %*% c(0,contrast_coefficients) ;
    predicted_power_two_waves_elaborate[effect_size_scenario,sample_size_scenario] <- 
      pnorm(sqrt(n_subjects*(log(OR)^2)/var1_contrast_two_waves_elaborate) -
              z_one_minus_half_alpha); 
    ##################################################################
    # Begin loop over simulations
    ##################################################################
    print("Begin simulation");
    for (this_sim in 1:n_sims) {
      simulation <- simulate_binary_smart_3_waves(n_subjects,
                                                  true_cell_conditional_means,
                                                  prob_R_given_A1_plus,
                                                  prob_R_given_A1_minus,
                                                  true_corr_matrix);
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
      params_ind[effect_size_scenario, sample_size_scenario,this_sim,] <- gee_ind$geese$beta;
      std_err_params_ind[effect_size_scenario, sample_size_scenario,this_sim,] <- sqrt(diag(gee_ind$geese$vbeta));
      gee_cov_beta_ind <- gee_ind$geese$vbeta;
      # Test contrast under working independence GEE;
      estimated_log_odds_by_time_and_regimen_ind <- cbind(
        design_matrices_by_regimen[,,"pp"]%*%gee_ind$geese$beta,
        design_matrices_by_regimen[,,"pm"]%*%gee_ind$geese$beta,
        design_matrices_by_regimen[,,"mp"]%*%gee_ind$geese$beta,
        design_matrices_by_regimen[,,"mm"]%*%gee_ind$geese$beta);
      colnames(estimated_log_odds_by_time_and_regimen_ind) <- c("pp","pm","mp","mm");
      contrasts_ind[effect_size_scenario, sample_size_scenario,this_sim] <-
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
      std_err_contrasts_ind[effect_size_scenario, sample_size_scenario,this_sim] <- sqrt(log_odds_ratio_var_ind);
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
      all_rho_hat_exch[effect_size_scenario,sample_size_scenario,this_sim] <- rho_hat_exch;
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
      params_exch[effect_size_scenario, sample_size_scenario,this_sim,] <- gee_exch$geese$beta;
      std_err_params_exch[effect_size_scenario, sample_size_scenario,this_sim,] <- sqrt(diag(gee_exch$geese$vbeta));
      gee_cov_beta_exch <- gee_exch$geese$vbeta;
      # Test contrast under working exchangeable GEE;
      estimated_log_odds_by_time_and_regimen_exch <- cbind(
        design_matrices_by_regimen[,,"pp"]%*%gee_exch$geese$beta,
        design_matrices_by_regimen[,,"pm"]%*%gee_exch$geese$beta,
        design_matrices_by_regimen[,,"mp"]%*%gee_exch$geese$beta,
        design_matrices_by_regimen[,,"mm"]%*%gee_exch$geese$beta);
      colnames(estimated_log_odds_by_time_and_regimen_exch) <- c("pp","pm","mp","mm");
      contrasts_exch[effect_size_scenario, sample_size_scenario,this_sim] <-
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
      std_err_contrasts_exch[effect_size_scenario, sample_size_scenario,this_sim] <- sqrt(log_odds_ratio_var_exch);
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
      params_ar1[effect_size_scenario, sample_size_scenario,this_sim,] <- gee_ar1$geese$beta;
      std_err_params_ar1[effect_size_scenario, sample_size_scenario,this_sim,] <- sqrt(diag(gee_ar1$geese$vbeta));
      gee_cov_beta_ar1 <- gee_ar1$geese$vbeta;
      # Test contrast under working AR-1 GEE;
      estimated_log_odds_by_time_and_regimen_ar1 <- cbind(
        design_matrices_by_regimen[,,"pp"]%*%gee_ar1$geese$beta,
        design_matrices_by_regimen[,,"pm"]%*%gee_ar1$geese$beta,
        design_matrices_by_regimen[,,"mp"]%*%gee_ar1$geese$beta,
        design_matrices_by_regimen[,,"mm"]%*%gee_ar1$geese$beta);
      colnames(estimated_log_odds_by_time_and_regimen_ar1) <- c("pp","pm","mp","mm");
      contrasts_ar1[effect_size_scenario, sample_size_scenario,this_sim] <-
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
      std_err_contrasts_ar1[effect_size_scenario, sample_size_scenario,this_sim] <- sqrt(log_odds_ratio_var_ar1);
    }
  }
}
finish_time <- Sys.time();
power_estimates <- cbind(
  predicted_power_simple = as.vector(t(predicted_power_Kidwell)),
  predicted_power_sharp = as.vector(t(predicted_power_sharp)),
  predicted_power_two_waves_approx = as.vector(t(predicted_power_two_waves_approx_first_way)),
  predicted_power_two_waves_approx_second_way = as.vector(t(predicted_power_two_waves_approx_second_way)),
  predicted_power_two_waves_pretty = as.vector(t(predicted_power_two_waves_pretty)),
  predicted_power_two_waves_sharp = as.vector(t(predicted_power_two_waves_elaborate)),
  observed_power_ind = as.vector(t(apply(abs(contrasts_ind) > 1.96* std_err_contrasts_ind,c(1,2),mean))),
  observed_power_ar1 = as.vector(t(apply(abs(contrasts_ar1) > 1.96* std_err_contrasts_ar1,c(1,2),mean))),
  observed_power_exch = as.vector(t(apply(abs(contrasts_exch) > 1.96* std_err_contrasts_exch,c(1,2),mean)))
);

print(round(power_estimates,4));
print(summary(all_rho_hat_exch));
print(finish_time-start_time);


data_for_analysis$d0 <- 1*(data_for_analysis$time==1);
data_for_analysis$d1 <- 1*(data_for_analysis$time==3)*(data_for_analysis$a1==1)*(data_for_analysis$a2==1);
data_for_analysis$d2 <- 1*(data_for_analysis$time==3)*(data_for_analysis$a1==1)*(data_for_analysis$a2==-1);
data_for_analysis$d3 <- 1*(data_for_analysis$time==3)*(data_for_analysis$a1==-1)*(data_for_analysis$a2==1);
data_for_analysis$d4 <- 1*(data_for_analysis$time==3)*(data_for_analysis$a1==-1)*(data_for_analysis$a2==-1);
data_for_check <- data_for_analysis[which(data_for_analysis$time!=2),];
gee_to_check_formula <- geeglm(formula=y~d0+d1+d2+d3+d4+0,
                               id=subject_id,  
                               family=binomial,
                               weights = known_weight,
                               data=data_for_check,
                               scale.fix=TRUE,
                               corstr = "independence");
m1 <- summary(gee_to_check_formula)$cov.scaled;
m2 <- diag(c(4,2*rep(2-prob_R,4)))%*%inv_bread_matrix_two_waves / n_subjects;
m3 <- 2*(2-prob_R)*inv_bread_matrix_two_waves / n_subjects;

coef_std_err_estimates <- sqrt(cbind(diag(m1),
                                     diag(m2),
                                     diag(m3)));
print(coef_std_err_estimates);
  

save.image(paste("Ran",as.integer(Sys.time()),".rdata",sep=""));
