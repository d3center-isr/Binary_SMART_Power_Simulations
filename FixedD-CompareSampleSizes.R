#####################################################
# Set up simulation
#####################################################
rm(list = ls());
library(bindata);
library(geepack);
source("FindNBySimulation.R");
source("SimulateBinarySmartThreeWaves.R");
set.seed(16802); 
#########################################################################
# Define parameters;
#########################################################################
target_power <- .80;
target_alpha <- .05;
prob_R_given_A1_minus <-   .6;
prob_R_given_A1_plus <- .7; 
prob_R <- (prob_R_given_A1_minus+prob_R_given_A1_plus)/2;
rho <- 0.3;  
effect_size_coefficients <- log(c(1.435,1.92,2.89)); 
true_corr_structures <- c("ind","AR1","exch");
##################################################################
# Specify which contrast to focus on;
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
all_answers <- NULL;
all_detailed_answers <- list();
run <- 0;
for (effect_size_scenario in 1:3) {
  for (true_corr_scenario in 1:3) {
    print(c(effect_size_scenario, true_corr_scenario))
    run <- run + 1;
    effect_size_coefficient <- effect_size_coefficients[effect_size_scenario];
    this_corr_structure_name <- true_corr_structures[true_corr_scenario];
    ###################################################
    # Define parameters (conditional means and response rates) to begin simulation;
    true_cell_conditional_means <- c(t1  = .70,  # time 1;
                                     t2p= .55,  # time 2, given A1=+1 
                                     t2m= .60,  # time 2, given A1=-1 
                                     t3p1=  plogis(1.0-.5*effect_size_coefficient),  # time 3, given A1=+1 and responder
                                     t3p0p= plogis(1.5-.5*effect_size_coefficient),  # time 3, given A1=+1 and nonresponder, given A2=+1
                                     t3p0m= plogis(1.5-.5*effect_size_coefficient),  # time 3, given A1=+1 and nonresponder, given A2=-1
                                     t3m1=  plogis(1.0+.5*effect_size_coefficient),  # time 3, given A1=-1 and responder
                                     t3m0p= plogis(1.5+.5*effect_size_coefficient),  # time 3, given A1=-1 and nonresponder, given A2=+1
                                     t3m0m= plogis(1.5+.5*effect_size_coefficient));  # time 3, given A1=-1 and nonresponder, given A2=-1 
    ###################################################
    # Define true values of marginal means and marginal estimands.;
    n_regimens <- 4;
    p <- true_cell_conditional_means;  # just to make the formulas below look shorter;
    true_cell_conditional_variances <- true_cell_conditional_means * (1-true_cell_conditional_means);
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
      c(pp=as.numeric(p["t3p1"]),
        pm=as.numeric(p["t3p1"]),
        mp=as.numeric(p["t3m1"]),
        mm=as.numeric(p["t3m1"]));
    regimens_true_final_nonresponder_probabilities <- 
      c(pp=as.numeric(p["t3p0p"]),
        pm=as.numeric(p["t3p0m"]),
        mp=as.numeric(p["t3m0p"]),
        mm=as.numeric(p["t3m0m"]));
    regimens_response_rates <- c(
      prob_R_given_A1_plus,
      prob_R_given_A1_plus,
      prob_R_given_A1_minus,
      prob_R_given_A1_minus);
    ###################################################
    # Calculate needed sample size using...
    # ... conservative approximation from Kidwell et al (2018)
    p1 <- regimens_true_final_marginal_probabilities[regimen1];
    p2 <- regimens_true_final_marginal_probabilities[regimen2];
    log_OR <- log(as.numeric((p1/(1-p1)) / (p2/(1-p2))));
    zQ <- qnorm(target_power);
    z_one_minus_half_alpha <- qnorm(1-target_alpha/2);
    var1_log_OR_Kidwell <- as.numeric(2*((2-prob_R_given_A1_plus)/(p1*(1-p1))+
                                           (2-prob_R_given_A1_minus)/(p2*(1-p2))));
    n_required_Kidwell <- as.numeric(ceiling(((zQ+z_one_minus_half_alpha)^2)*
                                               (var1_log_OR_Kidwell)/
                                               (log_OR^2)));  
    # ... one-wave sharp formula;
    V1 <- p1*(1-p1);
    V2 <- p2*(1-p2);
    V11 <- true_cell_conditional_variances[regimen1_responders];
    V10 <- true_cell_conditional_variances[regimen1_nonresponders];
    V21 <- true_cell_conditional_variances[regimen2_responders];
    V20 <- true_cell_conditional_variances[regimen2_nonresponders];
    var1_log_OR_sharp_one_wave <- (4*(1-prob_R_given_A1_plus)*V10 + 2*prob_R_given_A1_plus*V11)/(V1^2) + 
      (4*(1-prob_R_given_A1_minus)*V20 + 2*prob_R_given_A1_minus*V21)/(V2^2); 
    n_required_sharp_one_wave <- as.numeric(ceiling(((zQ+z_one_minus_half_alpha)^2)*
                                                      (var1_log_OR_sharp_one_wave)/
                                                      (log_OR^2)));  
    # two-wave approximate formula:
    if (this_corr_structure_name=="ind") {rho <- 0;}
    if (this_corr_structure_name=="AR1") {rho <- .3^2}
    if (this_corr_structure_name=="exch") {rho <- .3}; 
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
    var1_log_OR_two_waves_approx <- (4-2*prob_R) * c(0,contrast_coefficients) %*%
      inv_bread_matrix_two_waves %*% c(0,contrast_coefficients);
    n_required_two_waves_approx <- as.numeric(ceiling(((zQ+z_one_minus_half_alpha)^2)*
                                                        (var1_log_OR_two_waves_approx)/
                                                        (log_OR^2)));  
    # two-wave sharp formula;
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
    var1_log_OR_two_waves_sharp <-   c(0,contrast_coefficients) %*%
      inv_bread_matrix_two_waves %*% arrowhead_meat_matrix  %*%
      inv_bread_matrix_two_waves %*% c(0,contrast_coefficients) ;
    n_required_two_waves_sharp <- as.numeric(ceiling(((zQ+z_one_minus_half_alpha)^2)*
                                                       (var1_log_OR_two_waves_sharp)/
                                                       (log_OR^2)));  
    ###################################################
    # Use simulation estimate;
    if (this_corr_structure_name=="ind") {
      true_corr_matrix <- diag(rep(1,3));
    }
    if (this_corr_structure_name=="ar1") {
      true_corr_matrix <- rbind(c(1,rho,rho^2),
                                c(rho,1,rho),
                                c(rho^2,rho,1));
    }
    if (this_corr_structure_name=="exch") {
      true_corr_matrix <- (1-rho)*diag(rep(1,3)) + rho*matrix(1,3,3);
    }
    print(true_cell_conditional_means);
    simulation_answer <- find_n_by_simulation( n_to_try=seq(.9*n_required_two_waves_sharp, 1.1*n_required_Kidwell,length=20),
                                               n_sims_each=  500,
                                               true_cell_conditional_means=true_cell_conditional_means,
                                               prob_R_given_A1_plus=prob_R_given_A1_plus,
                                               prob_R_given_A1_minus=prob_R_given_A1_minus,
                                               true_corr_matrix=true_corr_matrix);
    all_detailed_answers[[run]] <- simulation_answer;
    ###################################################
    # Record answers;
    these_answers <- data.frame(
      odds_ratio = exp(abs(log_OR)),
      true_corr = this_corr_structure_name,
      Kidwell = n_required_Kidwell,
      sharp_one_wave = n_required_sharp_one_wave,
      two_waves_approx = n_required_two_waves_approx,
      two_waves_sharp = n_required_two_waves_sharp,
      simulated_ind = simulation_answer$n_ind,
      simulated_ar1 = simulation_answer$n_ar1,
      simulated_exch = simulation_answer$n_exch );
    all_answers <- rbind(all_answers,these_answers)
  }
}
finish_time <- Sys.time();
print(all_answers);





