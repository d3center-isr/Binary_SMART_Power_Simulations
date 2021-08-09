rm(list = ls());
set.seed(16801);
library(bindata);
library(geepack); 
source("SimulateBinarySmartTwoWaves.R");
#########################################################################
# Define parameters;
#########################################################################
target_power <- .80;
target_alpha <- .05;
n_subjects <- 400;
effect_size_coefficient <- 1.92;
rho <- 0.3;   
n_sims <- 2000;
gee_formula_2 <- y ~  time + time:a1 +time:a2 + time:a1:a2;
n_params <- 5;
n_regimens <- 4;
n_contrasts <- 4; 
regimens_2_wave_model <- rbind(pp=c(1,1,1,1,1),
                               pm=c(1,1,1,-1,-1),
                               mp=c(1,1,-1,1,-1),
                               mm=c(1,1,-1,-1,1));
colnames(regimens_2_wave_model) <- c("int","t","a1","a2","a1a2")
# Define data structures to hold results;
R_probs <- c(.1,.3,.5,.7,.9);
n_R_probs <- 5;
coefs_independence <- array(NA, c(n_R_probs, n_sims, n_params));
coefs_exchangeable <- array(NA, c(n_R_probs, n_sims, n_params));
std_errs_independence <- array(NA, c(n_R_probs, n_sims, n_params));
std_errs_exchangeable <- array(NA, c(n_R_probs, n_sims, n_params));
##################################################################
# Do main loop;
##################################################################
start_time <- Sys.time();
for (R_prob_scenario in 1:n_R_probs) {
  prob_R <- R_probs[R_prob_scenario];
  print(paste("Starting scenario",R_prob_scenario));
  ##################################################################################
  # Define parameters (conditional means and response rates) to begin simulation;
  ##################################################################################
  true_cell_conditional_means <- c(pre_0  = .80,  # baseline, future nonresponders;
                                   pre_1  = .60,  # baseline, future responders;
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
    pp=as.numeric(prob_R*p["post_p1"]+ 
                    (1-prob_R)*p["post_p0p"]),
    pm=as.numeric(prob_R*p["post_p1"]+ 
                    (1-prob_R)*p["post_p0m"]),
    mp=as.numeric(prob_R*p["post_m1"]+ 
                    (1-prob_R)*p["post_m0p"]),
    mm=as.numeric(prob_R*p["post_m1"]+ 
                    (1-prob_R)*p["post_m0m"]));
  ##################################################################
  # Begin loop over simulations
  ##################################################################
  print("Begin simulation");
  for (this_sim in 1:n_sims) {
    simulation <- simulate_binary_smart_2_waves(n_subjects,
                                                true_cell_conditional_means,
                                                prob_R,
                                                prob_R,
                                                rho);
    sim_data <- simulation$long_data;
    ###########################################################################
    # Start analysis:
    ############################################################################
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
    # Test contrast under working independence two-wave GEE; 
    gee_2_waves_ind <- geeglm(formula = gee_formula_2,  
                              id=subject_id,   
                              weights = known_weight,
                              data=data_for_analysis_2_waves,
                              corstr = "independence",
                              family=binomial()); 
    
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
    ################################################
    # Record results
    #################################################
    coefs_independence[R_prob_scenario,this_sim,] <- summary(gee_2_waves_ind)$coefficients[,"Estimate"];
    std_errs_independence[R_prob_scenario,this_sim,] <- summary(gee_2_waves_ind)$coefficients[,"Std.err"];
    coefs_exchangeable[R_prob_scenario,this_sim,] <- summary(gee_2_waves_exch)$coefficients[,"Estimate"];
    std_errs_exchangeable[R_prob_scenario,this_sim,] <- summary(gee_2_waves_exch)$coefficients[,"Std.err"];
  }
}
finish_time <- Sys.time();
print(apply(coefs_independence,c(3,1),sd));
print(apply(coefs_exchangeable,c(3,1),sd));
print(apply(std_errs_independence,c(3,1),mean));
print(apply(std_errs_exchangeable,c(3,1),mean));
save.image(paste("Ran-StudySEs-",as.integer(Sys.time()),".rdata",sep=""));