simulate_data_logistic_3_wave <- function( n_subjects, 
                                    E_Y0,
                                    theta_R,
                                    beta_Y1,
                                    beta_Y2) { 
  Y0 <- rbinom(n_subjects,1,E_Y0);
  A1 <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  A2NR <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  E_R_i <- plogis(theta_R["intercept"] +
                    theta_R["y0"]*Y0 +
                    theta_R["a1"]*A1);
  R <- rbinom(n_subjects,1,E_R_i);
  A2 <- (1-R)*A2NR;
  E_Y1 <- plogis(beta_Y1["intercept"] +
                   beta_Y1["y0"]*Y0 + 
                   beta_Y1["r"]*R+
                   beta_Y1["a1"]*A1+
                   beta_Y1["a2NR"]*A2 + 
                   beta_Y1["a1a2NR"]*A1*A2);
  Y1 <- rbinom(length(E_Y1),1,E_Y1); 
  E_Y2_given_Y1 <- plogis(beta_Y2["intercept"] +
                            beta_Y2["a1"]*A1 + 
                            beta_Y2["y1"]*Y1);
  Y2 <- rbinom(n=length(E_Y2_given_Y1),size=1,prob=E_Y2_given_Y1);
  sim_data_wide <- data.frame(id=1:n_subjects,
                              Y0,
                              A1,
                              R,
                              A2,
                              Y1,
                              Y2);
  return(sim_data_wide);
}