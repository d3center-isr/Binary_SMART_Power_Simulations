simulate_data_logistic <- function( n_subjects, 
                            E_Y0,
                            theta_R,
                            beta_Y2) { 
  Y0 <- rbinom(n_subjects,1,E_Y0);
  A1 <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  A2NR <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  E_R_i <- plogis(theta_R["intercept"] +
                    theta_R["y0"]*Y0 +
                    theta_R["a1"]*A1);
  R <- rbinom(n_subjects,1,E_R_i);
  A2 <- (1-R)*A2NR;
  E_Y2 <- plogis(beta_Y2["intercept"] +
                     beta_Y2["y0"]*Y0 + 
                     beta_Y2["r"]*R+
                     beta_Y2["a1"]*A1+
                     beta_Y2["a2NR"]*A2 + 
                     beta_Y2["a1a2NR"]*A1*A2);
  Y2 <- rbinom(length(E_Y2),1,E_Y2);
  sim_data_wide <- data.frame(id=1:n_subjects,
                              Y0,
                              A1,
                              R,
                              A2,
                              Y2);
  return(sim_data_wide);
}