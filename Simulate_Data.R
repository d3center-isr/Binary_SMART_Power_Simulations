simulate_data <- function (  n_subjects,
                             mu_pre,
                             beta_R_0,
                             beta_R_Ypre,
                             beta_R_A1,
                             beta_Ypost_0,
                             beta_Ypost_A1,
                             beta_Ypost_Ypre,
                             beta_Ypost_R,
                             beta_Ypost_A2) {
  # For simplicity, this function uses the  assumption
  # that A2 does not have an effect!; 
  Ypre <- rbinom(n_subjects,1,mu_pre);
  A1 <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  A2NR <- sample(c(1,-1),size=n_subjects,replace=TRUE); 
  E_R_i <- beta_R_0 +
    beta_R_Ypre*Ypre +
    beta_R_A1*A1;
  R <- rbinom(n_subjects,1,E_R_i);
  A2 <- (1-R)*A2NR;
  E_Ypost_i <- beta_Ypost_0 +
    beta_Ypost_Ypre*Ypre + 
    beta_Ypost_A1*A1+
    beta_Ypost_R*R+
    beta_Ypost_A2*A2;
  Ypost <- rbinom(n_subjects,1,E_Ypost_i);
  sim_data_wide <- data.frame(id=1:n_subjects,
                              Ypre,
                              A1,
                              R,
                              A2,
                              Ypost);
  return(sim_data_wide);
}