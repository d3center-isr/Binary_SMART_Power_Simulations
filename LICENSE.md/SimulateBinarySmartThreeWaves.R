simulate_binary_smart_3_waves <- function(n_subjects,
                                          true_cell_conditional_means,
                                          prob_R_given_A1_plus,
                                          prob_R_given_A1_minus,
                                          true_correlation_matrix) {
  #######################################################################################
  # Randomly generate data
  #######################################################################################
  # Step One: Using principal stratification (subgroups)
  #Subgroup	                        Response status 	Response status 
  #                                     if A_1=+1         if A_1=-1
  #1. Would always respond	               R=1	            R=1
  #2. Would only respond to a_1=+1         R=1            	R=0
  #3. Would only respond to a_1=-1	       R=0             	R=1
  #4. Would never respond	                 R=0	            R=0
  n4 <- round(min(n_subjects*(1-prob_R_given_A1_plus), n_subjects*(1-prob_R_given_A1_minus)));
  n1 <- round(n_subjects*(prob_R_given_A1_plus+prob_R_given_A1_minus-1))+n4;
  n2 <- round(n_subjects*(1-prob_R_given_A1_minus))-n4;
  n3 <- round(n_subjects*(1-prob_R_given_A1_plus))-n4;
  n4 <- n_subjects - (n1+n2+n3); # just to get rid of rounding error;
  subgroup <- sample(c(rep(1,n1),
                       rep(2,n2),
                       rep(3,n3),
                       rep(4,n4)));
  proportion_in_subgroups <- c(mean(subgroup==1),
                               mean(subgroup==2),
                               mean(subgroup==3),
                               mean(subgroup==4));
  # Step 2:  Apply randomization
  # Assign randomization to pick which one happens
  a1 <- sample(as.integer(c(-1,+1)),size=n_subjects, replace=TRUE);
  R <- rep(as.integer(NA), n_subjects);
  R[which(subgroup==1)] <- 1;
  R[which((subgroup==2)&(a1==-1))] <- 0;
  R[which((subgroup==2)&(a1==+1))] <- 1;
  R[which((subgroup==3)&(a1==-1))] <- 1;
  R[which((subgroup==3)&(a1==+1))] <- 0;
  R[which(subgroup==4)] <- 0;
  a2 <- vector("integer", n_subjects);
  a2[which(R==0)] <- sample(c(-1,+1),size=sum(R==0),replace=TRUE);
  a2[which(R==1)] <- 0;
  # Step 3: Simulate the outcomes;
  subject_level_design_matrix <- cbind(subject_id = 1:n_subjects,
                                       a1=a1,
                                       R=R,
                                       a2=a2);
  observation_level_design_matrix <- rbind(cbind(time=1,subject_level_design_matrix),
                                           cbind(time=2,subject_level_design_matrix),
                                           cbind(time=3,subject_level_design_matrix));
  all_corrs_conditional <- NULL;
  wide_mu <- matrix(NA,n_subjects,3);
  wide_y <- matrix(NA,n_subjects,3);
  time2_cell_name <- paste(ifelse(a1>0,"p","m"),
                           R,
                           sep="_");
  time3_cell_name <- paste(ifelse(a1>0,"p","m"),
                           R,
                           c("m","?","p")[a2+2],
                           sep="_");
  p <- true_cell_conditional_means;  # just to make the formulas below look shorter;
  for (this_time3_cell_name in unique(time3_cell_name)) { 
    these_subjects <- which(time3_cell_name == this_time3_cell_name); 
    if (this_time3_cell_name=="p_1_?") { mu_this_cell <- c(p["t1"],p["t2p"],p["t3p1"])}
    if (this_time3_cell_name=="p_0_p") { mu_this_cell <- c(p["t1"],p["t2p"],p["t3p0p"])}
    if (this_time3_cell_name=="p_0_m") { mu_this_cell <- c(p["t1"],p["t2p"],p["t3p0m"])}
    if (this_time3_cell_name=="m_1_?") { mu_this_cell <- c(p["t1"],p["t2m"],p["t3m1"])}
    if (this_time3_cell_name=="m_0_p") { mu_this_cell <- c(p["t1"],p["t2m"],p["t3m0p"])}
    if (this_time3_cell_name=="m_0_m") { mu_this_cell <- c(p["t1"],p["t2m"],p["t3m0m"])}
    y_these_subjects <- suppressWarnings(rmvbin(length(these_subjects),
                                                margprob=mu_this_cell,
                                                bincorr=true_correlation_matrix)); 
    wide_mu[these_subjects,] <- rep(1,length(these_subjects)) %o% mu_this_cell;
    wide_y[these_subjects,] <- y_these_subjects;
    y_corrs_conditional <- suppressWarnings(cor(y_these_subjects)[upper.tri(cor(y_these_subjects))]);
    all_corrs_conditional <- c(all_corrs_conditional, y_corrs_conditional);
  } 
  long_data <- data.frame(subject_id=rep(1:n_subjects,each=3),
                          a1=rep(a1,each=3),
                          R=rep(R,each=3),
                          a2=rep(a2,each=3),
                          time2_cell_name=rep(time2_cell_name,each=3),
                          time3_cell_name=rep(time3_cell_name,each=3),
                          time=rep(1:3,times=n_subjects),
                          y=as.vector(t(wide_y)),
                          mu=as.vector(t(wide_mu)));
  long_data$s1 <- 1*(long_data$time>1); # time spent in stage 1;
  long_data$s2 <- 1*(long_data$time>2); # time spent in stage 2;
  return(list(long_data=long_data,
              mean_corr_conditional=mean(all_corrs_conditional)));
}
 