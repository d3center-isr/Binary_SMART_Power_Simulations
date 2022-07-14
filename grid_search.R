target_mu <- .5;
target_rho <- .6;
X <- rbinom(n=100000,size=1,prob=target_mu);
minimand <- function(candidate_params) {
  beta0 <- candidate_params[1];
  beta1 <- candidate_params[2];
  E_Y_given_X <- plogis(beta0 + beta1*X);
  Y <- rbinom(n=length(E_Y_given_X),size=1,prob=E_Y_given_X);
  mu <- mean(Y);
  rho <- cor(X,Y);
  return((mu-target_mu)^2 + (rho-target_rho)^2);
}
minimizer <- optim(fn=minimand, par=c(0,0));
print(minimizer);