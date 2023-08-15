test_that("StaPLR",{
  set.seed(012)
  n <- 100
  cors <- seq(0.1,0.7,0.1)
  X <- matrix(NA, nrow=n, ncol=length(cors)+1)
  X[,1] <- rnorm(n)
  
  for(i in 1:length(cors)){
    X[,i+1] <- X[,1]*cors[i] + rnorm(n, 0, sqrt(1-cors[i]^2))
  }
  
  beta <- c(1,0,0,0,0,0,0,0)
  eta <- X %*% beta
  p <- exp(eta)/(1+exp(eta))
  y <- rbinom(n, 1, p)
  view_index <- rep(1:(ncol(X)/2), each=2)
  
  StaPLR_fit <- StaPLR(X, y, view_index)
  expect_equal(coef(StaPLR_fit)$meta[1:5], c(-2.091923,  2.532491,  0.000000,  1.569726, 0.000000), tolerance = 1e-03)
  })