###
# Modified from
# https://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html#:~:text=Hence%2C%20a%20sample%20from%20a,variable%20conditioned%20on%20the%20first.
###

setwd("~/WorkSpace/toroidalDiffusion/experiments/BivariateDist")

N <- 1000 # Number of random samples
set.seed(123)
# Target parameters for univariate normal distributions
rho <- -0.6
mu1 <- 1; s1 <- 2
mu2 <- 1; s2 <- 8

# Parameters for bivariate normal distribution
mu <- c(mu1,mu2) # Mean
sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2), 2) # Covariance matrix

library(tidyverse)

# Function to draw ellipse for bivariate normal data
ellipse_bvn <- function(bvn, x_string, y_string, size, alpha){
  require("ggplot2")
  p <- ggplot(bvn, aes_string(x = x_string, y = y_string)) +
    geom_point(size=size, alpha = alpha) +
    geom_density2d() +
    labs(title=paste0("mu1 = ", mu1, ", mu2 = ", mu2, ", sampling ", N, " points")) +
    theme_bw()
  return(p)
}

### 1. MASS package

library(MASS)
bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X1","bvn1_X2")
ellipse_bvn(as_tibble(bvn1), "bvn1_X1","bvn1_X2", size=.5, alpha = .3)


### 2. Cholesky decomposition of sigma (a positive definite matrix) yields
# a matrix M such that M times its transpose gives sigma back again.

M <- t(chol(sigma))
# M %*% t(M)
Z <- matrix(rnorm(2*N),2,N) # 2 rows, N/2 columns
bvn2 <- t(M %*% Z) + matrix(rep(mu,N), byrow=TRUE,ncol=2)
colnames(bvn2) <- c("bvn2_X1","bvn2_X2")
ellipse_bvn(as_tibble(bvn2), "bvn2_X1","bvn2_X2", size=.5, alpha = .3)


### 3. first simulating a point from the marginal distribution
# of one of the random variables and then simulating from the second
# random variable conditioned on the first.

rbvn<-function (n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, mu1, s1)
  X2 <- rnorm(n, mu2 + (s2/s1) * rho *
                (X1 - mu1), sqrt((1 - rho^2)*s2^2))
  cbind(X1, X2)
}

bvn3 <- rbvn(N,mu1,s1,mu2,s2,rho)
colnames(bvn3) <- c("bvn3_X1","bvn3_X2")
ellipse_bvn(as_tibble(bvn3), "bvn3_X1","bvn3_X2", size=.5, alpha = .3)

