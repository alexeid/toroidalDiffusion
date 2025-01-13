###
#

setwd("~/WorkSpace/toroidalDiffusion/experiments/BivariateDist")

library(tidyverse)

# contour plot for bivariate normal data
contour <- function(bvn, x_string, y_string, size, alpha){
  require("ggplot2")
  p <- ggplot(bvn, aes_string(x = x_string, y = y_string)) +
    geom_point(size=size, alpha = alpha) +
    geom_density2d() +
    labs(title=paste0("mu1 = ", round(mu1,2), ", mu2 = ", round(mu2,2), ", sampling ", N, " points")) +
    theme_bw()
  return(p)
}

wrap_to_torus <- function(df, RANGE = 2 * pi) {
  #if (ncol(df) != 2) stop("Input data frame must have exactly two columns.")
  # Compute modulo 2*pi to apply wrapping
  return(df %% RANGE)
}

N <- 1000 # Number of random samples
set.seed(777)
# 𝜙, 𝜓.  phi, psi
mu1 <- pi * 0.65; s1 <- 1.5
mu2 <- pi * 0.8; s2 <- 1.5
alpha <- c(1.0, 0.5, 0.5)
# A = (a1, a3 * s1/s2; a3 * s2/s1, a2)
A <- matrix(c(alpha[1], alpha[3]*s1/s2, alpha[3]*s2/s1, alpha[2]), 2)
Sigma = diag(s1^2, s2^2)
# TODO covariance matrix is A^−1 * Σ or 1/2 * A^−1 * Σ ?
covarM = 1/2 * solve(A) %*% Sigma

# X ∼ N (μ, 1/2 * A^−1 * Σ)

library(MASS)
bvn1 <- mvrnorm(N, mu = c(mu1, mu2), Sigma = covarM ) # from MASS package
colnames(bvn1) <- c("phi","psi")
p1 <- contour(as_tibble(bvn1), "phi","psi", size=.5, alpha = .3)
p1
ggsave("bivariate.pdf", p1, width = 6, height = 6)


# TODO wrapped or re-sample
bvn2 <- wrap_to_torus(bvn1)
# bvn2
p2 <- contour(as_tibble(bvn2), "phi","psi", size=.5, alpha = .3)
p2
ggsave("wrapped-bivariate.pdf", p2, width = 6, height = 6)

