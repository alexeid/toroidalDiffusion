###
#

setwd("~/WorkSpace/TDExperiments/Model/fig1")

library(tidyverse)

# contour plot for bivariate normal data
contour <- function(bvn, x_string, y_string, size, alpha){
  require("ggplot2")
  p <- ggplot(bvn, aes_string(x = x_string, y = y_string)) +
    geom_point(size=size, alpha = alpha) +
    #geom_density2d() +
    labs(title=paste0("mu1 = ", 0, ", mu2 = ", 0, ", sampling ", N, " points")) +
    scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi),
      labels = expression(-pi, -pi/2, 0, pi/2, pi)
    ) +
    scale_y_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi),
      labels = expression(-pi, -pi/2, 0, pi/2, pi)
    ) +
    theme_bw()
  return(p)
}

wrap_to_torus <- function(df, RANGE = 2 * pi) {
  #if (ncol(df) != 2) stop("Input data frame must have exactly two columns.")
  # Compute modulo 2*pi to apply wrapping
  return(df %% RANGE)
}

N <- 100000 # Number of random samples
psize=.1; palpha = .1
set.seed(777)
# 𝜙, 𝜓.  phi, psi
mu1 <- pi; s1 <- 1.5
mu2 <- pi; s2 <- 1.5
alpha <- c(1.5, 1.0, -0.5)
# A = (a1, a3 * s1/s2; a3 * s2/s1, a2)
A <- matrix(c(alpha[1], alpha[3]*s1/s2, alpha[3]*s2/s1, alpha[2]), 2)
Sigma = diag(s1^2, s2^2)
# covariance matrix is 1/2 * A^−1 * Σ
covarM = 0.5 * solve(A) %*% Sigma
# solved by alpha and sd
covarM2 = 0.5 / (alpha[1]*alpha[2]-alpha[3]*alpha[3]) *
  matrix(c(alpha[2]*s1*s1, -alpha[3]*s1*s2, -alpha[3]*s1*s2, alpha[1]*s2*s2), 2)
stopifnot(covarM == covarM2)

# X ∼ N (μ, 1/2 * A^−1 * Σ)

### 1. MASS package
library(MASS)
bvn1 <- mvrnorm(N, mu = c(mu1, mu2), Sigma = covarM ) # from MASS package
colnames(bvn1) <- c("phi","psi")
p1 <- contour(as_tibble(bvn1), "phi","psi", size=psize, alpha=palpha) +
  xlim(-5, 10) + ylim(-7.5, 12.5)
p1
ggsave("bivariate.pdf", p1, width = 6, height = 6)

# wrapped
bvnWrp1 <- wrap_to_torus(bvn1)
# bvn2
p2 <- contour(as_tibble(bvnWrp1), "phi","psi", size=psize, alpha=palpha)
p2
ggsave("wrapped-bivariate.pdf", p2, width = 6, height = 6)


### 3. first simulating a point from the marginal distribution
# of one of the random variables and then simulating from the second
# random variable conditioned on the first.

rbvn<-function (n, m1, s1, m2, s2, rho) {
  X1 <- rnorm(n, mu1, s1)
  X2 <- rnorm(n, mu2 + (s2/s1) * rho *
                (X1 - mu1), sqrt((1 - rho^2)*s2^2))
  cbind(X1, X2)
}

# transfer them into the expected parameters
# expected covariance matrix s1^2, s1*s2*rho, s1*s2*rho, s2^2
const = 0.5 / (alpha[1]*alpha[2]-alpha[3]*alpha[3])
rho = -alpha[3] * const / ( sqrt(const * alpha[2] ) * sqrt(const * alpha[1] ) )
s1e = sqrt(const * alpha[2] ) * s1
s2e = sqrt(const * alpha[1] ) * s2

bvn3 <- rbvn(N,mu1,s1e,mu2,s2e,rho)
colnames(bvn3) <- c("phi","psi")
p1 <- contour(as_tibble(bvn3), "phi","psi", size=psize, alpha=palpha) +
  xlim(-5, 10) + ylim(-7.5, 12.5)
p1
ggsave("bivariate2.pdf", p1, width = 6, height = 6)

# wrapped
bvnWrp3 <- wrap_to_torus(bvn3)
# transfer to -pi pi
bvnWrp3 <- bvnWrp3 - pi
# bvn2
p2 <- contour(as_tibble(bvnWrp3), "phi","psi", size=psize, alpha=palpha)
p2
ggsave("wrapped-bivariate2.pdf", p2, width = 6, height = 6)

