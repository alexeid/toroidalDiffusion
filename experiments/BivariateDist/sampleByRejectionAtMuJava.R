# check if Java implementation is same as R

setwd("~/WorkSpace/toroidalDiffusion/experiments/BivariateDist")

library(tidyverse)

# contour plot for bivariate normal data
contour <- function(bvn, x_string, y_string, size, alpha){
  require("ggplot2")
  p <- ggplot(bvn, aes_string(x = x_string, y = y_string)) +
    geom_point(size=size, alpha = alpha) +
    #geom_density2d() +
    labs(title=paste0("mu1 = ", round(mu1,2), ", mu2 = ",
                      round(mu2,2), ", sampling ", N, " points")) +
    theme_bw()
  return(p)
}

psize=.1; palpha = .1
mu1 <- pi * 0.65; s1 <- 1.5
mu2 <- pi * 0.8; s2 <- 1.5
alpha <- c(1.0, 0.5, 0.5)

angles <- read_delim("sampleByRejectionAtMu.txt", delim = "\t");
range(angles[,1:2])

N <- nrow(angles)
p1 <- contour(angles[,1:2], "phi","psi", size=psize, alpha=palpha)
p1
ggsave("sampleByRejectionAtMuJava.pdf", p1, width = 6, height = 6)


angles2 <- read_delim("WrappedBivariateNormalJava.txt", delim = "\t");
range(angles2[,1:2])

N <- nrow(angles2)
p2 <- contour(angles2[,1:2], "phi","psi", size=psize, alpha=palpha)
p2
ggsave("WrappedBivariateNormalJava.pdf", p2, width = 6, height = 6)




