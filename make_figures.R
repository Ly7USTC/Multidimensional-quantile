suppressPackageStartupMessages({library(ggplot2); library(MASS); library(matrixStats); library(Hmisc)})
set.seed(2019)
sv <- function(p, f) ggsave(file.path("figures", f), p, width = 5, height = 5, dpi = 150)

source("mqnorm.R"); source("mqemp.R"); source("mqdis.R")

Sig <- matrix(c(1, 0.5, 0.5, 1), 2)
sv(mqnorm(0.9, Sigma = Sig, Xlim = c(-4, 5), Ylim = c(-5, 5)) +
     ggtitle("Normal, q = 0.9, cone R²₊"), "normal_q09_R2plus.png")

C <- matrix(c(0, 1, 1, -1), byrow = TRUE, nrow = 2)
sv(mqnorm(0.9, Sigma = Sig, C = C, Xlim = c(-4, 5), Ylim = c(-5, 5)) +
     ggtitle("Normal, q = 0.9, cone spanned by (0,1), (1,-1)"), "normal_q09_cone.png")

D <- matrix(c(0, 0,  1, 2,  2, -1,  -1, 1,  3, 1), ncol = 2, byrow = TRUE)
sv(mqdis(0.5, D, we = c(0.3, 0.25, 0.2, 0.15, 0.1)) +
     ggtitle("Weighted discrete, q = 0.5"), "discrete_q05.png")

mu_1 <- c(0, 0);  Sig_1 <- diag(2)
mu_2 <- c(2, 1);  Sig_2 <- matrix(c(1, -0.5, -0.5, 1), 2)
source("CxLSprop.R")
sv(CxLSnorm(0.95, mu_1, Sig_1, mu_2, Sig_2, Lam = 0.5) +
     ggtitle("Mixture (black) vs components, q = 0.95"), "cxls_mixture_q095.png")
