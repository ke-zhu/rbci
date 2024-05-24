library(tidyverse)
library(parallel)
library(tictoc)
library(latex2exp)
source("fun.R")
set.seed(2024)


# 1 example in Luo et al. (2021) -----------------------------------------------------------------

# 1.1 data --------

n <- 8
n1 <- 4
Y0 <- c(0.14, 1.12, 0.80, 1.80, 0.90, 0.44, 1.13, 0.53)
Y1 <- Y0 + 1
z <- c(1, 1, 0, 1, 0, 0, 1, 0)
Y <- z * Y1 + (1 - z) * Y0
z_set <- map(1:choose(n, n1), ~ {
  zp <- rep(0, n)
  zp[combn(n, n1)[, .x]] <- 1
  zp
})


# 1.2 solve theta for one z_pi --------

zp <- z_set[[1]]
d <- zp - z
a <- t_stat(z, Y)^2 * s2_stat(zp, d, d) - diff_stat(zp, d)^2
b <- 2 * (t_stat(z, Y)^2 * s2_stat(zp, Y, d) - diff_stat(zp, Y) * diff_stat(zp, d))
c <- t_stat(z, Y)^2 * s2_stat(zp, Y, Y) - diff_stat(zp, Y)^2
root <- solve_quad(a, b, c)
if(!is.null(root)) {
  xmin <- 2 * min(root) - max(root)
  xmax <- 2 * max(root) - min(root)
}
x <- seq(xmin, xmax, length.out = 1000)
# plot(x, a * x^2 + b * x + c, type = "l")
# abline(h = 0)
y <- map_dbl(x, function(theta) {
  t_stat(zp, Y + (zp - z) * theta)
})
pdf(file = "output/theta_t.pdf", width = 10, height = 6)
plot(x, y, type = "l", 
     xlab = TeX(r'($\theta$)'),
     ylab = TeX(r'($t(z_{\pi},Y+\delta\theta)$)'),
     ylim = c(1.5, 4.5)
)
abline(h = t_stat(z, Y), lty = 2)
abline(v = root, lty = 2)
text(0.5, t_stat(z, Y) + 0.2, TeX(r'($t(z,Y)=3.85$)'))
text(root[1] + 0.3, 2, TeX(r'($\theta_{\pi,1}=1.33, J_{\pi,1}=1$)'))
text(root[2] + 0.3, 2, TeX(r'($\theta_{\pi,2}=2.27, J_{\pi,2}=-1$)'))
dev.off()


# 1.3 solve theta for all z_pi --------

theta_jump_set <- map_dfr(z_set, ~ solve_theta(z, ., Y)) %>% 
  arrange(theta) %>% 
  mutate(theta = round(theta, 10)) %>% 
  group_by(theta) %>% 
  summarise(jump = sum(jump)) # merge same theta
theta_set <- theta_jump_set$theta


# 1.4 construct p value function (">") --------

p0 <- pvalue(theta_set[1] - 1, z, Y, z_set, ">")
h <- 1 / length(z_set)
p_set <- c(p0, p0 + cumsum(theta_jump_set$jump) * h)
pfun <- stepfun(theta_set, p_set)

# check 1
tibble(
  x = c(theta_set[1] - 1, theta_set + 1e-5),
  p0 = map_dbl(x, ~ pvalue(., z, Y, z_set, ">")),
  p = p_set,
  d = abs(p0 - p),
  check = d < 1e-15
) %>% pull(check) %>% all

# check 2
tibble(
  x = c(
    theta_set[1] - 1, theta_set + 1e-5,
    runif(50, min(theta_set), max(theta_set))
  ),
  p0 = map_dbl(x, ~ pvalue(., z, Y, z_set, ">")),
  p = pfun(x),
  d = abs(p0 - p),
  check = d < 1e-15
) %>% pull(check) %>% all

# plot p value function
pdf(file = "output/pfun1.pdf", width = 10, height = 6)
plot(pfun, verticals = F, do.points = F, 
     xlab = TeX(r'($\theta$)'),
     ylab = TeX(r'($p^+(\theta)$)'),
     main = TeX(r'($H_1^{\theta+}:$ $Y_i(1) - Y_i(0) > \theta$)'))
abline(v = rbci(z, Y, z_set, H1 = ">")$ci_l, lty = 2)
text(rbci(z, Y, z_set, H1 = ">")$ci_l + 0.25, 0.3, TeX(r'($c_l=0.61$)'))
# x <- runif(100, min(theta_set) - 1, max(theta_set) + 1)
# points(x, map_dbl(x, ~ pvalue(., z, Y, z_set, ">")))
dev.off()


# 1.5 construct p value function ("<") --------

p0 <- pvalue(theta_set[1] - 1, z, Y, z_set, "<")
h <- 1 / length(z_set)
p_set <- c(p0, p0 - cumsum(theta_jump_set$jump) * h)
pfun <- stepfun(theta_set, p_set)

# check 1
tibble(
  x = c(theta_set[1] - 1, theta_set + 1e-5),
  p0 = map_dbl(x, ~ pvalue(., z, Y, z_set, "<")),
  p = p_set,
  d = abs(p0 - p),
  check = d < 1e-15
) %>% pull(check) %>% all

# check 2
tibble(
  x = c(
    theta_set[1] - 1, theta_set + 1e-5,
    runif(50, min(theta_set), max(theta_set))
  ),
  p0 = map_dbl(x, ~ pvalue(., z, Y, z_set, "<")),
  p = pfun(x),
  d = abs(p0 - p),
  check = d < 1e-15
) %>% pull(check) %>% all

# plot p value function
pdf(file = "output/pfun2.pdf", width = 10, height = 6)
plot(pfun, verticals = F, do.points = F, 
     xlab = TeX(r'($\theta$)'),
     ylab = TeX(r'($p^-(\theta)$)'),
     main = TeX(r'($H_1^{\theta-}:$ $Y_i(1) - Y_i(0) < \theta$)'))
# x <- c(runif(100, min(theta_set) - 1, max(theta_set) + 1))
# points(x, map_dbl(x, ~ pvalue(., z, Y, z_set, "<")))
dev.off()


# 1.6 pvalue & rbci --------

tic()
pvalue(0, z, Y, z_set, "!=")
toc()

tic()
rbci(z, Y, z_set)
toc()


# 1.7 repeated sampling performance --------

tic()
inf_res_small_n <- mclapply(z_set, function(z) {
  Y <- z * Y1 + (1 - z) * Y0
  tibble(rbci(z, Y, z_set), p_value = pvalue(0, z, Y0, z_set, H1 = "!="))
}, mc.cores = 4) %>% map_dfr(~.)
runtime_small_n <- toc()

(metrics_small_n <- inf_res_small_n %>% 
    summarise(
      CP = mean((ci_l <= tau) & (ci_u >= tau)), 
      `Type I` = mean(p_value <= 0.05)
    ))




# 2 simulation for n = 100 -----------------------------------------

# 2.1 data --------

RNGkind("L'Ecuyer-CMRG")
set.seed(2024)

n_fisher <- 10^4
n_rep <- 1000
n <- 100
n1 <- n / 2
Y0 <- rnorm(n)
tau <- 1
Y1 <- Y0 + tau
z <- rep(0, n)
z[sample(n, n1)] <- 1
Y <- z * Y1 + (1 - z) * Y0
if (choose(n, n1) <= n_fisher) {
  z_set <- map(1:choose(n, n1), ~ {
    zp <- rep(0, n)
    zp[combn(n, n1)[, .x]] <- 1
    zp
  })
} else {
  z_set <- map(1:n_fisher, ~ {
    zp <- rep(0, n)
    zp[sample(n, n1)] <- 1
    zp
  })
}


# 2.2 pvalue & rbci --------

plot_p_fun(z, Y, z_set, H1 = ">")
plot_p_fun(z, Y, z_set, H1 = "<")

tic()
pvalue(0, z, Y, z_set, "!=")
toc()

tic()
rbci(z, Y, z_set)
toc()


# 2.3 grid method --------

n_grid <- 100
init <- t.test(Y[z==1], Y[z==0], conf.level = 0.999)$conf.int
grid <- seq(init[1], init[2], length.out = n_grid)
tic()
p_grid <- map_dbl(grid, ~ {
  pvalue(., z, Y, z_set, ">")
})
runtime_grid <- toc()

(runtime_grid$toc - runtime_grid$tic) / n_grid * 2 * n_fisher / 3600

pfun <- stepfun(grid[-1], p_grid)
plot(pfun, verticals = F, do.points = F, 
     xlab = TeX(r'($\theta$)'),
     ylab = TeX(r'($p^+(\theta)$)'),
     main = TeX(r'($H_1^{\theta+}: Y_i(1) - Y_i(0) > \theta$)'))



# 2.4 repeated sampling performance --------

tic()
inf_res_large_n <- mclapply(1:n_rep, function(rep) {
  z <- rep(0, n)
  z[sample(n, n1)] <- 1
  Y <- z * Y1 + (1 - z) * Y0
  tibble(rbci(z, Y, z_set), p_value = pvalue(0, z, Y0, z_set, H1 = "!="))
}, mc.cores = 4) %>% map_dfr(~.)
runtime_large_n <- toc()


(metrics_large_n <- inf_res_large_n %>% 
    summarise(
      CP = mean((ci_l <= tau) & (ci_u >= tau)), 
      `Type I` = mean(p_value <= 0.05)
    ))


save(
  inf_res_small_n, runtime_small_n, metrics_small_n, 
  inf_res_large_n, runtime_large_n, metrics_large_n, 
  file = "output/metrics.RData"
)



