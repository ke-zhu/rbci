library(tidyverse)
library(parallel)
library(tictoc)
source("fun.R")
set.seed(2024)

case <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
n_rep <- 1000

n_seq <- c(50, 100, 500, 1000, 5000, 10000)
n_fisher_seq <- c(10^4, 2*10^4)

setup <- expand_grid(n = n_seq, n_fisher = n_fisher_seq) %>% 
  mutate(case = row_number(), .before = everything())

list2env(setup[case,], envir = .GlobalEnv)

# data
n1 <- n / 2
Y0 <- rnorm(n)
tau <- 1
Y1 <- Y0 + tau
z <- rep(0, n)
z[sample(n, n1)] <- 1
Y <- z * Y1 + (1 - z) * Y0

# reference assignments
z_set <- map(1:n_fisher, ~ {
  zp <- rep(0, n)
  zp[sample(n, n1)] <- 1
  zp
})

writeLines("", str_glue("output/prog_{n}_{n_fisher}.txt"))
# inference
inf_res <- mclapply(1:n_rep, function(rep) {
  z <- rep(0, n)
  z[sample(n, n1)] <- 1
  Y <- z * Y1 + (1 - z) * Y0
  # FRT
  tic()
  p_value <- pvalue(0, z, Y0, z_set, "!=")
  t_p <- toc()
  # RBCI
  tic()
  rbci_res <- rbci(z, Y, z_set)
  t_ci <- toc()
  # output
  cat(rep, "\n", file = str_glue("output/prog_{n}_{n_fisher}.txt"), append = TRUE)
  tibble(
    rbci_res,
    p_value,
    time_p = t_p$toc - t_p$tic,
    time_ci = t_ci$toc - t_ci$tic
  )
}, mc.cores = n_cores) %>% map_dfr(~.)

# summary
sim_res <- inf_res %>% 
  summarise(
    CP = mean((ci_l <= tau) & (ci_u >= tau)), 
    `Type I` = mean(p_value <= 0.05),
    time_p = mean(time_p),
    time_ci = mean(time_ci),
  ) %>% 
  mutate(n, n_fisher, .before = everything())


save(sim_res, file = str_glue("output/metrics_supp_{case}.RData"))
