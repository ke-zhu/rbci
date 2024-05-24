# solve theta -------------------------------------------------------------

diff_stat <- function(z, Y) {
  mean(Y[z == 1]) - mean(Y[z == 0])
}

s2_stat <- function(z, Y, d) {
  cov(Y[z == 1], d[z == 1]) / sum(z == 1) +
    cov(Y[z == 0], d[z == 0]) / sum(z == 0)
}

t_stat <- function(z, Y) {
  diff_stat(z, Y) / sqrt(s2_stat(z, Y, Y))
}

solve_quad <- function(a, b, c) {
  d2 <- b^2 - 4 * a * c
  if (d2 < 0) {
    NULL
  } else if (d2 == 0) {
    -b / (2 * a)
  } else {
    (-b + c(-1, 1) * sqrt(d2)) / (2 * a)
  }
}

solve_theta <- function(z, zp, Y) {
  if (identical(z, zp)) {
    return(NULL)
  } else {
    d <- zp - z
    t_obs <- t_stat(z, Y)
    a <- t_obs^2 * s2_stat(zp, d, d) - diff_stat(zp, d)^2
    b <- 2 * (t_obs^2 * s2_stat(zp, Y, d) - diff_stat(zp, Y) * diff_stat(zp, d))
    c <- t_obs^2 * s2_stat(zp, Y, Y) - diff_stat(zp, Y)^2
    root <- solve_quad(a, b, c)
    if (length(root) == 2) {
      e <- abs(root[2] - root[1]) / 2
    } else {
      e <- 1
    }
    jump <- map_int(root, function(theta) {
      t_l <- t_stat(zp, Y + (zp - z) * (theta - e))
      t_r <- t_stat(zp, Y + (zp - z) * (theta + e))
      (sign(t_r - t_obs) - sign(t_l - t_obs)) / 2
    })
    tibble(theta = root, jump) %>% filter(jump != 0)
  }
}


# rbci --------------------------------------------------------------

pvalue <- function(theta, z, Y, z_set, H1 = "!=") {
  stat_obs <- t_stat(z, Y)
  stat_rep <- sapply(z_set, function(zp) {
    d <- zp - z
    Y_tilde <- Y + d * theta
    t_stat(zp, Y_tilde)
  })
  if (H1 == ">") {
    mean(stat_rep >= stat_obs)
  } else if (H1 == "<") {
    mean(stat_rep <= stat_obs)
  } else if (H1 == "!=") {
    2 * min(
      mean(stat_rep >= stat_obs),
      mean(stat_rep <= stat_obs)
    )
  }
}

plot_p_fun <- function(z, Y, z_set, H1 = ">") {
  theta_jump_set <- map_dfr(z_set, ~ solve_theta(z, ., Y)) %>% 
    arrange(theta) %>% 
    mutate(theta = round(theta, 10)) %>% 
    group_by(theta) %>% 
    summarise(jump = sum(jump)) # merge same theta
  tt <- theta_jump_set$theta
  J <- theta_jump_set$jump
  h <- 1 / length(z_set)
  if (H1 == ">") {
    # construct p value function (">")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, ">")
    p_set <- c(p0, p0 + cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    plot(pfun, verticals = F, do.points = F, 
         xlab = TeX(r'($\theta$)'),
         ylab = TeX(r'($p^+(\theta)$)'),
         main = TeX(r'($H_1^{\theta+}: Y_i(1) - Y_i(0) > \theta$)'))
  } else if (H1 == "<") {
    # construct p value function ("<")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, "<")
    p_set <- c(p0, p0 - cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    plot(pfun, verticals = F, do.points = F, 
         xlab = TeX(r'($\theta$)'),
         ylab = TeX(r'($p^-(\theta)$)'),
         main = TeX(r'($H_1^{\theta-}:$ $Y_i(1) - Y_i(0) < \theta$)'))
  }
}

rbci <- function(z, Y, z_set, sig_level = 0.05, H1 = "!=") {
  theta_jump_set <- map_dfr(z_set, ~ solve_theta(z, ., Y)) %>% 
    arrange(theta) %>% 
    mutate(theta = round(theta, 10)) %>% 
    group_by(theta) %>% 
    summarise(jump = sum(jump)) # merge same theta
  tt <- theta_jump_set$theta
  J <- theta_jump_set$jump
  h <- 1 / length(z_set)
  ci_u <- Inf
  ci_l <- -Inf
  if (H1 == ">") {
    # construct p value function (">")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, ">")
    p_set <- c(p0, p0 + cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    # squeeze 
    tt_plus <- c(head(tt, 1) - 1, tt)
    i <- 1
    while (pfun(tt_plus[i]) <= sig_level) {
      if ((i + 1) <= length(tt_plus)) {
        ci_l <- tt_plus[i + 1]
        i <- i + 1
      } else {
        ci_l <- Inf
        break
      }
    }
  } else if (H1 == "<") {
    # construct p value function ("<")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, "<")
    p_set <- c(p0, p0 - cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    # squeeze
    tt_plus <- c(tail(tt, 1) + 1, rev(tt))
    i <- 1
    while (pfun(tt_plus[i]) <= sig_level) {
      if ((i + 1) <= length(tt_plus)) {
        ci_u <- tt_plus[i + 1]
        i <- i + 1
      } else {
        ci_u <- -Inf
        break
      }
    }
  } else if (H1 == "!=") {
    # construct p value function (">")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, ">")
    p_set <- c(p0, p0 + cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    # squeeze 
    tt_plus <- c(head(tt, 1) - 1, tt)
    i <- 1
    while (pfun(tt_plus[i]) <= sig_level / 2) {
      if ((i + 1) <= length(tt_plus)) {
        ci_l <- tt_plus[i + 1]
        i <- i + 1
      } else {
        ci_l <- Inf
        break
      }
    }
    # construct p value function ("<")
    p0 <- pvalue(tt[1] - 1, z, Y, z_set, "<")
    p_set <- c(p0, p0 - cumsum(J) * h)
    pfun <- stepfun(tt, p_set)
    # squeeze
    tt_plus <- c(tail(tt, 1) + 1, rev(tt))
    i <- 1
    while (pfun(tt_plus[i]) <= sig_level / 2) {
      if ((i + 1) <= length(tt_plus)) {
        ci_u <- tt_plus[i + 1]
        i <- i + 1
      } else {
        ci_u <- -Inf
        break
      }
    }
  }
  tibble(ci_l, ci_u)
}
