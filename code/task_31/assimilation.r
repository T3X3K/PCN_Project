library(ggplot2) |> suppressMessages()
library(igraph) |> suppressMessages()
library(dplyr) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()

schelling_constrained_assimilation_1D <- function(L = 1000, r0 = 0.1, m = 0.2,
                                                  p_assimilation = 0.05,
                                                  max_iters = 10000) {
  
  
  n_vac <- round(r0 * L)
  n_pos <- round((1 - r0) * (1 + m)/2 * L)
  n_neg <- round((1 - r0) * (1 - m)/2 * L)
  total_filled <- n_vac + n_pos + n_neg
  if (total_filled != L) {
    n_neg <- n_neg + (L - total_filled)
  }
  
  
  state <- sample(c(rep(0, n_vac), rep(1, n_pos), rep(-1, n_neg)))
  
  get_neighbors <- function(i) {
    nb <- c(i - 1, i + 1)
    nb[nb >= 1 & nb <= L]
  }

  unhappy <- function(i, state) {
    if (state[i] == 0) return(FALSE)
    nb <- get_neighbors(i)
    occ <- nb[state[nb] != 0]
    if (length(occ) == 0) return(FALSE)
    same <- sum(state[occ] == state[i])
    frac <- same / length(occ)
    return(frac < 0.5)
  }

  unhappy_agents <- which(sapply(1:L, function(i) unhappy(i, state)))
  vacancies <- which(state == 0)

  max_iters <- 500  
  no_change_steps <- 0
  last_unhappy_count <- Inf
  
  iter <- 0
  while (length(unhappy_agents) > 0 && iter < max_iters && no_change_steps < 10) {
    iter <- iter + 1

    
    agent <- sample(unhappy_agents, 1)
    vac <- sample(vacancies, 1)
    
    nb <- get_neighbors(vac)
    occ <- nb[state[nb] != 0]
    happy <- if (length(occ) == 0) TRUE else {
      same <- sum(state[occ] == state[agent])
      (same / length(occ)) >= 0.5
    }

    if (happy) {
      state[vac] <- state[agent]
      state[agent] <- 0
      vacancies[vacancies == vac] <- agent
    }

    
    unhappy_agents <- which(sapply(1:L, function(i) unhappy(i, state)))
    
    
    for (i in 1:L) {
      if (state[i] == 0) next
      nb <- get_neighbors(i)
      occ <- nb[state[nb] != 0]
      if (length(occ) == 0) next
      opposite <- sum(state[occ] == -state[i])
      if (opposite > length(occ)/2 && runif(1) < p_assimilation) {
        state[i] <- -state[i]
      }
    }

    current_unhappy_count <- length(unhappy_agents)
    if (current_unhappy_count == last_unhappy_count) {
      no_change_steps <- no_change_steps + 1
    } else {
      no_change_steps <- 0
    }
    last_unhappy_count <- current_unhappy_count

  }
  
  
  minority <- ifelse(n_pos < n_neg, 1, -1)
  majority <- -minority
  
  minority_sites <- which(state == minority)
  majority_sites <- which(state == majority)
  
  u_min <- if (length(minority_sites) == 0) 0 else
    sum(sapply(minority_sites, function(i) unhappy(i, state))) / length(minority_sites)
  
  u_maj <- if (length(majority_sites) == 0) 0 else
    sum(sapply(majority_sites, function(i) unhappy(i, state))) / length(majority_sites)
  
  return(c(minority = u_min, majority = u_maj))
}


p_assim_vals <- c(0, 0.00005, 0.0005, 0.001, 0.01)
r0_seq <- 10^seq(-3, 0, length.out = 20)
num_reps <- 20

results <- data.frame()

for (p_assim in p_assim_vals) {
  for (r0 in r0_seq) {
    reps <- replicate(num_reps, {
      res <- schelling_constrained_assimilation_1D(
        L = 10000,
        r0 = r0,
        m = 0.2,
        p_assimilation = p_assim
      )
      res["minority"]  
    })
    avg <- mean(reps)
    results <- rbind(results, data.frame(
      r0 = r0,
      p_assim = p_assim,
      u_min = avg
    ))
  }
}

write.csv(results, file = "assimilation.csv", row.names = FALSE)

ggplot(results, aes(x = r0, y = u_min, color = factor(p_assim))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Minority Unhappiness vs Vacancy Rate (r₀)",
    x = expression(r[0]),
    y = "Unhappiness (Minority)",
    color = "Assimilation rate"
  ) +
  theme_minimal()
ggsave("segregation_plot_assimilation.png")
