library(ggplot2) |> suppressMessages()
library(igraph) |> suppressMessages()
library(dplyr) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()

r0_seq <- 10^seq(-3, 0, length.out = 20)  
leave_vals <- c(0.0, 0.01, 0.05)
enter_vals <- c(0.0, 0.01, 0.05)
num_reps <- 20  
L <- 10000
m <- 0.2
max_steps <- 10000

schelling_constrained_migration_1D <- function(L = 1000, r0 = 0.1, m = 0.2, 
                                               p_leave_majority = 0.01, p_enter_minority = 0.01,
                                               max_total_steps = 10000) {
  
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
  
  
  for (step in 1:max_total_steps) {
    
    if (length(unhappy_agents) == 0) break  

    agent <- sample(unhappy_agents, 1)
    vac <- sample(vacancies, 1)
    
    
    nb <- get_neighbors(vac)
    occ <- nb[state[nb] != 0]
    if (length(occ) == 0) {
      happy <- TRUE
    } else {
      same <- sum(state[occ] == state[agent])
      happy <- same / length(occ) >= 0.5
    }
    
    if (happy) {
      
      state[vac] <- state[agent]
      state[agent] <- 0
      vacancies[vacancies == vac] <- agent
    }
    
    
    unhappy_agents <- which(sapply(1:L, function(i) unhappy(i, state)))

    
    majority_indices <- which(state == 1)
    leavers <- majority_indices[runif(length(majority_indices)) < p_leave_majority]
    state[leavers] <- 0
    vacancies <- union(vacancies, leavers)  

    
    vacancy_indices <- which(state == 0)
    arrivals <- vacancy_indices[runif(length(vacancy_indices)) < p_enter_minority]
    state[arrivals] <- -1

    
    if (length(unhappy_agents) == 0) {
      break  
    }
  }
  
  
  minority_sites <- which(state == -1)
  majority_sites <- which(state == 1)

  unhappy_min <- if (length(minority_sites) == 0) 0 else
    sum(sapply(minority_sites, function(i) unhappy(i, state))) / length(minority_sites)

  unhappy_maj <- if (length(majority_sites) == 0) 0 else
    sum(sapply(majority_sites, function(i) unhappy(i, state))) / length(majority_sites)

  return(c(minority = unhappy_min, majority = unhappy_maj))
}



results <- data.frame()


for (leave in leave_vals) {
  for (enter in enter_vals) {
    for (r0 in r0_seq) {
      reps <- replicate(num_reps, {
        res <- schelling_constrained_migration_1D(
          L = L, r0 = r0, m = m,
          p_leave_majority = leave,
          p_enter_minority = enter,
          max_total_steps = max_steps
        )
        res["minority"]  
      })
      avg <- mean(reps)
      results <- rbind(results, data.frame(
        r0 = r0,
        u_min = avg,
        p_leave = leave,
        p_enter = enter,
        label = paste0("leave=", leave, ", enter=", enter)
      ))
    }
  }
}

results$group <- interaction(results$p_leave, results$p_enter, sep = ", ")

ggplot(results, aes(x = r0, y = u_min, 
                    color = group, 
                    linetype = group, 
                    shape = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Minority Unhappiness vs Vacancy Rate (r₀)",
    x = expression(r[0]),
    y = "Unhappiness (Minority)",
    color = "leave, enter",
    shape = "leave, enter",
    linetype = "leave, enter"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

write.csv(results, file = "migration_minority.csv", row.names = FALSE)
ggsave("segregation_plot_migration_minority.png")

rm(leave, enter, r0, reps, avg, results)
gc()
