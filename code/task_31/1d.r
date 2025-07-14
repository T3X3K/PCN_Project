library(ggplot2) |> suppressMessages()
library(igraph) |> suppressMessages()
library(dplyr) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()


schelling_constrained_1D <- function(L = 1000, r0 = 0.1, m = 0.2, max_iters = 1e5) {

  # population sizes
  n_vac <- round(r0 * L)
  n_pos <- round((1 - r0) * (1 + m)/2 * L)
  n_neg <- round((1 - r0) * (1 - m)/2 * L)
  total_filled <- n_vac + n_pos + n_neg
  if (total_filled != L) {
    n_neg <- n_neg + (L - total_filled)
  }

  # initial random state
  state <- sample(c(rep(0, n_vac), rep(1, n_pos), rep(-1, n_neg)))

  # neighbors
  get_neighbors <- function(i) {
    nb <- c(i - 1, i + 1)
    nb[nb >= 1 & nb <= L]
  }

  # unhappy
  unhappy <- function(i, state) {
    if (state[i] == 0) return(FALSE)
    nb <- get_neighbors(i)
    occ <- nb[state[nb] != 0]
    if (length(occ) == 0) return(FALSE)
    same <- sum(state[occ] == state[i])
    frac <- same / length(occ)
    return(frac < 0.5)
  }

  # initial unhappy
  unhappy_agents <- which(sapply(1:L, function(i) unhappy(i, state)))
  vacancies <- which(state == 0)

  it <- 0
  while (length(unhappy_agents) > 0 && it < max_iters) {
    it <- it + 1

    agent <- sample(unhappy_agents, 1)
    vac <- sample(vacancies, 1)

    # test move
    nb <- get_neighbors(vac)
    occ <- nb[state[nb] != 0]
    if (length(occ) == 0) {
      happy <- TRUE
    } else {
      same <- sum(state[occ] == state[agent])
      happy <- (same / length(occ)) >= 0.5
    }

    if (happy) {
      # move
      state[vac] <- state[agent]
      state[agent] <- 0
      vacancies[vacancies == vac] <- agent
      moved <- TRUE
    } else {
      moved <- FALSE
    }

    if (moved) {
      # update affected local areas
      affected <- unique(c(agent, vac, get_neighbors(agent), get_neighbors(vac)))
      affected <- affected[affected >= 1 & affected <= L]
      unhappy_flags <- sapply(affected, function(i) unhappy(i, state))
      unhappy_agents <- unique(c(
        setdiff(unhappy_agents, agent),      # remove old agent
        affected[unhappy_flags]              # add new unhappies
      ))
    } else {
      # agent did not move, keep it in the unhappy list
      # leave unhappy_agents unchanged
    }
  }

  # measure unhappiness of minority
  minority <- ifelse(n_pos < n_neg, 1, -1)
  minority_sites <- which(state == minority)
  if (length(minority_sites) == 0) {
    unhappy_final <- 0
  } else {
    unhappy_flags <- sapply(minority_sites, function(i) unhappy(i, state))
    unhappy_final <- sum(unhappy_flags) / length(minority_sites)
  }

  return(unhappy_final)
}



r0_seq <- 10^seq(-3, 0, length=20)
u_values <- sapply(r0_seq, function(r0) {
  mean(replicate(20, schelling_constrained_1D(L=100000, r0=r0, m=0.2)))
})
# print(paste("r0:", r0_seq, "u_values:", u_values))


r0_seq_04 <- 10^seq(-3, 0, length=20)
u_values_04 <- sapply(r0_seq_04, function(r0) {
  mean(replicate(20, schelling_constrained_1D(L=100000, r0=r0, m=0.4)))
})
# print(paste("r0:", r0_seq_04, "u_values:", u_values_04))

png("segregation_plot_1D_1e5.png")
# Combine all y values to set ylim automatically
all_y <- c(u_values, u_values_04)
plot(r0_seq, u_values, type="b", col="blue", log="x", xaxt="n", xlab="r0", ylab="Unhappiness (u)", 
	 main="Unhappiness vs r0", ylim=range(all_y, na.rm=TRUE))
lines(r0_seq_04, u_values_04, type="b", col="red", pch=2)
ticks <- 10^seq(floor(log10(min(c(r0_seq, r0_seq_04)))), ceiling(log10(max(c(r0_seq, r0_seq_04)))), by=1)
axis(1, at=ticks, labels=parse(text=paste0("10^", log10(ticks))))
legend("topright", legend=c("m=0.2", "m=0.4"), col=c("blue", "red"), pch=c(1,2), lty=1)
dev.off()

# Clean up variables from previous section to free memory
rm(r0_seq, u_values, r0_seq_04, u_values_04, all_y, ticks)
gc()
