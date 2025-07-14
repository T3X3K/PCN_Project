library(ggplot2) |> suppressMessages()
library(igraph) |> suppressMessages()
library(dplyr) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()



henry_model_asymmetric <- function(n = 300, m = 2, p_in = 0.1, p_out = 0.01,
                                   minority_frac = 0.3,
                                   max_steps = 300) {
  block_sizes <- rep(n / m, m)
  pref_matrix <- matrix(p_out, nrow = m, ncol = m)
  diag(pref_matrix) <- p_in
  
  g <- sample_sbm(n, pref.matrix = pref_matrix, block.sizes = block_sizes)
  
  # Assign asymmetric types
  V(g)$type <- sample(c(1, -1), n, replace = TRUE, prob = c(1 - minority_frac, minority_frac))
  
  homophily <- numeric(max_steps)
  minority_stress <- numeric(max_steps)
  
  for (step in 1:max_steps) {
    types <- V(g)$type
    edge_ends <- ends(g, E(g))
    dissim_idx <- which(types[edge_ends[, 1]] != types[edge_ends[, 2]])
    dissim_edges <- E(g)[dissim_idx]
    
    if (length(dissim_edges) == 0) {
      homophily <- homophily[1:(step - 1)]
      minority_stress <- minority_stress[1:(step - 1)]
      break
    }
    
    # Rewire a dissimilar edge
    e <- sample(dissim_edges, 1)
    ends_e <- ends(g, e)
    u <- ends_e[1]
    v <- ends_e[2]
    u_type <- V(g)$type[u]
    
    candidates <- which(V(g)$type == u_type &
                        !are.connected(g, u, 1:n) &
                        1:n != u)
    
    if (length(candidates) > 0) {
      new_v <- sample(candidates, 1)
      g <- delete_edges(g, e)
      g <- add_edges(g, c(u, new_v))
    }
    
    # Track homophily
    edge_ends <- ends(g, E(g))
    same <- sum(types[edge_ends[, 1]] == types[edge_ends[, 2]])
    homophily[step] <- same / ecount(g)
    
    # Track minority stress
    minority_nodes <- which(types == -1)
    stress_vals <- sapply(minority_nodes, function(v) {
      nbs <- neighbors(g, v)
      if (length(nbs) == 0) return(0)
      mean(types[nbs] != -1)
    })
    minority_stress[step] <- mean(stress_vals)
  }
  
  return(list(graph = g,
              homophily = homophily,
              stress = minority_stress,
              types = V(g)$type))
}

set.seed(42)

result <- henry_model_asymmetric(n = 400, minority_frac = 0.3, max_steps = 1000)

plot(result$homophily, type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1), xlab = "Step", ylab = "Value",
     main = "Segregation Dynamics Over Time (SBM)")
lines(result$stress, col = "steelblue", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("Homophily", "Minority Stress"),
       col = c("darkred", "steelblue"),
       lty = c(1, 2), lwd = 2, bty = "n")


library(igraph)

henry_model_asymmetric_2D <- function(n = 256, minority_frac = 0.3, max_steps = 300) {
  side_length <- sqrt(n)
  if (side_length != floor(side_length)) {
    stop("n must be a perfect square to form a square lattice.")
  }

  # Create a 2D lattice graph
  g <- erdos.renyi.game(n, p.or.m = 0.02, type = "gnp", directed = FALSE)



  # Get actual number of vertices in g
  node_count <- vcount(g)

  # Assign types to vertices based on actual graph size
  V(g)$type <- sample(c(1, -1), node_count, replace = TRUE,
                      prob = c(1 - minority_frac, minority_frac))

  homophily <- numeric(max_steps)
  minority_stress <- numeric(max_steps)

  for (step in 1:max_steps) {
    types <- V(g)$type
    edge_ends <- ends(g, E(g))
    dissim_idx <- which(types[edge_ends[, 1]] != types[edge_ends[, 2]])
    dissim_edges <- E(g)[dissim_idx]

    if (length(dissim_edges) == 0) {
  if (step > 1) {
    homophily <- homophily[1:(step - 1)]
    minority_stress <- minority_stress[1:(step - 1)]
  } else {
    homophily <- c(NA)
    minority_stress <- c(NA)
  }
  break
}


    # Rewire a dissimilar edge
    e <- sample(dissim_edges, 1)
    ends_e <- ends(g, e)
    u <- as.integer(ends_e[1])
    v <- as.integer(ends_e[2])
    u_type <- V(g)$type[u]

    candidates <- which(V(g)$type == u_type &
                        !are.connected(g, u, V(g)) &
                        V(g) != u)

    if (length(candidates) > 0) {
      new_v <- sample(candidates, 1)
      g <- delete_edges(g, e)
      g <- add_edges(g, c(u, new_v))
    }

    # Track homophily
    edge_ends <- ends(g, E(g))
    same <- sum(types[edge_ends[, 1]] == types[edge_ends[, 2]])
    homophily[step] <- same / ecount(g)

    # Track minority stress
    minority_nodes <- which(types == -1)
    stress_vals <- sapply(minority_nodes, function(v) {
      nbs <- neighbors(g, v)
      if (length(nbs) == 0) return(0)
      mean(types[nbs] != -1)
    })
    minority_stress[step] <- mean(stress_vals)
  }

  return(list(graph = g,
              homophily = homophily,
              stress = minority_stress,
              types = V(g)$type))
}


set.seed(42)
result <- henry_model_asymmetric_2D(n = 400, minority_frac = 0.3, max_steps = 1000)


plot(result$homophily, type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1), xlab = "Step", ylab = "Value",
     main = "Segregation Dynamics (ER)")
lines(result$stress, col = "steelblue", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("Homophily", "Minority Stress"),
       col = c("darkred", "steelblue"),
       lty = c(1, 2), lwd = 2, bty = "n")

library(igraph)

henry_model_asymmetric_2D <- function(n = 256, minority_frac = 0.3, max_steps = 300) {
  side_length <- sqrt(n)
  if (side_length != floor(side_length)) {
    stop("n must be a perfect square to form a square lattice.")
  }

  # Create a 2D lattice graph
  g <- watts.strogatz.game(1, size = n, nei = 5, p = 0.05)



  # Get actual number of vertices in g
  node_count <- vcount(g)

  # Assign types to vertices based on actual graph size
  V(g)$type <- sample(c(1, -1), node_count, replace = TRUE,
                      prob = c(1 - minority_frac, minority_frac))

  homophily <- numeric(max_steps)
  minority_stress <- numeric(max_steps)

  for (step in 1:max_steps) {
    types <- V(g)$type
    edge_ends <- ends(g, E(g))
    dissim_idx <- which(types[edge_ends[, 1]] != types[edge_ends[, 2]])
    dissim_edges <- E(g)[dissim_idx]

    if (length(dissim_edges) == 0) {
  if (step > 1) {
    homophily <- homophily[1:(step - 1)]
    minority_stress <- minority_stress[1:(step - 1)]
  } else {
    homophily <- c(NA)
    minority_stress <- c(NA)
  }
  break
}


    # Rewire a dissimilar edge
    e <- sample(dissim_edges, 1)
    ends_e <- ends(g, e)
    u <- as.integer(ends_e[1])
    v <- as.integer(ends_e[2])
    u_type <- V(g)$type[u]

    candidates <- which(V(g)$type == u_type &
                        !are.connected(g, u, V(g)) &
                        V(g) != u)

    if (length(candidates) > 0) {
      new_v <- sample(candidates, 1)
      g <- delete_edges(g, e)
      g <- add_edges(g, c(u, new_v))
    }

    # Track homophily
    edge_ends <- ends(g, E(g))
    same <- sum(types[edge_ends[, 1]] == types[edge_ends[, 2]])
    homophily[step] <- same / ecount(g)

    # Track minority stress
    minority_nodes <- which(types == -1)
    stress_vals <- sapply(minority_nodes, function(v) {
      nbs <- neighbors(g, v)
      if (length(nbs) == 0) return(0)
      mean(types[nbs] != -1)
    })
    minority_stress[step] <- mean(stress_vals)
  }

  return(list(graph = g,
              homophily = homophily,
              stress = minority_stress,
              types = V(g)$type))
}

set.seed(42)
result <- henry_model_asymmetric_2D(n = 400, minority_frac = 0.3, max_steps = 1000)


plot(result$homophily, type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1), xlab = "Step", ylab = "Value",
     main = "Segregation Dynamics (WS)")
lines(result$stress, col = "steelblue", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("Homophily", "Minority Stress"),
       col = c("darkred", "steelblue"),
       lty = c(1, 2), lwd = 2, bty = "n")


library(igraph)

henry_model_asymmetric_2D <- function(n = 256, minority_frac = 0.3, max_steps = 300) {
  side_length <- sqrt(n)
  if (side_length != floor(side_length)) {
    stop("n must be a perfect square to form a square lattice.")
  }

  # Create a 2D lattice graph
  g <- barabasi.game(n, m = 3, directed = FALSE)




  # Get actual number of vertices in g
  node_count <- vcount(g)

  # Assign types to vertices based on actual graph size
  V(g)$type <- sample(c(1, -1), node_count, replace = TRUE,
                      prob = c(1 - minority_frac, minority_frac))

  homophily <- numeric(max_steps)
  minority_stress <- numeric(max_steps)

  for (step in 1:max_steps) {
    types <- V(g)$type
    edge_ends <- ends(g, E(g))
    dissim_idx <- which(types[edge_ends[, 1]] != types[edge_ends[, 2]])
    dissim_edges <- E(g)[dissim_idx]

    if (length(dissim_edges) == 0) {
  if (step > 1) {
    homophily <- homophily[1:(step - 1)]
    minority_stress <- minority_stress[1:(step - 1)]
  } else {
    homophily <- c(NA)
    minority_stress <- c(NA)
  }
  break
}


    # Rewire a dissimilar edge
    e <- sample(dissim_edges, 1)
    ends_e <- ends(g, e)
    u <- as.integer(ends_e[1])
    v <- as.integer(ends_e[2])
    u_type <- V(g)$type[u]

    candidates <- which(V(g)$type == u_type &
                        !are.connected(g, u, V(g)) &
                        V(g) != u)

    if (length(candidates) > 0) {
      new_v <- sample(candidates, 1)
      g <- delete_edges(g, e)
      g <- add_edges(g, c(u, new_v))
    }

    # Track homophily
    edge_ends <- ends(g, E(g))
    same <- sum(types[edge_ends[, 1]] == types[edge_ends[, 2]])
    homophily[step] <- same / ecount(g)

    # Track minority stress
    minority_nodes <- which(types == -1)
    stress_vals <- sapply(minority_nodes, function(v) {
      nbs <- neighbors(g, v)
      if (length(nbs) == 0) return(0)
      mean(types[nbs] != -1)
    })
    minority_stress[step] <- mean(stress_vals)
  }

  return(list(graph = g,
              homophily = homophily,
              stress = minority_stress,
              types = V(g)$type))
}

set.seed(42)
result <- henry_model_asymmetric_2D(n = 400, minority_frac = 0.3, max_steps = 1000)


plot(result$homophily, type = "l", lwd = 2, col = "darkred",
     ylim = c(0, 1), xlab = "Step", ylab = "Value",
     main = "Segregation Dynamics (BA)")
lines(result$stress, col = "steelblue", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("Homophily", "Minority Stress"),
       col = c("darkred", "steelblue"),
       lty = c(1, 2), lwd = 2, bty = "n")

