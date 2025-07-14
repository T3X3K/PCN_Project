library(ggplot2) |> suppressMessages()
library(igraph) |> suppressMessages()
library(dplyr) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()


similarity_probability <- function(g, i, j, n_features) {
    sum(V(g)$label[[i]] == V(g)$label[[j]])/n_features}


return_active_pairs <- function(g, n_features) {
  active_pairs <- 0
  for (v in V(g)) {
    for (n in neighbors(g, v)) {
      sim <- similarity_probability(g, v, n, n_features)
      if (sim > 0 && sim < 1) {
        active_pairs <- active_pairs + 1
      }
    }
  }
  return(active_pairs)
}

count_cultural_regions <- function(g) {
  labels <- sapply(V(g)$label, function(x) paste(x, collapse = "-"))
  regions <- 0
  seen <- rep(FALSE, vcount(g))

  for (i in seq_along(labels)) {
    if (!seen[i]) {
      same_label_nodes <- which(labels == labels[i])
      subgraph <- induced_subgraph(g, same_label_nodes)
      comps <- components(subgraph)
      regions <- regions + comps$no
      seen[same_label_nodes] <- TRUE
    }
  }
  return(regions)
}


run <- function(g, n_traits, n_features) {

    V(g)$label <- lapply(1:vcount(g), function(x) sample(1:n_traits, n_features, replace=TRUE))
    active_pairs <- 1
    i <- 0
    check_interval <- 1000  # Check active pairs every N steps

    while (active_pairs != 0) {
      e2 <- sample(V(g), 1)
      e1 <- sample(neighbors(g, e2), 1)
      similarity <- similarity_probability(g, e1, e2, n_features)

      if (runif(1) < similarity && similarity != 1) {
        diff_features <- which(V(g)$label[[e1]] != V(g)$label[[e2]])
        if (length(diff_features) > 0) {
          feature_to_change <- sample(diff_features, 1)
          V(g)$label[[e2]][feature_to_change] <- V(g)$label[[e1]][feature_to_change]
        }
      }

      i <- i + 1
      if (i %% check_interval == 0) {
        active_pairs <- return_active_pairs(g, n_features)
        # print(paste("Step:", i, "Active pairs:", active_pairs))
      }
    }


    # Final output
    num_regions <- count_cultural_regions(g)

    return(num_regions)
}


lattice_dimensions <- 2
lattice_side_length <- 10
g <- make_lattice(length=lattice_side_length, dim=lattice_dimensions)

N_traits <- c(5, 10, 15)
N_features <- c(5, 10, 15)
# N_traits <- c(5, 10)#, 10, 15)
# N_features <- c(5, 10)#, 10, 15)

# results <- expand.grid(n_traits = N_traits, n_features = N_features) %>%
#     mutate(avg_regions = map2_dbl(n_traits, n_features, function(traits, features) {
#         mean(replicate(1, run(make_lattice(length = lattice_side_length, dim = lattice_dimensions), traits, features)))
#     }))

results <- expand.grid(n_traits = N_traits, n_features = N_features) %>%
    mutate(
        sims = map2(n_traits, n_features, function(traits, features) {
            replicate(10, run(make_lattice(length = lattice_side_length, dim = lattice_dimensions), traits, features))
        }),
        avg_regions = map_dbl(sims, mean),
        sd_regions = map_dbl(sims, sd)
    ) %>%
    select(-sims)

print(results)