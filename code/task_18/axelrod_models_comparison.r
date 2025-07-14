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


run <- function(g, n_traits, n_features, check_interval = 1000, max_steps = 200000) {
  V(g)$label <- lapply(1:vcount(g), function(x) sample(1:n_traits, n_features, replace=TRUE))
  active_pairs <- 1
  i <- 0

  active_trace <- c()
  region_trace <- c()

  while (active_pairs != 0 && i < max_steps) {
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
      regions <- count_cultural_regions(g)
      active_trace <- c(active_trace, active_pairs)
      region_trace <- c(region_trace, regions)
    }
  }

  return(list(active_trace = active_trace, region_trace = region_trace))
}



lattice_dimensions <- 2
lattice_side_length <- 10
n_traits <- 15
n_features <- 5


n_reps <- 10
check_interval <- 1000

all_data <- data.frame()

for (i in 1:n_reps) {
  for (model in c("lattice", "ba", "ws", "sbm")) {
    # Create network
    g <- switch(model,
      lattice = make_lattice(c(10, 10)),
      ba      = sample_pa(100, m = 2, directed = FALSE),
      ws      = sample_smallworld(1, 100, nei = 2, p = 0.1),
      sbm     = sample_sbm(100, pref.matrix = matrix(c(0.9, 0.1, 0.1, 0.9), 2), block.sizes = c(50, 50))
    )

    result <- run(g, n_traits, n_features, check_interval)

    trace_len <- length(result$region_trace)
    df <- data.frame(
      time = 1:trace_len * check_interval,
      region = result$region_trace,
      active = result$active_trace,
      model = model,
      run_id = i
    )

    all_data <- bind_rows(all_data, df)
  }
}


summary_df <- all_data %>%
  group_by(model, time) %>%
  summarise(
    mean_regions = mean(region, na.rm = TRUE),
    sd_regions = sd(region, na.rm = TRUE),
    mean_active = mean(active, na.rm = TRUE),
    sd_active = sd(active, na.rm = TRUE)
  )

ggplot(summary_df, aes(x = time, y = mean_regions, color = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_regions - sd_regions, ymax = mean_regions + sd_regions, fill = model), alpha = 0.2, color = NA) +
  labs(title = "Cultural Regions Over Time", y = "Number of Cultural Regions", x = "Time") +
  theme_minimal()

ggsave("cultural_regions_over_time.png", width = 10, height = 6)
write.csv(summary_df, "summary_df.csv", row.names = FALSE)
