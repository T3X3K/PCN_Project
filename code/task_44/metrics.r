library(igraph)
library(tidyverse)

data_folder <- "data/"
edge_files <- list.files(data_folder, pattern = "^weighted_edges.*\\.csv$", full.names = TRUE)

results_list <- list()

for (f in edge_files) {
  df <- read_csv(f, show_col_types = FALSE)
  colnames(df) <- c("nodeID_from", "nodeID_to", "user_country", "country_ISO3", "SCI")

  country_name <- unique(df$user_country)[1]
  cat("Processing:", country_name, "\n")

  g <- graph_from_data_frame(df %>% select(nodeID_from, nodeID_to, SCI), directed = FALSE)

  comps <- components(g)
  g_lcc <- induced_subgraph(g, which(comps$membership == which.max(comps$csize)))

  if (vcount(g_lcc) < 5) next

  deg <- degree(g_lcc)
#   btw <- betweenness(g_lcc, normalized = TRUE)
btw <- betweenness(g_lcc, weights = 1 / E(g_lcc)$SCI, normalized = TRUE)


  results_list[[country_name]] <- tibble(
    country = country_name,
    mean_degree = mean(deg, na.rm = TRUE),
    mean_betweenness = mean(btw, na.rm = TRUE)
  )
}

results_df <- bind_rows(results_list)

ggplot(results_df, aes(x = mean_degree, y = mean_betweenness, label = country)) +
  geom_point(color = "tomato", size = 2.5) +
  geom_text(vjust = -0.6, size = 2.5) +
  labs(
    title = "Mean Degree vs Mean Betweenness Centrality (per country)",
    x = "Mean Degree",
    y = "Mean Betweenness Centrality (normalized)"
  ) +
  theme_minimal()


ggsave("mean_degree_vs_betweenness.png", width = 8, height = 6)





library(igraph)
library(tidyverse)

data_folder <- "data/"
edge_files <- list.files(data_folder, pattern = "^weighted_edges.*\\.csv$", full.names = TRUE)

centrality_list <- list()

analyze_country_network <- function(df, country_name) {
  g <- graph_from_data_frame(df %>% select(nodeID_from, nodeID_to, SCI), directed = FALSE)

  clust <- components(g)
  g_lcc <- induced_subgraph(g, which(clust$membership == which.max(clust$csize)))

  if (vcount(g_lcc) < 5) return(NULL)

  eig <- eigen_centrality(g_lcc, weights = E(g_lcc)$SCI)$vector
  sci_totals <- graph.strength(g_lcc, weights = E(g_lcc)$SCI)

  centrality_df <- data.frame(
    node = names(eig),
    eigenvector = eig,
    total_SCI = sci_totals,
    country = country_name
  )

  return(centrality_df)
}

for (f in edge_files) {
  df <- read_csv(f, show_col_types = FALSE)
  colnames(df) <- c("nodeID_from", "nodeID_to", "user_country", "country_ISO3", "SCI")
  country_name <- unique(df$user_country)[1]

  result <- analyze_country_network(df, country_name)

  if (!is.null(result)) {
    centrality_list[[country_name]] <- result
  }

  rm(df, result)
  gc()
}

centrality_df <- bind_rows(centrality_list)

centrality_means <- centrality_df %>%
  group_by(country) %>%
  summarise(
    mean_eigenvector = mean(eigenvector, na.rm = TRUE),
    mean_total_SCI = mean(total_SCI, na.rm = TRUE)
  )

ggplot(centrality_means, aes(x = mean_total_SCI, y = mean_eigenvector, label = country)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  labs(
    title = "Mean Eigenvector Centrality vs Mean Total SCI",
    x = "Mean Total SCI (per node)",
    y = "Mean Eigenvector Centrality"
  ) +
  theme_minimal()

ggsave('eig.png')
