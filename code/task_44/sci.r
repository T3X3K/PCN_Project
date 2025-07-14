library(dplyr) |> suppressMessages()
library(tidyr) |> suppressMessages()
library(terra) |> suppressMessages()
library(geodata) |> suppressMessages()
library(eurostat) |> suppressMessages()
library(sf) |> suppressMessages()
library(tigris) |> suppressMessages()
library(igraph) |> suppressMessages()
library(ggplot2) |> suppressMessages()
library(purrr) |> suppressMessages()
library(stringr) |> suppressMessages()
library(magrittr) |> suppressMessages()
library(ggraph) |> suppressMessages()
library(readr) |> suppressMessages()
library(countrycode) |> suppressMessages()


path_data_fb <- 'gadm1_nuts3_counties-gadm1_nuts3_counties - FB Social Connectedness Index - October 2021.tsv'
df <- read.table(file = path_data_fb, sep = "\t", header = TRUE) 
levels <- read.table(file = 'gadm1_nuts3_counties_levels.csv', sep = ",", header = TRUE)

levels <- levels[levels$level != 'county', ]
levels <- levels[levels$level != 'country', ]



add_iso_id <- function(df, level_col = "level", key_col = "key") {
    
    if (is.character(df)) {
        stop("Input must be a data.frame or tibble, not a character vector.")
    }
    df <- df %>%
        mutate(
            iso = case_when(
                !!sym(level_col) == "country" ~ !!sym(key_col),
                !!sym(level_col) %in% c("gadm1", "gadm2") ~ substr(!!sym(key_col), 1, 3),
                TRUE ~ NA_character_
            ),
            id = case_when(
                !!sym(level_col) == "country" ~ 0L,
                !!sym(level_col) %in% c("gadm1", "gadm2") ~ as.integer(gsub("\\D", "", !!sym(key_col))),
                TRUE ~ NA_integer_
            )
        ) %>%
        mutate(
            name = countrycode::countrycode(iso, origin = "iso3c", destination = "country.name.en")
        )
    return(df)
}

levels <- add_iso_id(levels)
levels <- levels %>%
  group_by(name) %>%
  mutate(id = row_number()) %>%  
  ungroup()


nuts3_shapefile <- eurostat::get_eurostat_geospatial(nuts_level = 3, year = 2024, output_class = "sf")
nuts3_idx <- which(levels$level == "nuts3")
nuts3_keys <- levels$key[nuts3_idx]

matched <- match(nuts3_keys, nuts3_shapefile$id)
cntr_codes <- nuts3_shapefile$CNTR_CODE[matched]
cntr_codes[cntr_codes == "EL"] <- "GR"
cntr_codes[cntr_codes == "UK"] <- "GB"

iso3c_codes <- countrycode::countrycode(cntr_codes, origin = "iso2c", destination = "iso3c")
levels$iso[nuts3_idx] <- iso3c_codes
levels$name[nuts3_idx] <- countrycode::countrycode(iso3c_codes, origin = "iso3c", destination = "country.name.en")


levels$id[nuts3_idx] <- sub("^[A-Z]+", "", nuts3_keys)

df$user_country <- levels$name[match(df$user_loc, levels$key)]
df$fr_country <- levels$name[match(df$fr_loc, levels$key)]
df$user_id <- levels$id[match(df$user_loc, levels$key)]
df$fr_id <- levels$id[match(df$fr_loc, levels$key)]
df$user_iso <- levels$iso[match(df$user_loc, levels$key)]

cleaned_df <- df %>%
  filter(!is.na(user_country)) %>%
  filter(!is.na(fr_country)) %>%
  filter(user_country == fr_country) %>%
  filter(user_country %in% levels$name)


country_sci <- cleaned_df %>%
  group_by(user_country) %>%
  summarise(total_sci = sum(scaled_sci, na.rm = TRUE)) %>%
  arrange(desc(total_sci)) %>%
  slice_head(n = 100)

top_countries <- country_sci$user_country




gadm1 <- st_read("../gadm28/gadm28_adm1.shp")

gadm1_centroids <- gadm1 %>%
  st_make_valid() %>%
  mutate(
    key = paste0(ISO, "_", ID_1),
    centroid = st_point_on_surface(geometry),
    longitude = st_coordinates(centroid)[, 1],
    latitude = st_coordinates(centroid)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(key, latitude, longitude)

levels <- levels %>%
  left_join(bind_rows(gadm1_centroids), by = "key")
rm(gadm1, gadm1_centroids)
gc()

gadm2 <- st_read("../gadm28/gadm28_adm2.shp")

gadm2_centroids <- gadm2 %>%
  st_make_valid() %>%
  mutate(
    key = paste0(ISO, ID_1, "_", ID_2),
    centroid = st_point_on_surface(geometry),
    longitude = st_coordinates(centroid)[, 1],
    latitude = st_coordinates(centroid)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(key, latitude, longitude)

levels <- levels %>%
  left_join(bind_rows(gadm2_centroids), by = "key")
rm(gadm2, gadm2_centroids)
gc()


nuts3_sf <- get_eurostat_geospatial(
  output_class = "sf",
  resolution = "60",
  nuts_level = 3,
  year = 2021
)

nuts3_centroids <- nuts3_sf %>%
  st_make_valid() %>%
  mutate(
    key = NUTS_ID,
    centroid = st_point_on_surface(geometry),
    longitude = st_coordinates(centroid)[, 1],
    latitude = st_coordinates(centroid)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(key, latitude, longitude)


levels <- levels %>%
  left_join(bind_rows(nuts3_centroids), by = "key")

rm(nuts3_sf, nuts3_centroids)
gc()



final_df <- cleaned_df %>%
  filter(user_country %in% top_countries)

write_csv(final_df, "top100_cleaned_sci_data.csv")
write_csv(country_sci, "top100_country_sci_totals.csv")

selected_df <- final_df %>%
  select(user_id, fr_id, user_country, user_iso) %>%
  setNames(c('nodeID_from','nodeID_to','country_name','country_ISO3'))

selected_df_weighted <- final_df %>%
  select(user_id, fr_id, user_country, user_iso, scaled_sci) %>%
  setNames(c('nodeID_from','nodeID_to','country_name','country_ISO3', 'SCI'))

levels_with_coords <- levels %>%
  filter(name %in% top_countries) %>%
  select(id, key, latitude, longitude, name) %>%
  setNames(c('nodeID','nodeLabel','latitude','longitude', 'country'))


write_csv(selected_df, "top100_selected_sci_data.csv")
write.csv(levels_with_coords, "user_locations_with_coordinates.csv")


print(head(selected_df))
print(unique(selected_df$country_name))
print(nrow(selected_df))


selected_df <- final_df %>%
  filter(user_country %in% top_countries) %>%
  select(user_id, fr_id, user_country, user_iso) %>%
  setNames(c('nodeID_from','nodeID_to','country_name','country_ISO3'))

selected_df_weighted <- final_df %>%
  filter(user_country %in% top_countries) %>%
  select(user_id, fr_id, user_country, user_iso, scaled_sci) %>%
  setNames(c('nodeID_from','nodeID_to','country_name','country_ISO3', 'SCI'))

levels_with_coords <- levels %>%
  filter(name %in% top_countries) %>%
  select(id, key, latitude, longitude, name) %>%
  setNames(c('nodeID','nodeLabel','latitude','longitude', 'country'))



for (country in unique(selected_df$country_name)) {
  df_country <- selected_df %>% filter(country_name == country)
  print(df_country)
  safe_country <- gsub("[^A-Za-z0-9]", "_", country)
  filename <- paste0("data/edges_", safe_country, ".csv")
  write_csv(df_country, filename)
}


for (country in unique(levels_with_coords$country)) {
  df_country <- levels_with_coords %>% 
                filter(country == country) %>%
                select(nodeID, nodeLabel, latitude, longitude)  
  safe_country <- gsub("[^A-Za-z0-9]", "_", country)  
  filename <- paste0("data/nodes_", safe_country, ".csv")
  write_csv(df_country, filename)
}


for (country in unique(selected_df_weighted$country_name)) {
  df_country <- selected_df_weighted %>% filter(country_name == country)
  safe_country <- gsub("[^A-Za-z0-9]", "_", country)
  filename <- paste0("data/weighted_edges", safe_country, ".csv")
  write_csv(df_country, filename)  
}
