# This file: - Creates tables with parameter values (medians & 90% CI) for all key metrics

# Prep -------------------------------------------------------------------------

fit <- readRDS("./outputs/fits/fit-main.RDS")

# Extract the draws from all parameters
df.draws <- fit$draws(format = "df")

# Remove columns with "true" in their names, because they estimate true event
#    times for each individual
df.draws <- df.draws[, !grepl("true_first", names(df.draws))]
df.draws <- df.draws[, !grepl("true_peak", names(df.draws))]
df.draws <- df.draws[, !grepl("true_last", names(df.draws))]

# Remove the log likelihood
df.draws <- df.draws[, !grepl("lp_", names(df.draws))]

# Remove lab parameters (shown elsewhere)
df.draws <- df.draws[, !grepl("lab", names(df.draws))]

# Extract the medians & quantiles
df.summary <- summarise_draws(df.draws, 
                              median, 
                              ~quantile(.x, probs = c(0.05, 0.95)))

# Round the median & quantiles to 2 digits
df.summary$median_rounded <- round(df.summary$median, digits = 2)
df.summary$q5_rounded <- round(df.summary$`5%`, digits = 2)
df.summary$q95_rounded <- round(df.summary$`95%`, digits = 2)
df.summary$quantiles_combined <- paste0(df.summary$median_rounded,
                                        " [", df.summary$q5_rounded,
                                        ", ", df.summary$q95_rounded,
                                        "]")

# The below functions are only used in this file, so we keep them here
prep_param_estimates <- function(df, metric) {
  
  # Subset to only the parameters for the given metric
  df <- df[str_detect(df$variable, metric), ]
  
  # Remove the adjusted parameters (only used for fitting)
  df <- df[!str_detect(df$variable, "adj"), ]
  
  # Remove the metric header for all parameters
  df$variable <- str_remove_all(df$variable, metric)
  
  # Remove median name and underscores
  df$variable <- str_remove_all(df$variable, "median")
  df$variable <- str_remove_all(df$variable, "_")
  
  # Fix the shape parameter name
  df$variable[str_detect(df$variable, "shape")] <- str_remove_all(df$variable[str_detect(df$variable, "shape")], "intercept")
  
  # Fix the SD parameter names
  df$variable[str_detect(df$variable, "obsinterceptsd")] <- str_replace(df$variable[str_detect(df$variable, "obsinterceptsd")], "obsinterceptsd", "SD Observed Titer Intercept")
  df$variable[str_detect(df$variable, "truesd")] <- str_replace(df$variable[str_detect(df$variable, "truesd")], "truesd", "SD True Titer")
  
  # Extract the cofactor type
  df$cofactor <- sub("\\[.*", "", df$variable)
  
  # Make the cofactor names capitalized
  df$cofactor[!str_detect(df$cofactor, "SD")] <- str_to_sentence(df$cofactor[!str_detect(df$cofactor, "SD")])
  
  # Get the organ group category for each cofactor type
  cofactor_names <- c("Intercept", "Shape", "Location", "Sex", 
                      "SD Observed Titer Intercept", "SD True Titer")
  df$location_id[df$cofactor %in% cofactor_names] <- 
    sub(".*\\[([0-9])[^0-9]?.*", "\\1", df$variable[df$cofactor %in% cofactor_names])
  df$location_id[df$cofactor %notin% cofactor_names] <- 
    sub(".*,\\s*([0-9])[^0-9]?.*", "\\1", df$variable[df$cofactor %notin% cofactor_names])
  
  # Group by organ type for different panels
  #    Note route & dose are more specific, so they need to be handled differently
  cofactor_names2 <- c("Route", "Dose")
  df$organ_group <- df$location_id
  df$organ_group[df$cofactor %in% cofactor_names2 &  # URT tissues
                           df$location_id %in% c(1, 2)] <- 1
  df$organ_group[df$cofactor %in% cofactor_names2 &  # LRT tissues
                           df$location_id %in% c(3, 4)] <- 2
  df$organ_group[df$cofactor %in% cofactor_names2 &  # GI tissues
                           df$location_id %in% c(5, 6)] <- 3
  
  # Get the cofactor level for categorical cofactors
  cofactor_names3 <- c("Route", "Dose", "Age", "Species", "Assay", "Time")
  df$cofactor_level <- 1 # set as 1 for other cofactors
  df$cofactor_level[df$cofactor %in% cofactor_names3] <-
    sub(".*\\[([0-9])[^0-9]?.*", "\\1", df$variable[df$cofactor %in% cofactor_names3])
  
  # Set cofactor names
  df$cofactor_name <- " "
  df$cofactor_name[df$cofactor == "Age" & df$cofactor_level == 1] <- "Juvenile"
  df$cofactor_name[df$cofactor == "Age" & df$cofactor_level == 2] <- "Adult"
  df$cofactor_name[df$cofactor == "Age" & df$cofactor_level == 3] <- "Geriatric"
  df$cofactor_name[df$cofactor == "Species" & df$cofactor_level == 1] <- "RM"
  df$cofactor_name[df$cofactor == "Species" & df$cofactor_level == 2] <- "CM"
  df$cofactor_name[df$cofactor == "Species" & df$cofactor_level == 3] <- "AGM"
  df$cofactor_name[df$cofactor == "Assay" & df$cofactor_level == 1] <- "Total RNA"
  df$cofactor_name[df$cofactor == "Assay" & df$cofactor_level == 2] <- "gRNA"
  df$cofactor_name[df$cofactor == "Assay" & df$cofactor_level == 3] <- "sgRNA"
  df$cofactor_name[df$cofactor == "Assay" & df$cofactor_level == 4] <- "Culture"
  df$cofactor_name[df$cofactor == "Route" & df$cofactor_level == 1] <- "IN"
  df$cofactor_name[df$cofactor == "Route" & df$cofactor_level == 2] <- "IT"
  df$cofactor_name[df$cofactor == "Route" & df$cofactor_level == 3] <- "IN+IT"
  df$cofactor_name[df$cofactor == "Route" & df$cofactor_level == 4] <- "AE"
  df$cofactor_name[df$cofactor == "Route" & df$cofactor_level == 5] <- "IG"
  df$cofactor_name[df$cofactor == "Dose" & df$cofactor_level == 1] <- "Nose"
  df$cofactor_name[df$cofactor == "Dose" & df$cofactor_level == 2] <- "Throat"
  df$cofactor_name[df$cofactor == "Dose" & df$cofactor_level == 3] <- "Trachea"
  df$cofactor_name[df$cofactor == "Dose" & df$cofactor_level == 4] <- "Lung"
  df$cofactor_name[df$cofactor == "Dose" & df$cofactor_level == 5] <- "Upper GI"
  df$cofactor_name[df$cofactor == "Time" & df$cofactor_level == 1] <- "Not Inoculated"
  df$cofactor_name[df$cofactor == "Time" & df$cofactor_level == 2] <- "Inoculated"
  
  # Reset SD cofactor names
  df$cofactor_name[df$cofactor == "SD Observed Titer Intercept"] <- "Observed Titer Intercept"
  df$cofactor_name[df$cofactor == "SD True Titer"] <- "True Titer"
  df$cofactor[df$cofactor %in% c("SD Observed Titer Intercept", "SD True Titer")] <- "SD"
  
  # Factor cofactor names for the correct ordering in the figure
  df$cofactor_name <- factor(df$cofactor_name,
                             levels = c("AGM", "CM", "RM",
                                        "Geriatric", "Adult", "Juvenile",
                                        "Culture", "sgRNA", "gRNA",
                                        "Total RNA", 
                                        "Upper GI", "Lung", "Trachea",
                                        "Throat", "Nose",
                                        "IG", "AE", "IN+IT", "IT", "IN",
                                        "Not Inoculated", "Inoculated",
                                        "Observed Titer Intercept",
                                        "True Titer",
                                        " "))
  
  # Combine text for a couple cells
  df.collapse <- df %>%
    arrange(cofactor, cofactor_name, organ_group, location_id) %>%
    group_by(cofactor, organ_group, cofactor_name) %>%  # Include other grouping vars if needed
    summarise(
      quantiles_combined = paste(quantiles_combined, collapse = "  |  "),
      .groups = "drop"
    )
  
  # Set the organ group names for plotting
  df.collapse$organ_name[df.collapse$organ_group == 1] <- "URT"
  df.collapse$organ_name[df.collapse$organ_group == 2] <- "LRT"
  df.collapse$organ_name[df.collapse$organ_group == 3] <- "GI"
  df.collapse$organ_name <- factor(df.collapse$organ_name,
                               levels = c("URT", "LRT", "GI"))
  
  # Change cofactor label from time to titer for time to undetectability
  if (metric == "last") {
    df.collapse$cofactor[df.collapse$cofactor == "Time"] <- "Titer"
  }
  
  
  # Set the factor for the cofactor groups for plotting
  df.collapse$cofactor <- factor(df.collapse$cofactor,
                             levels = c("Shape", "Intercept", "Location",
                                        "Route", "Dose", "Age", "Sex",
                                        "Species", "Assay", "Time", "Titer",
                                        "SD"))
  

  return(df.collapse)
  
}

plot_param_estimates <- function(df) {
  
  table <- ggplot(df) +
    geom_tile(aes(x = 1, y = cofactor_name), 
              fill = "white", color = "black", 
              linewidth = 0.5) +
    geom_text(aes(x = 1, y = cofactor_name, label = quantiles_combined)) +
    facet_grid(cofactor ~ organ_name, scales = "free", space = "free") +
    coord_cartesian(expand = FALSE, clip = "off") +
    labs(x = "", y = "Cofactor") +
    theme(axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          strip.background = element_rect(fill = "grey92", color = "black"),
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold", angle = 0, hjust = 0))
  
  return(table)
}


# Table S4: Probability of Positivity ------------------------------------------

df.percent <- prep_param_estimates(df.summary, "percent")
tblS4 <- plot_param_estimates(df.percent); tblS4

ggsave('./outputs/tables/tblS4-parameters-percent-positivity.png',
       plot = tblS4,
       width = 11, 
       height = 7,
       dpi = 600)


# Table S5: Time to Detectability ------------------------------------------

df.first <- prep_param_estimates(df.summary, "first")
tblS5 <- plot_param_estimates(df.first); tblS5

ggsave('./outputs/tables/tblS5-parameters-detectability.png',
       plot = tblS5,
       width = 11, 
       height = 7,
       dpi = 600)


# Table S6: Time to Peak -------------------------------------------------------

df.peak <- prep_param_estimates(df.summary, "peak")

# Pull out parameters that are just set to zero (for indexing reasons in the model)
df.peak$quantiles_combined <- str_replace(df.peak$quantiles_combined, 
                                          "^0\\s*\\[0\\s*,\\s*0\\s*\\]\\s*\\|\\s*", 
                                          "NA | ")
df.peak$quantiles_combined <- str_replace(df.peak$quantiles_combined, 
                                          "^\\s*0\\s*\\[0\\s*,\\s*0\\s*\\]$", 
                                          "NA")

# Plot
tblS6 <- plot_param_estimates(df.peak); tblS6

ggsave('./outputs/tables/tblS6-parameters-peak-time.png',
       plot = tblS6,
       width = 11, 
       height = 7,
       dpi = 600)


# Table S7: Peak Titer ---------------------------------------------------------

df.titer <- prep_param_estimates(df.summary, "titer")

# Pull out parameters that are just set to zero (for indexing reasons in the model)
df.titer$quantiles_combined <- str_replace(df.titer$quantiles_combined, 
                                          "^0\\s*\\[0\\s*,\\s*0\\s*\\]\\s*\\|\\s*", 
                                          "NA | ")
df.titer$quantiles_combined <- str_replace(df.titer$quantiles_combined, 
                                          "^\\s*0\\s*\\[0\\s*,\\s*0\\s*\\]$", 
                                          "NA")

# Plot
tblS7 <- plot_param_estimates(df.titer); tblS7

ggsave('./outputs/tables/tblS7-parameters-peak-titer.png',
       plot = tblS7,
       width = 11, 
       height = 7,
       dpi = 600)


# Table S8: Time to Undetectability --------------------------------------------

df.last <- prep_param_estimates(df.summary, "last")

# Pull out parameters that are just set to zero (for indexing reasons in the model)
df.last$quantiles_combined <- str_replace(df.last$quantiles_combined, 
                                           "^0\\s*\\[0\\s*,\\s*0\\s*\\]\\s*\\|\\s*", 
                                           "NA | ")
df.last$quantiles_combined <- str_replace(df.last$quantiles_combined, 
                                           "^\\s*0\\s*\\[0\\s*,\\s*0\\s*\\]$", 
                                           "NA")

# Plot
tblS8 <- plot_param_estimates(df.last); tblS8

ggsave('./outputs/tables/tblS8-parameters-undetectability.png',
       plot = tblS8,
       width = 11, 
       height = 7,
       dpi = 600)

