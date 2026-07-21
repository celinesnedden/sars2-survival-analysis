# This file: - Generates tables with the predictions for each metric & tissue location, for all cofactors combined
# Dependencies: 

# Prep -------------------------------------------------------------------------

# Load predictions
df <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Subset to doses of interest
df <- subset(df, dose_total %in% c(4, 7) & assay_idx %in% c(1, 4))

# Set all cofactor names
df <- assign_all_names(df)

# Convert to long format for plotting
df.long <- df %>%
  pivot_longer(
    cols = all_of(c("percent_positive", "first_pos_median", "peak_median",     
                    "titer_mean", "last_median", "auc", "duration_median")),
    names_to = "metric_name",
    values_to = "value"
  )


# Set elaborated metric names
df.long$metric_name[df.long$metric_name == "percent_positive"] <- "Probability of positivity"
df.long$metric_name[df.long$metric_name == "first_pos_median"] <- "Time to detectability"
df.long$metric_name[df.long$metric_name == "peak_median"] <- "Time to peak titer"
df.long$metric_name[df.long$metric_name == "titer_mean"] <- "Peak titer"
df.long$metric_name[df.long$metric_name == "last_median"] <- "Time to undetectability"
df.long$metric_name[df.long$metric_name == "duration_median"] <- "Duration"
df.long$metric_name[df.long$metric_name == "auc"] <- "AUC"
df.long$metric_name <- factor(df.long$metric_name, 
                              levels = c("Probability of positivity",
                                         "Time to detectability",
                                         "Time to peak titer",
                                         "Peak titer",
                                         "Time to undetectability",
                                         "Duration", "AUC"))

# Calculate median & CI across routes, doses, locations, assays, and metrics
df.summarize <- df.long %>%
  group_by(location_name, route_name, dose_total, assay_name, metric_name) %>%
  summarise(
    n = sum(!is.na(value)),
    q5 = round(quantile(value, 0.05, na.rm = TRUE), digits = 2),
    median = round(median(value, na.rm = TRUE), digits = 2),
    q95 = round(quantile(value, 0.95, na.rm = TRUE), digits =2),
    .groups = "drop"
  )

# Express positivity as a percent
df.summarize$median[df.summarize$metric_name == "Probability of positivity"] <- 
  100 * df.summarize$median[df.summarize$metric_name == "Probability of positivity"]
df.summarize$q5[df.summarize$metric_name == "Probability of positivity"] <- 
  100 * df.summarize$q5[df.summarize$metric_name == "Probability of positivity"]
df.summarize$q95[df.summarize$metric_name == "Probability of positivity"] <- 
  100 * df.summarize$q95[df.summarize$metric_name == "Probability of positivity"]


# Remove Upper GI metrics, except prob positivity and time to detectability
df.summarize <- subset(df.summarize, !(metric_name %notin% c("Probability of positivity", "Time to detectability") &
                                         location_name == "Upper GI"))


plot_prediction_tables <- function(df, metrics) {
  
  df <- subset(df, metric_name %in% metrics)
  
  
  tbl <- ggplot(df) + 
    geom_tile(aes(x = route_name, y = as.character(dose_total)), 
              fill = "white", color = "black", linewidth = 0.25) + 
    geom_text(aes(x = route_name, y = as.character(dose_total), 
                  label = paste0(median, " [", 
                                 q5, ", ",
                                 q95, "]")),
              size = 2.5) + 
    scale_y_discrete(limits = c("7", "4"), labels = c(expression(10^7), expression(10^4))) +
    facet_grid(assay_name + location_name ~ metric_name, 
               scales = "free", space = "free") +
    coord_cartesian(expand = FALSE, clip = "off") +
    labs(x = "Exposure Route", y = "Exposure Dose (pfu)") + 
    theme(text = element_text(size = 9),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey92", color = "black"),
          strip.text.x = element_text(face = "bold"),
          strip.text.y = element_text(face = "bold", angle = 0, hjust = 0)); tbl
  
  return(tbl)
  
}


# Table S9: Prob. of Positivity -----------------------------------------

tblS9 <- plot_prediction_tables(df.summarize, c("Probability of positivity")); tblS9

ggsave('./outputs/tables/tblS9-predictions-percent-positivity.png',
       plot = tblS9,
       width = 6, 
       height = 6,
       dpi = 600)


# Table S10: Time to detectability ---------------------------------------------

tblS10 <- plot_prediction_tables(df.summarize, c("Time to detectability")); tblS10

ggsave('./outputs/tables/tblS10-predictions-detectability.png',
       plot = tblS10,
       width = 6, 
       height = 6,
       dpi = 600)


# Table S11: Peak Time ---------------------------------------------------------

tblS11 <- plot_prediction_tables(df.summarize, c("Time to peak titer")); tblS11

ggsave('./outputs/tables/tblS11-predictions-peak-time.png',
       plot = tblS11,
       width = 6, 
       height = 5.2,
       dpi = 600)


# Table S12: Peak Titer ---------------------------------------------------------

tblS12 <- plot_prediction_tables(df.summarize, c("Peak titer")); tblS12

ggsave('./outputs/tables/tblS12-predictions-peak-titer.png',
       plot = tblS12,
       width = 6, 
       height = 5.2,
       dpi = 600)


# Table S13: Time to undetectability -------------------------------------------

tblS13 <- plot_prediction_tables(df.summarize, c("Time to undetectability")); tblS13

ggsave('./outputs/tables/tblS13-predictions-undetectability.png',
       plot = tblS13,
       width = 6.5, 
       height = 5.2,
       dpi = 600)


# Table S14: Duration ----------------------------------------------------------

tblS14 <- plot_prediction_tables(df.summarize, c("Duration")); tblS14

ggsave('./outputs/tables/tblS14-predictions-duration.png',
       plot = tblS14,
       width = 6.5, 
       height = 5.2,
       dpi = 600)


# Table S15: AUC ---------------------------------------------------------------

tblS15 <- plot_prediction_tables(df.summarize, c("AUC")); tblS15

ggsave('./outputs/tables/tblS15-predictions-AUC.png',
       plot = tblS15,
       width = 6.5, 
       height = 5.2,
       dpi = 600)



