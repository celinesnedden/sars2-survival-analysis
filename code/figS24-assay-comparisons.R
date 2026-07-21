# This file: - Compares the trajectories across assay types


# Prep predictions -------------------------------------------------------------

# Load predictions
df <- fread( "./outputs/predictions/pred-across-cofactors-1000.csv")

# Subset to dose of 10^4 pfu from IN+IT individuals; female, adult, RM
df <- subset(df, dose_total == 4 & route_idx == 3 & sp_idx == 1 & age_idx == 2 & sex_idx == 0)

# Set all cofactor names
df <- assign_all_names(df)

# Convert into a long dataframe for plotting
df.long <- df %>%
  pivot_longer(cols = c(percent_positive, first_pos_median, peak_median, 
                        titer_mean, duration_median, last_median, auc),
               names_to = "metric", values_to = "value")


# Subset to relevant metrics
df.long <- subset(df.long, metric %in% c("percent_positive", "first_pos_median", "peak_median", "last_median"))
df.long <- subset(df.long, !(metric %in% c("peak_median", "last_median") & tissue_name == "Upper GI"))

# Rename the metrics
df.long$metric[df.long$metric == "percent_positive"] <- "Probability of positivity"
df.long$metric[df.long$metric == "first_pos_median"] <- "Time to detectability"
df.long$metric[df.long$metric == "peak_median"] <- "Time to peak titer"
df.long$metric[df.long$metric == "last_median"] <- "Time to undetectability"                                                       
    
df.long$metric <- factor(df.long$metric, levels = c("Probability of positivity",
                                                    "Time to detectability",
                                                    "Time to peak titer",
                                                    "Time to undetectability" ))
                                                     
# Plot -------------------------------------------------------------------------

## Event times -----------------------------------------------------------------

fig.times <- ggplot(subset(df.long, metric != "Probability of positivity")) +
  geom_density_ridges(data = subset(df.long,  metric == "Time to detectability"),  # Plotting % below so the legend will show it
                      aes(x = value, y = assay_name, fill = "Probability of positivity"),
                      alpha = 1, scale = 1, rel_min_height = 0.00001,
                      linewidth = 0.5) +
  geom_density_ridges(aes(x = value, y = assay_name, fill = metric),
                      alpha = 1, scale = 1, rel_min_height = 0.00001) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  coord_cartesian(xlim = c(0, 18)) +
  scale_fill_manual(values = c("Probability of positivity" = "#4A8396",
                               "Time to detectability" = "#4F609C",
                               "Time to peak titer" = "#9E4472",
                               "Time to undetectability" = "#628D56"),
                    limits = c("Probability of positivity",
                               "Time to detectability",
                               "Time to peak titer",
                               "Time to undetectability" ),
                    drop = "FALSE") +
  facet_wrap(tissue_name ~., nrow = 1) +
  labs(x = "Days post infection", y = "Detection Assay") +
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background  = element_blank(),
        strip.text = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "white",
                                         size = 0.2, linetype = "solid"),
        legend.margin = margin(c(5, 5, 5, 2)),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.3, color = "grey88"), 
        panel.grid.minor = element_line(linewidth = 0.2, color = "grey88")); fig.times


## Probabilty of positivity ----------------------------------------------------

fig.percent <- ggplot(subset(df.long, metric == "Probability of positivity")) +
  geom_density_ridges(aes(x = value * 100, y = assay_name, fill = metric),
                      alpha = 0.8, scale = 1, rel_min_height = 0.0001) +
  scale_x_continuous(breaks = seq(0, 100, 25)) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_fill_manual(values = c("Probability of positivity" = "#4A8396",
                               "Time to detectability" = "#4F609C",
                               "Time to peak titer" = "#9E4472",
                               "Time to undetectability" = "#628D56")) +
  facet_wrap(tissue_name ~ ., nrow = 1) +
  labs(x = "Probability of positivity (%)", y = "Detection Assay") +
  theme(text = element_text(size = 14),
        legend.position = "none",
        legend.title = element_blank(),
        strip.background  = element_rect(colour = "black", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "white",
                                         size = 0.2, linetype = "solid"),
        legend.margin = margin(c(5, 5, 5, 2)),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.3, color = "grey88"), 
        panel.grid.minor = element_line(linewidth = 0.2, color = "grey88")); fig.percent


## Combine ---------------------------------------------------------------------

fig.comb <- fig.percent + fig.times + plot_layout(nrow = 2)


## Save ---------------------------------------------------------------------

ggsave('./outputs/figures/figS24-assay-comparisons.png',
       plot = fig.comb,
       width = 12, 
       height = 5,
       dpi = 600)

