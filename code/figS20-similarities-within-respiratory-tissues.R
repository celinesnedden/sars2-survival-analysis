# This file: - Compares trajectories between the nose & throat, and the trachea & lung

# Prep -------------------------------------------------------------------------

# Load predictions
df <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Subset to dose of 10^4 pfu
df <- subset(df, dose_total == 4 & age_idx == 2 & sex_idx == 0 & sp_idx == 1)

# Set location and route names, with factors
df <- assign_all_names(df)

# Dataframe, with trajectories grouped by cofactor & sample number
df.traj <- df
df.traj$titer_mean <- abs(df.traj$titer_mean)
df.traj$group <- paste0(df.traj$dose_total, df.traj$route_idx,
                        df.traj$sp_idx, df.traj$age_idx,
                        df.traj$sex_idx, df.traj$assay_idx,
                        df.traj$sample_num, df.traj$tissue_idx)

# Get quantiles
df.traj.quantiles <- df.traj %>%
  group_by(tissue_idx, route_name, assay_name, assay_idx, organ_group) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 


# Plot ---------------------------------------------------------------------

# Choose a random sample from all of the trajectories
random_samps <- sample(unique(df.traj$sample_num), 100)

# Plot
fig.compare <- 
  ggplot(subset(df.traj, tissue_idx %notin% c(5, 6) & 
                  assay_idx %in% c(1, 4) &
                  sample_num %in% random_samps)) +
  geom_segment(aes(x = first_pos_median,
                   xend = peak_median,
                   y = 0, yend = titer_mean,
                   color = as.character(tissue_idx)),
               linewidth = 0.05, alpha = 0.2) +
  geom_segment(aes(x = peak_median,
                   xend = last_median,
                   yend = 0, y = titer_mean,
                   color = as.character(tissue_idx)),
               linewidth = 0.05, alpha = 0.2) +
  geom_segment(data = subset(df.traj.quantiles, 
                             tissue_idx %notin% c(5, 6) &
                               assay_idx %in% c(1, 4)),
               aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = as.character(tissue_idx)),
               linewidth = 0.5, alpha = 1) +
  geom_segment(data = subset(df.traj.quantiles, 
                             tissue_idx %notin% c(5, 6) &
                               assay_idx %in% c(1, 4)),
               aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = as.character(tissue_idx)),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 4)) +
  coord_cartesian(xlim = c(0, 16)) +
  facet_grid(factor(organ_group, levels = c("URT", "LRT", "GI")) + assay_name ~ 
               factor(route_name, levels = c("IN", "IT", "IN + IT", "AE", "IG")),
             scales = "free") +
  labs(x = "Days post infection", y = "Viral titer\n(log10 total RNA copies or pfu)",
       color = "Tissue") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "2" = df.color$Throat,
                                "3" = df.color$Trachea,
                                "4" = df.color$Lung),
                     labels = c("Nose", "Throat", "Trachea", "Lung")) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.1))) + 
  theme(text = element_text(size = 11),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.8, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text()); fig.compare


# Save ---------------------------------------------------------------------

ggsave('./outputs/figures/figS20-trajectory-similarities-within-respiratory-tissues.png',
       plot = fig.compare,
       width = 7, 
       height = 4.5,
       dpi = 600)

