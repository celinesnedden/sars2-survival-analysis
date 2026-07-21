# This file: - Plots trajectories across full dose set, for culture and all routes


# Prep -------------------------------------------------------------------------

# Load predictions
df <- fread( "./outputs/predictions/pred-across-doses-450.csv")

# Subset to culture data, from relevant tissues, & fixed cofactor set
df.dose.fixed <- subset(df, tissue_idx %in% c(1, 4, 6) & assay_idx == 4 &
                          sp_idx == 1 & age_idx == 2 & sex_idx == 0)

# Set groups by sample number, tissue, and assay, for plotting
df.dose.fixed$group <- paste0(df.dose.fixed$sample_num, "-", 
                              df.dose.fixed$tissue_idx, "-",
                               df.dose.fixed$assay_idx)

# Calculate medians & quantiles for plotting
df.dose.quantiles <- df.dose.fixed %>%
  group_by(dose_total, tissue_idx, route_idx, assay_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 


# Assign route names
df.dose.quantiles <- assign_route_names(df.dose.quantiles)


# Plot -------------------------------------------------------------------------

## Nose -------------------------------------------------------------------------

fig.traj.nose <- 
  ggplot(subset(df.dose.quantiles, tissue_idx == 1)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_grid(route_name ~ "Nose") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#D8CDDF", high = df.color$Nose) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(2, 2, 0, 2),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.traj.nose


## Lung -------------------------------------------------------------------------

fig.traj.lung <- 
  ggplot(subset(df.dose.quantiles, tissue_idx == 4)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_grid(route_name ~ "Lung") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FBEFC6", high = df.color$Lung) +
  theme(text = element_text(size = 11),
        legend.position =  "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title = element_text(size = 11),
        axis.title.y = element_blank(),
        plot.margin = margin(0, 2, 0, 2),
        strip.text.y = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.traj.lung


## Lower GI -------------------------------------------------------------------------

fig.traj.gi <- 
  ggplot(subset(df.dose.quantiles, tissue_idx == 6)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_grid(route_name ~ "Lower GI") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FFDAD6", high = df.color$Lower.GI) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.traj.gi


## Legend -----------------------------------------------------------------------

fig.traj.legend <- 
  ggplot(df.dose.quantiles) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14)) +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "Days post infection", y = "Viral titer (log10)", 
       color = "Dose\n(log10 pfu)") +
  scale_color_gradient(low = "grey88", high = "grey22") +
  theme(text = element_text(size = 9),
        legend.position = c(0.3, 0.8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_blank()); fig.traj.legend

# Extract the legend
g <- ggplotGrob(fig.traj.legend)
legend <- g$grobs[[which(g$layout$name == "guide-box-inside")]] # legend at the bottom
fig.legend <- ggdraw(legend); fig.legend # wrap as ggplot compatible item


# Combine --------------------------------------------------------------- 

fig.traj <- fig.traj.nose + 
  fig.traj.lung + fig.traj.gi + 
  fig.legend + 
  plot_layout(ncol = 4, widths = c(1, 1, 1, 1)); fig.traj


# Save --------------------------------------------------------------- 

ggsave("./outputs/figures/figS28-dose-effects-trajectories-all-routes-culture.png",
       plot = fig.traj,
       width = 6, 
       height = 5.2,
       dpi = 600)
