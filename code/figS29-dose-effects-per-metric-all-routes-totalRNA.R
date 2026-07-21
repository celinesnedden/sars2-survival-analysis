# This file: - Plots as lines the relationships between dose & main metrics, for total RNA and all routes


# Prep -------------------------------------------------------------------------

# Load predictions
df <- fread( "./outputs/predictions/pred-across-doses-450.csv")

# Subset to culture data, from relevant tissues, & fixed cofactor set
df.dose.fixed <- subset(df, tissue_idx %in% c(1, 4, 6) & assay_idx == 1 &
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


# Panel A: Probability of positivity -------------------------------------------

fig.dose.prob.ci <- ggplot(df.dose.quantiles) + 
  geom_ribbon(aes(x = dose_total, 
                  ymin = percent_positive_q5 * 100, 
                  ymax = percent_positive_q95 * 100,
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, 
                y = percent_positive_median * 100,
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Days post infection") +
  facet_grid(route_name ~ "Probability of positivity") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.prob.ci



# Panel B: Time to detectability ----------------------------------------------------

fig.dose.first.ci <- ggplot(df.dose.quantiles) + 
  geom_ribbon(aes(x = dose_total, 
                  ymin = log2(first_pos_median_q5), 
                  ymax = log2(first_pos_median_q95),
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, 
                y = log2(first_pos_median_median),
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 3, 1),
                     labels = 2^seq(-2, 3, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Days post infection") +
  facet_grid(route_name ~ "Time to detectability") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.first.ci



# Panel C: Time to peak --------------------------------------------------------

fig.dose.peak <- ggplot(df.dose.quantiles) + 
  geom_ribbon(aes(x = dose_total, 
                  ymin = log2(peak_median_q5), 
                  ymax = log2(peak_median_q95),
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, 
                y = log2(peak_median_median),
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 3, 1),
                     labels = 2^seq(-2, 3, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Days post infection") +
  facet_grid(route_name ~ "Time to peak titer") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.peak


# Panel D: Peak titer ------------------------------------------------------------

fig.dose.titer <-ggplot(df.dose.quantiles) + 
  geom_ribbon(aes(x = dose_total, 
                  ymin = titer_mean_q5, 
                  ymax = titer_mean_q95,
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, y = titer_mean_median,
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Viral titer (log10 total RNA copies)") +
  facet_grid(route_name ~ "Peak titer") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.title.y = element_text(size = 9),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.titer



# Panel E: Time to undetectability ---------------------------------------------

fig.dose.last <- ggplot(df.dose.quantiles) +  
  geom_ribbon(aes(x = dose_total, 
                  ymin = log2(last_median_q5), 
                  ymax = log2(last_median_q95),
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, y = log2(last_median_median),
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 5, 1),
                     labels = 2^seq(-2, 5, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Days post infection") +
  facet_grid(route_name ~ "Time to undetectability") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.last


# Panel F: Duration ------------------------------------------------------------

fig.dose.duration <- ggplot(df.dose.quantiles) + 
  geom_ribbon(aes(x = dose_total, 
                  ymin = log2(duration_median_q5), 
                  ymax = log2(duration_median_q95),
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, y = log2(duration_median_median),
                color = as.character(tissue_idx)),
            alpha = 1, linewidth = 1) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-2, 5, 1),
                     labels = 2^seq(-2, 5, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Days") +
  facet_grid(route_name ~ "Duration of infection") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI),
                     labels = c("Nose", "Lung", "Lower GI")) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI),
                    labels = c("Nose", "Lung", "Lower GI")) +
  theme(text = element_text(size = 10),
        legend.position =  c(0.25, 0.8),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.dose.duration


# Legend -----------------------------------------------------------------------

fig.for.legend <- ggplot(df.dose.quantiles) +
  geom_point(aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = as.character(tissue_idx)), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.2,
             size = 1.5,
             color = "black", alpha = 1) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI),
                    labels = c("Nose", "Lung", "Lower GI")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = c("Nose", "Throat",
                              "Trachea", "Lung", 
                              "Upper GI",
                              "Lower GI")) +
  labs(y = "Probability of positivity (%)", fill = "Tissue") +
  coord_cartesian(clip = "on") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2))) +
  theme(text = element_text(size = 11),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.1, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "transparent", color = NA),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(1, 1, 10, 1),
        strip.placement = "outside",
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_blank()); fig.for.legend

# Extract the legend
g <- ggplotGrob(fig.for.legend)
legend <- g$grobs[[which(g$layout$name == "guide-box-top")]] # legend at the bottom
fig.legend.color <- ggdraw(legend); fig.legend.color # wrap as ggplot compatible item


# Combine ---------------------------------------------------------------

fig.dose.ed <-  fig.dose.prob.ci + labs(tag = "a") + 
  fig.dose.first.ci + labs(tag = "b") + 
  fig.dose.peak + labs(tag = "c") + 
  fig.dose.titer + labs(tag = "d") + 
  fig.dose.last + labs(tag = "e") + 
  (fig.dose.duration + labs(tag = "f") + theme(legend.position = "none")) + 
  plot_layout(nrow = 1); fig.dose.ed

fig.dose.ed <- (fig.legend.color / fig.dose.ed) + 
  plot_layout(heights = c(0.05, 1))


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS29-dose-effects-per-metric-all-routes-totalRNA.png",
       plot = fig.dose.ed,
       width = 12, 
       height = 8,
       dpi = 600)
