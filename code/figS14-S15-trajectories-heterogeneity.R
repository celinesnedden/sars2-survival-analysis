# This file: - Generates time-series for all possible cofactor combinations
#            - Stratifies by whether the cofactor combination is included in the dataset used for fitting or extrapolated
#            - Computes the maximum difference in each metric for each tissue & detection assay among cofactor combinations


# Prep -------------------------------------------------------------------------

# Load the data passed to Stan 
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS")

# Load the predictions
df.pred <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")
df.pred <- subset(df.pred, dose_total %notin% c(1.2, 7.4))

# Assign group names 
df.pred$group <- paste0(df.pred$dose_total, df.pred$route_idx,
                        df.pred$sp_idx, df.pred$age_idx,
                        df.pred$sex_idx, #df.pred$assay_idx,
                        df.pred$sample_num, df.pred$tissue_idx)

# Assign location names
df.pred <- assign_tissue_names(df.pred)
df.pred$tissue_name <- factor(df.pred$tissue_name,
                                levels = rev(c("Nose", "Throat", "Trachea",
                                               "Lung", "Upper GI", "Lower GI")))
df.traj <- df.pred

# Create a column to distinguish among cofactor combinations (i.e., all cofactors)
df.traj$group <- paste(df.traj$dose_total, df.traj$route_idx,
                       df.traj$sp_idx, df.traj$age_idx,
                       df.traj$sex_idx, #df.traj$assay_idx,
                       df.traj$tissue_idx, sep = "-")

# Calculate medians & percentiles by group, tissue, and assay type
df.cof.medians <- df.traj %>% 
  group_by(tissue_idx, assay_idx, group, 
           dose_total, route_idx, sp_idx, age_idx, sex_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median), 
                   list(
                     median = ~median(.),
                     q_low = ~quantile(., probs = 0.025, na.rm = TRUE),
                     q_high = ~quantile(., probs = 0.975, na.rm = TRUE))))

# Set assay names
df.cof.medians <- assign_assay_names(df.cof.medians)

# Check which cofactor combinations exist in the database specifically
dat.stan$group <- paste(floor(log10(dat.stan$dose_total)),
                        dat.stan$route,
                        dat.stan$sp, dat.stan$age,
                        dat.stan$sex, #dat.stan$assay,
                        dat.stan$tissue_location, sep = "-")
groups_in_data <- unique(dat.stan$group[!str_detect(dat.stan$group, "9999")])

# Flag which groups are extrapolations
df.cof.medians$in_dataset <- 0
df.cof.medians$in_dataset[df.cof.medians$group %in% groups_in_data] <- 1


# Figure S14 -------------------------------------------------------------------

## Panel A: Culture, Probability -----------------------------------------------

# The probability of positivity
fig.percent <- ggplot() +
  geom_point(data = subset(df.cof.medians, assay_idx == 4 & in_dataset == 0),
             aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = "Not in dataset"), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 1.5,
             color = "black", alpha = 0.15) +
  geom_point(data = subset(df.cof.medians, assay_idx == 4 & in_dataset == 1),
               aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = as.character(tissue_idx)), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 1.5,
             color = "black", alpha = 0.8) +
  geom_point(data = subset(df.cof.medians, assay_idx == 4 & route_idx == 1 &
                             dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                             sp_idx == 1),
               aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 2,
             color = "black", fill = "black",
             alpha = 1) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "2" = df.color$Throat,
                               "3" = df.color$Trachea,
                               "4" = df.color$Lung,
                               "5" = df.color$Upper.GI,
                               "6" = df.color$Lower.GI,
                               "Not in dataset" = "grey")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = c("Nose", "Throat",
                              "Trachea", "Lung", 
                              "Upper GI",
                              "Lower GI")) +
  labs(y = "Probability of\nculture positivity (%)", x = "Tissue Sampled") +
  guides(fill = "none") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.percent


## Panel B: Culture, Trajectories ----------------------------------------------

# Trajectories
fig.traj <- ggplot(subset(df.cof.medians, assay_idx == 4 & tissue_idx != 5 & in_dataset == 1)) +
  geom_segment(data = subset(df.cof.medians, assay_idx == 4 & tissue_idx != 5 & in_dataset == 0),
                 aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   group = group),
               color = "grey", 
               alpha = 0.1, linewidth = 0.2) + 
  geom_segment(data = subset(df.cof.medians, assay_idx == 4 & tissue_idx != 5 & in_dataset == 0),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   group = group),
               color = "grey",
               alpha = 0.1, linewidth = 0.2) +
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = as.character(tissue_idx),
                   group = group),
               alpha = 10.8, linewidth = 0.2) + 
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = as.character(tissue_idx),
                   group = group),
               alpha = 0.8, linewidth = 0.2) +
  geom_segment(data = subset(df.cof.medians, assay_idx == 4 & route_idx == 1 &
                               dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                               sp_idx == 1 & tissue_idx != 5),
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   #color = as.character(tissue_idx),
                   group = group),
               color = "black",
               alpha = 0.8, linewidth = 0.5) + 
  geom_segment(data = subset(df.cof.medians, assay_idx == 4 & route_idx == 1 &
                               dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                               sp_idx == 1 & tissue_idx != 5),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   #color = as.character(tissue_idx),
                   group = group),
               color = "black",
               alpha = 1, linewidth = 0.5) +
  facet_wrap(.~ tissue_idx, nrow = 1) + 
  scale_color_manual(values = c("1" = df.color$Nose,
                                "2" = df.color$Throat,
                                "3" = df.color$Trachea,
                                "4" = df.color$Lung,
                                "5" = df.color$Upper.GI,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "2" = df.color$Throat,
                               "3" = df.color$Trachea,
                               "4" = df.color$Lung,
                               "5" = df.color$Upper.GI,
                               "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 16)) +
  labs(x = "Days post infection", y = "Viral titer\n(log10 pfu)") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"), 
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x = element_blank()); fig.traj



## Legend ----------------------------------------------------------------------

fig.for.legend <- ggplot(df.cof.medians) +
  geom_point(aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = as.character(tissue_idx),
                 alpha = as.character(in_dataset)), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.2,
             size = 1.5,
             color = "black") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "2" = df.color$Throat,
                               "3" = df.color$Trachea,
                               "4" = df.color$Lung,
                               "5" = df.color$Upper.GI,
                               "6" = df.color$Lower.GI),
                    labels = c("Nose", "Throat", "Trachea",
                               "Lung", "Upper GI", "Lower GI")) +
  scale_alpha_manual(values = c("0" = 0.2,
                                "1" = 1),
                    labels = c("No", "Yes")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = c("Nose", "Throat",
                              "Trachea", "Lung", 
                              "Upper GI",
                              "Lower GI")) +
  labs(y = "Probability of positivity (%)", fill = "Tissue", alpha = "Cofactor Combination In Dataset?") +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2)),
         alpha = guide_legend(nrow = 1, override.aes = list(size = 2, fill = c("grey", df.color$Nose),
                                                            alpha = 1))
         ) +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
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

fig.legend <- get_legend(fig.for.legend)


## Panel C: Total RNA, Probability ---------------------------------------------

fig.percent.RNA <- ggplot() +
  geom_point(data = subset(df.cof.medians, assay_idx == 1 & in_dataset == 0),
             aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = "Not in dataset"), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 1.5,
             color = "black", alpha = 0.15) +
  geom_point(data = subset(df.cof.medians, assay_idx == 1 & in_dataset == 1),
             aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100,
                 fill = as.character(tissue_idx)), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 1.5,
             color = "black", alpha = 0.8) +
  geom_point(data = subset(df.cof.medians, assay_idx == 1 & route_idx == 1 &
                             dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                             sp_idx == 1),
             aes(x = as.character(tissue_idx), 
                 y = percent_positive_median * 100), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.1,
             size = 2,
             color = "black", fill = "black",
             alpha = 1) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "2" = df.color$Throat,
                               "3" = df.color$Trachea,
                               "4" = df.color$Lung,
                               "5" = df.color$Upper.GI,
                               "6" = df.color$Lower.GI,
                               "Not in dataset" = "grey")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(labels = c("Nose", "Throat",
                              "Trachea", "Lung", 
                              "Upper GI",
                              "Lower GI")) +
  labs(y = "Probability of\ntotal RNA positivity (%)", x = "Tissue Sampled") +
  guides(fill = "none") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.percent.RNA


## Panel D: Total RNA, Trajectories --------------------------------------------

fig.traj.RNA <- ggplot(subset(df.cof.medians, assay_idx == 1 & tissue_idx != 5 & in_dataset == 1)) +
  geom_segment(data = subset(df.cof.medians, assay_idx == 1 & tissue_idx != 5 & in_dataset == 0),
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   group = group),
               color = "grey", 
               alpha = 0.1, linewidth = 0.2) + 
  geom_segment(data = subset(df.cof.medians, assay_idx == 1 & tissue_idx != 5 & in_dataset == 0),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   group = group),
               color = "grey",
               alpha = 0.1, linewidth = 0.2) +
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = as.character(tissue_idx),
                   group = group),
               alpha = 10.8, linewidth = 0.2) + 
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = as.character(tissue_idx),
                   group = group),
               alpha = 0.8, linewidth = 0.2) +
  geom_segment(data = subset(df.cof.medians, assay_idx == 1 & route_idx == 1 &
                               dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                               sp_idx == 1 & tissue_idx != 5),
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   #color = as.character(tissue_idx),
                   group = group),
               color = "black",
               alpha = 0.8, linewidth = 0.5) + 
  geom_segment(data = subset(df.cof.medians, assay_idx == 1 & route_idx == 1 &
                               dose_total == 4 & sex_idx == 0 & age_idx == 2 & 
                               sp_idx == 1 & tissue_idx != 5),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   #color = as.character(tissue_idx),
                   group = group),
               color = "black",
               alpha = 1, linewidth = 0.5) +
  facet_wrap(.~ tissue_idx, nrow = 1) + 
  scale_color_manual(values = c("1" = df.color$Nose,
                                "2" = df.color$Throat,
                                "3" = df.color$Trachea,
                                "4" = df.color$Lung,
                                "5" = df.color$Upper.GI,
                                "6" = df.color$Lower.GI)) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "2" = df.color$Throat,
                               "3" = df.color$Trachea,
                               "4" = df.color$Lung,
                               "5" = df.color$Upper.GI,
                               "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 16)) +
  labs(x = "Days post infection", y = "Viral titer\n(log10 total RNA copies)") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"), 
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x = element_blank()); fig.traj.RNA



## Combine ---------------------------------------------------------------------

figS14 <- (as_ggplot(fig.legend) + 
                (fig.percent + labs(tag = "a") +
                   fig.traj + labs(tag = "b") + plot_layout(nrow = 1, widths = c(1, 5))) + 
                (fig.percent.RNA + labs(tag = "c") +
                   fig.traj.RNA + labs(tag = "d") + plot_layout(nrow = 1, widths = c(1, 5)))) +
                  plot_layout(nrow = 3, heights = c(0.1, 1, 1)); figS14


## Save ------------------------------------------------------------------------

ggsave("./outputs/figures/figS14-trajectories-across-cofactor-combinations.png",
       plot = figS14,
       width = 9.5, 
       height = 4,
       dpi = 600)


# Figure S15 -------------------------------------------------------------------

## Prep -----------------------------------------------------------------------

# Calculate the maximum differences among all possible cofactor combos
df.diffs.all <- df.cof.medians %>%
  select(-contains("q_low"), -contains("q_high")) %>%
  pivot_longer(
    cols = -c(tissue_idx, assay_idx, assay_name, in_dataset, group, dose_total, route_idx, sp_idx, age_idx, sex_idx),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(tissue_idx, assay_name, assay_idx, metric) %>%
  summarise(max_diff = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
            .groups = "drop")

# Calculate the maximum differences among cofactor combos in data only
df.diffs.data <- subset(df.cof.medians, in_dataset == 1) %>%
  select(-contains("q_low"), -contains("q_high")) %>%
  pivot_longer(
    cols = -c(tissue_idx, assay_idx, assay_name, in_dataset, group, dose_total, route_idx, sp_idx, age_idx, sex_idx),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(tissue_idx, assay_name, assay_idx, metric) %>%
  summarise(max_diff = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
            .groups = "drop")

# Flag the computation type of each dataframe
df.diffs.all$type <- "All cofactor combinations"
df.diffs.data$type <- "Only cofactor combinations in dataset"

# Combine them into one dataframe for plotting
df.diffs <- rbind(df.diffs.all, df.diffs.data)

# Change the names of the metrics for plotting
df.diffs$metric[df.diffs$metric == "percent_positive_median"] <- "Probability of positivity"
df.diffs$metric[df.diffs$metric == "first_pos_median_median"] <- "Time to detectability"
df.diffs$metric[df.diffs$metric == "peak_median_median"] <- "Time to peak titer"
df.diffs$metric[df.diffs$metric == "titer_mean_median"] <- "Peak titer"
df.diffs$metric[df.diffs$metric == "last_median_median"] <- "Time to undetectability"
df.diffs$metric <- factor(df.diffs$metric, levels = c("Time to undetectability",
                                                      "Peak titer",
                                                      "Time to peak titer",
                                                      "Time to detectability",
                                                      "Probability of positivity"))

# Assign tissue names
df.diffs <- assign_tissue_names(df.diffs)

# Remove all upper GI metrics except positivity & time to detectability
df.diffs <- subset(df.diffs, !(tissue_idx == 5 & metric %notin% c("Probability of positivity", 
                                                                  "Time to detectability")))

# Rescale differences in positivity to out of 100
df.diffs$max_diff[df.diffs$metric == "Probability of positivity"] <- 100 * df.diffs$max_diff[df.diffs$metric == "Probability of positivity"]
df.diffs$max_diff[df.diffs$metric == "Probability of positivity"] <- round(df.diffs$max_diff[df.diffs$metric == "Probability of positivity"], 0)

## Plot ------------------------------------------------------------------------

fig.diffs <- ggplot(subset(df.diffs, assay_idx %in% c(1, 4)), 
                    aes(x = tissue_name, y = metric)) +
  geom_tile(aes(fill = tissue_name), alpha = 0.2, color = "black", linewidth = 0.25) + 
  geom_text(aes(label = round(max_diff, 1))) +
  labs(x = "Tissue Sampled", y = "Metric") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("Nose" = df.color$Nose,
                               "Throat" = df.color$Throat,
                               "Trachea" = df.color$Trachea,
                               "Lung" = df.color$Lung,
                               "Upper GI" = df.color$Upper.GI,
                               "Lower GI" = df.color$Lower.GI)) +
  coord_cartesian(clip = "off") +
  facet_grid(assay_name ~ type) +
  theme(text = element_text(size = 12),
        legend.position = "none",
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1)); fig.diffs


## Save ---------------------------------------------------------------------

ggsave("./outputs/figures/figS15-metric-heterogeneity-quantified.png",
       plot = fig.diffs,
       width = 7, 
       height = 4,
       dpi = 600)


