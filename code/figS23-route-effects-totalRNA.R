# This file: - Plots the effects of route on all metrics, when monitoring via total RNA PCR 


# Prep predictions -------------------------------------------------------------

# Load predictions
df <- fread( "./outputs/predictions/pred-across-cofactors-1000.csv")

# Subset to dose of 10^4 pfu, total RNA data, from relevant tissues, & fixed cofactor set
df <- subset(df, dose_total == 4 & tissue_idx %in% c(1, 4, 6) & assay_idx == 1 &
               sp_idx == 1 & age_idx == 2 & sex_idx == 0)

# Set all cofactor names
df <- assign_all_names(df)

# Change order of route factors
df$route_name <- factor(df$route_name, levels = rev(levels(df$route_name)))

# Classify tissues as exposed or not
df$exposed <- "Not exposed"
df$exposed[df$route_idx %in% c(1, 3, 4) & df$tissue_name == "Nose"] <- "Exposed"
df$exposed[df$route_idx %in% c(1, 3, 4) & df$tissue_name == "Throat"] <- "Exposed"
df$exposed[df$route_idx %in% c(2, 3, 4) & df$tissue_name == "Trachea"] <- "Exposed"
df$exposed[df$route_idx %in% c(2, 3, 4) & df$tissue_name == "Lung"] <- "Exposed"
df$exposed[df$route_idx %in% c(5) & df$tissue_name == "Upper GI"] <- "Exposed"
df$exposed[df$route_idx %in% c(5) & df$tissue_name == "Lower GI"] <- "Exposed"

# Classify as URT exposed or not
df$urt_exposed <- "Not URT exposed"
df$urt_exposed[df$route_idx %in% c(1, 3, 4)] <- "URT exposed"

# Add combined tissue & exposed category
df$tissue_exposed <- paste(df$tissue_name, df$exposed)

# Get the quantiles for each route, dose, tissue, and assay type
df.quantiles <- df %>%
  group_by(route_name, exposed, urt_exposed, tissue_exposed, 
           dose_total, tissue_idx, tissue_name, organ_group,
           assay_idx, assay_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 


# Plot -------------------------------------------------------------------------

## Panel A: Probability --------------------------------------------------------

fig.percent <- ggplot(df) +
  geom_density_ridges(aes(x = percent_positive * 100, 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          fill = tissue_exposed,
                          linetype = urt_exposed),
                      scale = 1, rel_min_height = 0.002) +
  labs(x = "Probability (%)",
       y = "Exposure Routes") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, #"#5D3D73"
                               "Throat Exposed" = df.color$Throat,
                               "Trachea Exposed" = df.color$Trachea,
                               "Lung Exposed" = df.color$Lung, 
                               "Upper GI Exposed" = df.color$Upper.GI, 
                               "Lower GI Exposed" = df.color$Lower.GI, 
                               "Nose Not exposed" = df.color$Nose, #"#5D3D73"
                               "Throat Not exposed" = df.color$Throat,
                               "Trachea Not exposed" = df.color$Trachea,
                               "Lung Not exposed" = df.color$Lung, # df.color$Nose,
                               "Upper GI Not exposed" = df.color$Upper.GI,
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.45)) +
  scale_y_discrete(expand = expansion(add = c(.15, 0.9))) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 25)) +
  facet_wrap(.~ "Probability of positivity") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text.y = element_text(vjust = -0.3, size = 11),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(),
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
        strip.text = element_text(size = 11)); fig.percent


## Panel B: Detectability ------------------------------------------------------

fig.first <- ggplot(df) +
  geom_density_ridges(aes(x = log2(first_pos_median), 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          linetype = urt_exposed,
                          fill = tissue_exposed),
                      scale = 1, rel_min_height = 0.005) +
  labs(x = "Days post infection") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, 
                               "Lung Exposed" = df.color$Lung, 
                               "Lower GI Exposed" = df.color$Lower.GI,
                               "Nose Not exposed" = df.color$Nose, 
                               "Lung Not exposed" = df.color$Lung, 
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.45)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_y_discrete(expand = expansion(add = c(.15, 1))) +
  scale_x_continuous(limits = c(-3.5, 4), breaks = seq(-2, 4, 2),
                     labels = 2^seq(-2, 4, 2)) +
  facet_wrap(.~ "Time to detectability") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
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
        strip.text = element_text(size = 11)); fig.first



## Panel C: Peak time ----------------------------------------------------------

fig.peak <- ggplot(df) +
  geom_density_ridges(aes(x = log2(peak_median), 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          linetype = urt_exposed,
                          fill = tissue_exposed),
                      scale = 1, rel_min_height = 0.005) +
  labs(x = "Median days post infection", y = "Exposure Route") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, 
                               "Lung Exposed" = df.color$Lung, 
                               "Lower GI Exposed" = df.color$Lower.GI,
                               "Nose Not exposed" = df.color$Nose, 
                               "Lung Not exposed" = df.color$Lung, 
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.45)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_y_discrete(expand = expansion(add = c(.15, 1))) +
  scale_x_continuous(breaks = seq(0, 4, 1),
                     labels = 2^seq(0, 4, 1)) +
  facet_wrap(.~ "Time to peak titer") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.ticks.y = element_blank(),
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
        strip.text = element_text(size = 11)); fig.peak



## Panel D: Peak titer ---------------------------------------------------------

fig.titer <- ggplot(subset(df, !(route_idx == 5 & tissue_idx == 4))) +
  geom_density_ridges(aes(x = titer_mean, 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          linetype = urt_exposed,
                          fill = tissue_exposed),
                      scale = 1, rel_min_height = 0.005) +
  labs(x = "Mean viral titer (log10 pfu)") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, 
                               "Lung Exposed" = df.color$Lung, 
                               "Lower GI Exposed" = df.color$Lower.GI,
                               "Nose Not exposed" = df.color$Nose, 
                               "Lung Not exposed" = df.color$Lung, 
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.45)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_y_discrete(expand = expansion(add = c(.15, 1))) +
  #scale_x_continuous(breaks = seq(0, 10, 2),
  #                   labels = c(expression(paste(10^0)), 
  #                              expression(paste(10^2)), 
  #                              expression(paste(10^4)),  
  #                              expression(paste(10^6)),
  #                              expression(paste(10^8)),
  #                              expression(paste(10^10)))) +
  facet_wrap(.~ "Peak titer") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 11)); fig.titer



## Panel E: Undetectability ----------------------------------------------------

fig.last <- ggplot(df) +
  geom_density_ridges(aes(x = log2(last_median), 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          linetype = urt_exposed,
                          fill = tissue_exposed),
                      scale = 1, rel_min_height = 0.005) +
  labs(x = "Days post infection", y = "Exposure Route") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, 
                               "Lung Exposed" = df.color$Lung, 
                               "Lower GI Exposed" = df.color$Lower.GI,
                               "Nose Not exposed" = df.color$Nose, 
                               "Lung Not exposed" = df.color$Lung, 
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.4)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_y_discrete(expand = expansion(add = c(.15, 1))) +
  scale_x_continuous(breaks = seq(0, 5, 1),
                     labels = 2^seq(0, 5, 1)) +
  facet_wrap(.~ "Time to undetectability") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.ticks.y = element_blank(),
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
        strip.text = element_text(size = 11)); fig.last


## Panel F: Duration -----------------------------------------------------------

fig.duration <- ggplot(df) +
  geom_density_ridges(aes(x = log2(duration_median), 
                          y = route_name,
                          alpha = exposed,
                          color = exposed,
                          linetype = urt_exposed,
                          fill = tissue_exposed),
                      scale = 1, rel_min_height = 0.005) +
  labs(x = "Days", y = "Exposure Route") +
  scale_fill_manual(values = c("Nose Exposed" = df.color$Nose, 
                               "Lung Exposed" = df.color$Lung, 
                               "Lower GI Exposed" = df.color$Lower.GI,
                               "Nose Not exposed" = df.color$Nose, 
                               "Lung Not exposed" = df.color$Lung, 
                               "Lower GI Not exposed" = df.color$Lower.GI)) +
  scale_color_manual(values = c("Exposed" = "black",
                                "Not exposed" = "grey40")) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.4)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  scale_y_discrete(expand = expansion(add = c(.15, 1))) +
  scale_x_continuous(breaks = seq(0, 5, 1),
                     labels = 2^seq(0, 5, 1)) +
  facet_wrap(.~ "Duration of infection") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
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
        strip.text = element_text(size = 11)); fig.duration


## Panel G: URT exposed ---------------------------------------------------------------------

fig.traj.RNA.exposed <- 
  ggplot(df.quantiles) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   linetype = urt_exposed,
                   color = as.character(tissue_idx)),
               linewidth = 0.8, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   linetype = urt_exposed ,
                   color = as.character(tissue_idx)),
               linewidth = 0.8, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 2)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(clip = "off") +
  facet_wrap(factor(tissue_name, 
                    levels = c("Nose", "Lung", "Lower GI")) ~ .,
             scales = "free_x") +
  labs(x = "Days post infection", y = "Viral titer (log10)") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "2" = df.color$Throat,
                                "3" = df.color$Trachea,
                                "4" = df.color$Lung,
                                "5" = df.color$Upper.GI,
                                "6" = df.color$Lower.GI)) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash")) +
  guides(color = "none", linetype = guide_legend(override.aes = list(linewidth = 0.5))) +
  theme(text = element_text(size = 11),
        legend.position = c(0.18+0.33, 0.88),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_blank()); fig.traj.RNA.exposed



## Panel H: Trajectories -------------------------------------------------------

# Trajectories for each route in the nose, lung, and lower GI

# Create specific dataframe 
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

# Change factor order
df.traj.quantiles$route_name <- factor(df.traj.quantiles$route_name,
                                       levels = c("IN", "IT", "IN + IT", "AE", "IG"))

# Change factor order
df.traj$route_name <- factor(df.traj$route_name,
                             levels = rev(c("IN", "IT", "IN + IT", "AE", "IG")))

# Randomly sample some trajectories 
set.seed(4444)
random_samps <- sample(unique(df.traj$sample_num), size = 100)

# Plot
fig.traj.indiv.RNA <- 
  ggplot(subset(df.traj, sample_num %in% random_samps)) +
  geom_segment(aes(x = first_pos_median,
                   xend = peak_median,
                   y = 0, yend = titer_mean,
                   color = as.character(tissue_idx)),
               linewidth = 0.05, alpha = 0.15) +
  geom_segment(aes(x = peak_median,
                   xend = last_median,
                   yend = 0, y = titer_mean,
                   color = as.character(tissue_idx)),
               linewidth = 0.05, alpha = 0.15) +
  geom_segment(data = df.traj.quantiles,
               aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = as.character(tissue_idx)),
               linewidth = 0.5, alpha = 1) +
  geom_segment(data = df.traj.quantiles,
               aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = as.character(tissue_idx)),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 12), ylim = c(0, 9)) +
  facet_wrap(factor(route_name, levels = c("IN", "IT", "IN + IT", "AE", "IG")) ~., 
             ncol = 1, strip.position = "top") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)") +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
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
        strip.text = element_text(size = 11)); fig.traj.indiv.RNA


## Legend --------------------------------------------------------------

fig.for.legend <- ggplot(subset(df, tissue_idx %in% c(1, 4, 6))) +
  geom_point(aes(x = as.character(tissue_idx), 
                 y = percent_positive * 100,
                 fill = as.character(tissue_idx),
                 alpha = exposed), 
             shape = 21, position = position_jitter(w = 0.2, h = 0),
             stroke = 0.2,
             size = 1.5,
             color = "black") +
  geom_line(aes(x = as.character(tissue_idx), 
                y = percent_positive * 100,
                linetype = urt_exposed), 
            color = "black") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI),
                    labels = c("Nose", "Lung", "Lower GI")) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_alpha_manual(values = c("Exposed" = 0.95,
                                "Not exposed" = 0.45),
                     labels = c("Exposed tissue", "Not exposed tissue")) +
  scale_linetype_manual(values = c("URT exposed" = "solid",
                                   "Not URT exposed" = "twodash"),
                        limits = c("URT exposed", "Not URT exposed")) +
  scale_x_discrete(labels = c("Nose",  "Lung", "Lower GI")) +
  labs(y = "Probability of positivity (%)", fill = "",
       alpha = "", linetype = "") +
  coord_cartesian(clip = "off") +
  #guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2)),
  #       alpha = guide_legend(override.aes = list(size = 2, fill = "black"))) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3)),
         alpha = "none", linetype = "none") +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.5, "line"),
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
legend <- g$grobs[[which(g$layout$name == "guide-box-bottom")]] # legend at the bottom
fig.legend <- ggdraw(legend); fig.legend # wrap as ggplot compatible item




# Combine-----------------------------------------------------------------------

fig <- ((((fig.percent + fig.first + 
             fig.peak + theme(axis.title.x = element_text()) + labs(x = "Days post infection") +
             fig.titer + 
             fig.last  + 
             fig.duration +
             plot_layout(ncol = 2)) / fig.traj.RNA.exposed) + plot_layout(heights = c(4, 1))) |  
          fig.traj.indiv.RNA)  + plot_layout(widths = c(2.5, 1)) +
  plot_annotation(tag_levels = "a"); fig

fig <- cowplot::plot_grid(fig.legend, fig, ncol = 1, rel_heights = c(0.1, 2)); fig


# Save -------------------------------------------------------------------------

ggsave('./outputs/figures/figS23-route-effects-totalRNA.png',
       plot = fig,
       width = 7.5, 
       height = 10,
       dpi = 600)

