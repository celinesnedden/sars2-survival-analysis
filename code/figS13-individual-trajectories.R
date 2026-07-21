# This file:        - Creates individual time series (i.e., samples from full posteriors, 
#                       not just medians)


# Prep -------------------------------------------------------------------------

# Load color palette for plotting
df.color <- assign_colors()

# Load the model fit
fit <- readRDS('./outputs/fits/fit-main.RDS')

# Get individual samples from a fixed cofactor set
df.indiv <- get_individual_times(fit, 180,
                                 route_options = 3,
                                 tissue_options = c(1:4, 6),
                                 dose_options = 4,
                                 age_options = 2,
                                 species_options = 1,
                                 sex_options = 0,
                                 assay_options = 4)

# Convert to calendar time
df.indiv$last_sample <- df.indiv$first_sample + df.indiv$peak_sample + df.indiv$last_sample
df.indiv$peak_sample <- df.indiv$first_sample + df.indiv$peak_sample 


# Get samples of the median for the same fixed cofactor set
df.median <- get_metrics_across_cofactors(fit, 1000,
                                          route_options = 3,
                                          tissue_options = c(1:4, 6),
                                          dose_options = 4,
                                          age_options = 2,
                                          species_options = 1,
                                          sex_options = 0,
                                          assay_options = 4,
                                          rescale_doses = TRUE)

# Convert to calendar time
df.median$duration_median <- df.median$peak_median + df.median$last_median
df.median$last_median <- df.median$first_pos_median + df.median$peak_median + df.median$last_median
df.median$peak_median <- df.median$first_pos_median + df.median$peak_median

# Get individual samples from all covariate sets
df.indiv.all <- get_individual_times(fit, 1,
                                     route_options = 1:5,
                                     tissue_options = c(1:4, 6),
                                     dose_options = c(4, 7),
                                     age_options = 1:3,
                                     species_options = 1:3,
                                     sex_options = 0:1,
                                     assay_options = 4)

# Convert to calendar time
df.indiv.all$last_sample <- df.indiv.all$first_sample + df.indiv.all$peak_sample + df.indiv.all$last_sample
df.indiv.all$peak_sample <- df.indiv.all$first_sample + df.indiv.all$peak_sample 

# Assign all covariate names
df.indiv <- assign_all_names(df.indiv)
df.indiv.all <- assign_all_names(df.indiv.all)
df.median <- assign_all_names(df.median)

# Get the median of the medians
df.median.median <- df.median %>%
  group_by(route_name, dose_total, tissue_idx, tissue_name,
           assay_idx, assay_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 


# Panel A: Median of Medians + Median Samples ------------------------------

samp_nums <- sample(unique(df.median$sample_num), 180)

fig.median <- ggplot(subset(df.median, sample_num %in% samp_nums)) + 
  geom_segment(aes(x = first_pos_median, 
                   xend = peak_median,
                   y = 0,
                   yend = abs(titer_mean),
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.1) + 
  geom_segment(aes(x = peak_median, 
                   xend = last_median,
                   y = abs(titer_mean),
                   yend = 0,
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.1) +
  geom_segment(data = df.median.median,
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = abs(titer_mean_median),
                   color = as.character(tissue_idx)),
               alpha = 1, linewidth = 1) + 
  geom_segment(data = df.median.median,
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = abs(titer_mean_median),
                   yend = 0,
                   color = as.character(tissue_idx)),
               alpha = 1, linewidth = 1) +
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
  scale_y_continuous(breaks = seq(0, 12, 3))  +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 8)) +
  facet_grid("Samples of the inferred\nmedian for a fixed\ncofactor set" ~ tissue_name) +
  labs(x = "Days post infection", y = "Viral titer\n(log10 pfu)") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"), 
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = NA)); fig.median



# Panel B: Median of Medians + Individual Samples ------------------------------

fig.indivs <- ggplot(df.indiv) + 
  geom_segment(aes(x = first_sample, 
                   xend = peak_sample,
                   y = 0,
                   yend = abs(titer_sample),
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.2) + 
  geom_segment(aes(x = peak_sample, 
                   xend = last_sample,
                   y = abs(titer_sample),
                   yend = 0,
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.2) +
  geom_segment(data = df.median.median,
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = abs(titer_mean_median),
                   color = as.character(tissue_idx)),
               alpha = 1, linewidth = 1) + 
  geom_segment(data = df.median.median,
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = abs(titer_mean_median),
                   yend = 0,
                   color = as.character(tissue_idx)),
               alpha = 1, linewidth = 1) +
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
  scale_y_continuous(breaks = seq(0, 12, 3))  +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 8)) +
  facet_grid("Individual samples\nfor a fixed cofactor set" ~ tissue_name) +
  labs(x = "Days post infection", y = "Viral titer\n(log10 pfu)") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"), 
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.title.x = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = NA)); fig.indivs



# Panel C: Individual Samples across Covariates ------------------------------

fig.indiv.all <- ggplot(df.indiv.all) + 
  geom_segment(aes(x = first_sample, 
                   xend = peak_sample,
                   y = 0,
                   yend = abs(titer_sample),
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.2) + 
  geom_segment(aes(x = peak_sample, 
                   xend = last_sample,
                   y = abs(titer_sample),
                   yend = 0,
                   color = as.character(tissue_idx),
                   group = sample_num),
               alpha = 0.1, linewidth = 0.2) +
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
  scale_y_continuous(breaks = seq(0, 12, 3))  +
  coord_cartesian(xlim = c(0, 20), ylim = c(0, 8)) +
  facet_grid("Individual samples\nacross cofactor sets" ~ tissue_name) +
  labs(x = "Days post infection", y = "Viral titer") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"), 
        strip.placement = "outside",
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_blank(),
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
        strip.background = element_rect(fill = "white", color = NA)); fig.indiv.all



# Combine ----------------------------------------------------------------------

fig.comb <- fig.median + fig.indivs +  
  plot_layout(nrow = 2); fig.comb



# Save -------------------------------------------------------------------------

ggsave('./outputs/figures/figS13-individual-vs-median-trajectories.png',
       plot = fig.comb,
       width = 8, 
       height = 2.9,
       dpi = 600)

