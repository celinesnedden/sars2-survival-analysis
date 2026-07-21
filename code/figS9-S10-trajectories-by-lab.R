# This file: - Compares predictions across labs for a given cofactor set

# Prep -------------------------------------------------------------------------

# Load the model fit
fit <- readRDS('./outputs/fits/fit-main.RDS')

# Load the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS")

# Get the unique assay-lab pairs and their routes & doses & tissue locations
df.groups <- data.frame(lab = dat.stan$lab,
                       assay = dat.stan$assay,
                       tissue = dat.stan$tissue_location,
                       organ = dat.stan$organ_location,
                       route = dat.stan$route,
                       dose_total = dat.stan$dose_total,
                       dose_nose = dat.stan$dose_nose,
                       dose_throat = dat.stan$dose_throat,
                       dose_trachea = dat.stan$dose_trachea,
                       dose_lung = dat.stan$dose_lung,
                       dose_gi = dat.stan$dose_gi,
                       age = dat.stan$age,
                       sex = dat.stan$sex,
                       species = dat.stan$species)

# Flag labs with unknown assay
labs_with_unknown <- unique(df.groups$lab[df.groups$assay == -9999])

# When various cofactors are unknown, we set to the most common
df.groups$assay[df.groups$assay == -9999] <- 1
df.groups$age[df.groups$age == -9999] <- 2
df.groups$sex[df.groups$sex == -9999] <- 1
df.groups$tissue[df.groups$tissue == -9999 & df.groups$organ == 1] <- 1
df.groups$tissue[df.groups$tissue == -9999 & df.groups$organ == 3] <- 6

# Get the distinct combinations
df.groups <- df.groups %>%
  distinct(lab, assay, tissue, route, dose_total, dose_nose, dose_throat,
           dose_trachea, dose_lung, dose_gi, age, sex, species)

num_rows <- dim(df.groups)[1]


# Fig S9 -----------------------------------------------------------------------

## Get predictions -------------------------------------------------------------

# Set the number of draws for each
n_draws <- 1000

# Generate predictions for the first lab, to set up the dataframe
df.pred.all <- get_metrics_for_specific_labs(fit, n_draws,
                                             seed = 5589,
                                             tissue = df.groups$tissue[1],
                                             dose_total = log10(df.groups$dose_total[1]),
                                             dose_nose = df.groups$dose_nose[1],
                                             dose_throat = df.groups$dose_throat[1],
                                             dose_trachea = df.groups$dose_trachea[1],
                                             dose_lung = df.groups$dose_lung[1],
                                             dose_gi = df.groups$dose_gi[1],
                                             route = df.groups$route[1],
                                             sex = df.groups$sex[1],
                                             species = df.groups$species[1], 
                                             age = df.groups$age[1],
                                             assay = df.groups$assay[1],
                                             lab = df.groups$lab[1])

# Loop over all the other rows
for (row_num in 2:dim(df.groups)[1]) {
  
  cat("Running row number ", row_num, "out of ", num_rows, "./")
  df.pred.next <- get_metrics_for_specific_labs(fit, n_draws,
                                                seed = 5589,
                                                tissue = df.groups$tissue[row_num],
                                                dose_total = log10(df.groups$dose_total[row_num]),
                                                dose_nose = df.groups$dose_nose[row_num],
                                                dose_throat = df.groups$dose_throat[row_num],
                                                dose_trachea = df.groups$dose_trachea[row_num],
                                                dose_lung = df.groups$dose_lung[row_num],
                                                dose_gi = df.groups$dose_gi[row_num],
                                                route = df.groups$route[row_num],
                                                sex = df.groups$sex[row_num],
                                                species = df.groups$species[row_num], 
                                                age = df.groups$age[row_num],
                                                assay = df.groups$assay[row_num],
                                                lab = df.groups$lab[row_num])
  
  df.pred.all <- rbind(df.pred.all, df.pred.next)
}

# Put all metrics on raw DPI scale not delay scale
df.pred.all$duration_median <- df.pred.all$peak_median + df.pred.all$last_median
df.pred.all$last_median <- df.pred.all$first_pos_median + df.pred.all$peak_median + df.pred.all$last_median
df.pred.all$peak_median <- df.pred.all$first_pos_median + df.pred.all$peak_median

# Get the medians for each lab from these samples
df.pred.medians <- df.pred.all %>%
  group_by(tissue_idx, assay_idx, lab_idx, dose_total, route_idx, sp_idx, age_idx, sex_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median), 
                   list(
                     median = ~median(.),
                     q_low = ~quantile(., probs = 0.025, na.rm = TRUE),
                     q_high = ~quantile(., probs = 0.975, na.rm = TRUE))))

# Set a group
df.pred.medians$group <- paste(df.pred.medians$lab_idx, df.pred.medians$tissue_idx, 
                               df.pred.medians$assay_idx, df.pred.medians$dose_total, 
                               df.pred.medians$route_idx, df.pred.medians$sp_idx,
                               df.pred.medians$age_idx, df.pred.medians$sex_idx, 
                               sep = "-")

# Flag unknown assays for plotting
df.pred.medians$assay_idx[df.pred.medians$lab_idx %in% labs_with_unknown] <- -9999

# Set all names
df.pred.medians <- assign_route_names(df.pred.medians)
df.pred.medians <- assign_tissue_names(df.pred.medians)
df.pred.medians <- assign_assay_names(df.pred.medians, long = FALSE)


## Get counts ------------------------------------------------------------------

df.medians.counts <- subset(df.pred.medians, tissue_idx != 5) %>%
  group_by(tissue_name, route_name, assay_name, assay_idx) %>%
  summarize(counts = n(), .groups = "drop")

df.medians.counts$assay_idx[df.medians.counts$assay_idx == -9999] <- 5


## Timing --------------------------------------------------------------------

fig.time <- ggplot(subset(df.pred.medians, tissue_idx != 5)) + 
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_name,
                   group = group), 
               alpha = 0.6) +
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_name,
                   #color = tissue_name,
                   group = group), alpha = 0.6) +
  geom_text(data = df.medians.counts,
            aes(x = 22, y = 9.8 - 1.2*assay_idx, color = assay_name,
                label = counts),
            size = 2.5, fontface = "bold") +
  labs(x = "Days post infection", y = "Viral titer (log10)", color = "Detection Assay") +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks = seq(0, 45, 6)) +
  coord_cartesian(clip = "off") +
  facet_grid(route_name ~ tissue_name) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.key.size = unit(0.7, "line"), 
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.time


## Percent --------------------------------------------------------------------

df.pred.medians$tissue_assay <- df.pred.medians$tissue_idx
df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 1] <- df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 1] - 0.3
df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 2] <- df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 2] - 0.15
df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 4] <- df.pred.medians$tissue_assay[df.pred.medians$assay_idx == 4] + 0.15
df.pred.medians$tissue_assay[df.pred.medians$assay_idx == -9999] <- df.pred.medians$tissue_assay[df.pred.medians$assay_idx == -9999] + 0.3

fig.percent <- ggplot(df.pred.medians) + 
  geom_point(aes(x = tissue_assay,
                 y = percent_positive_median * 100,
                 color = assay_name,
                 group = group), 
             alpha = 0.5) +
  scale_x_continuous(breaks = seq(1, 6, 1),
                     labels = c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI")) +
  labs(x = "Tissue Sampled", y = "Probability of positivity (%)") +
  coord_cartesian(clip = "off", ylim = c(0, 100)) +
  facet_grid(route_name ~ .) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_blank()); fig.percent 



## Legend --------------------------------------------------------------------

fig.for.legend <- ggplot(subset(df.pred.medians, tissue_idx != 5)) + 
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_name,
                   group = group), 
               alpha = 0.6) +
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_name,
                   #color = tissue_name,
                   group = group), alpha = 0.6) +
  labs(x = "Days post infection", y = "Viral titer (log10)", color = "Detection Assay") +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks = seq(0, 45, 6)) +
  facet_grid(route_name ~ tissue_name) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(text = element_text(size = 11),
        legend.position = "top",
        legend.key.size = unit(0.7, "line"), 
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.for.legend


fig.legend <- get_legend(fig.for.legend); fig.legend


## Combine ---------------------------------------------------------------------

fig.comb <- (as_ggplot(fig.legend)) + (fig.percent + labs(tag = "a") + 
  fig.time + labs(tag = "b") + 
  plot_layout(nrow = 1, widths = c(1, 4))) + 
  plot_layout(nrow = 2, heights = c(0.08, 1)); fig.comb


## Save ------------------------------------------------------------------------

ggsave("./outputs/figures/figS9-trajectories-by-lab-specific-cofactors.png",
       plot = fig.comb,
       width = 9, 
       height = 5.8,
       dpi = 600)



# Fig S10 -----------------------------------------------------------------------

## Get predictions -------------------------------------------------------------

# Set the number of draws for each
n_draws <- 1000

# Assuming the same dose distribution and route for all
dose_total <- 4
assign_dose_distribution(dose_total)
dose_nose <- get(paste0("route", 3, ".dose.nose")) / max_dose_nose
dose_throat <- get(paste0("route", 3, ".dose.throat")) / max_dose_throat
dose_trachea <- get(paste0("route", 3, ".dose.trachea")) / max_dose_trachea
dose_lung <- get(paste0("route", 3, ".dose.lung")) / max_dose_lung
dose_gi <- get(paste0("route", 3, ".dose.gi")) / max_dose_gi

# Generate predictions for the first lab, to set up the dataframe
df.pred.standard <- get_metrics_for_specific_labs(fit, n_draws,
                                                  seed = 5589,
                                                  tissue = df.groups$tissue[1],
                                                  dose_total = dose_total,
                                                  dose_nose = dose_nose,
                                                  dose_throat = dose_throat,
                                                  dose_trachea = dose_trachea,
                                                  dose_lung = dose_lung,
                                                  dose_gi = dose_gi,
                                                  route = 3,
                                                  sex = 1,
                                                  species = 1, 
                                                  age = 2,
                                                  assay = df.groups$assay[1],
                                                  lab = df.groups$lab[1])

# Loop over all the other rows
for (row_num in 2:dim(df.groups)[1]) {
  
  if (row_num %% 100 == 0) {
    cat("Running row number ", row_num, "out of ", num_rows, "./")
  }
  
  df.pred.next <- get_metrics_for_specific_labs(fit, n_draws,
                                                seed = 5589,
                                                tissue = df.groups$tissue[row_num],
                                                dose_total = dose_total,
                                                dose_nose = dose_nose,
                                                dose_throat = dose_throat,
                                                dose_trachea = dose_trachea,
                                                dose_lung = dose_lung,
                                                dose_gi = dose_gi,
                                                route = 3,
                                                sex = 1,
                                                species = 1, 
                                                age = 2,
                                                assay = df.groups$assay[row_num],
                                                lab = df.groups$lab[row_num])
  
  df.pred.standard <- rbind(df.pred.standard, df.pred.next)
}

# Put all metrics on raw DPI scale not delay scale
df.pred.standard$duration_median <- df.pred.standard$peak_median + df.pred.standard$last_median
df.pred.standard$last_median <- df.pred.standard$first_pos_median + df.pred.standard$peak_median + df.pred.standard$last_median
df.pred.standard$peak_median <- df.pred.standard$first_pos_median + df.pred.standard$peak_median

# Get the medians for each lab from these samples
df.pred.standard.medians <- df.pred.standard %>%
  group_by(tissue_idx, assay_idx, lab_idx, dose_total, route_idx, sp_idx, age_idx, sex_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median), 
                   list(
                     median = ~median(.),
                     q_low = ~quantile(., probs = 0.025, na.rm = TRUE),
                     q_high = ~quantile(., probs = 0.975, na.rm = TRUE))))

# Set a group
df.pred.standard.medians$group <- paste(df.pred.standard.medians$lab_idx, df.pred.standard.medians$tissue_idx, 
                               df.pred.standard.medians$assay_idx, df.pred.standard.medians$dose_total, 
                               df.pred.standard.medians$route_idx, df.pred.standard.medians$sp_idx,
                               df.pred.standard.medians$age_idx, df.pred.standard.medians$sex_idx, 
                               sep = "-")

# Flag unknown assays for plotting
df.pred.standard.medians$assay_idx[df.pred.standard.medians$lab_idx %in% labs_with_unknown] <- -9999

# Set all names
df.pred.standard.medians <- assign_route_names(df.pred.standard.medians)
df.pred.standard.medians <- assign_tissue_names(df.pred.standard.medians)
df.pred.standard.medians <- assign_assay_names(df.pred.standard.medians, long = FALSE)
#df.pred.standard.medians$assay_name <- factor(df.pred.standard.medians$assay_name, 
#                                              levels = rev(levels(df.pred.standard.medians$assay_name)))


## Get counts ------------------------------------------------------------------

df.medians.standard.counts <- subset(df.pred.standard.medians, tissue_idx != 5) %>%
  group_by(tissue_name, route_name, assay_name, assay_idx) %>%
  summarize(counts = n(), .groups = "drop")

df.medians.standard.counts$assay_idx[df.medians.standard.counts$assay_idx == -9999] <- 5


## Timing --------------------------------------------------------------------

fig.time.standard <- ggplot(subset(df.pred.standard.medians, tissue_idx != 5)) + 
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_name,
                   group = group), 
               alpha = 0.2) +
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_name,
                   group = group), alpha = 0.2) +
  geom_text(data = df.medians.standard.counts,
            aes(x = 15, y = 6, color = assay_name,
                label = counts),
            size = 2.5, fontface = "bold") +
  labs(x = "Days post infection", y = "Viral titer (log10)", 
       color = "Detection Assay") +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks = seq(0, 45, 6)) +
  coord_cartesian(clip = "off") +
  facet_grid(assay_name ~ tissue_name) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.key.size = unit(0.7, "line"), 
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.time.standard


## Percent --------------------------------------------------------------------

fig.percent.standard <- ggplot(subset(df.pred.standard.medians)) + 
  geom_point(aes(x = tissue_name,
                 y = percent_positive_median * 100,
                 color = assay_name,
                 group = group), 
             alpha = 0.2) +
  labs(x = "Tissue Sampled", y = "Probability of positivity (%)") +
  facet_grid(assay_name ~ .) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_blank()); fig.percent.standard


## Legend --------------------------------------------------------------------

fig.for.legend <- ggplot(subset(df.pred.medians, tissue_idx != 5)) + 
  geom_segment(aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_name,
                   group = group), 
               alpha = 0.6) +
  geom_segment(aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_name,
                   #color = tissue_name,
                   group = group), alpha = 0.6) +
  labs(x = "Days post infection", y = "Viral titer (log10)", color = "Detection Assay") +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks = seq(0, 45, 6)) +
  facet_grid(route_name ~ tissue_name) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(text = element_text(size = 11),
        legend.position = "top",
        legend.key.size = unit(0.7, "line"), 
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black")); fig.for.legend


fig.legend <- get_legend(fig.for.legend); fig.legend


## Combine ---------------------------------------------------------------------

fig.comb.standard <- (as_ggplot(fig.legend)) + 
  (fig.percent.standard + labs(tag = "a") + 
     fig.time.standard + labs(tag = "b") + 
     plot_layout(nrow = 1, widths = c(1, 4))) + 
  plot_layout(nrow = 2, heights = c(0.08, 1)); fig.comb.standard


## Save ------------------------------------------------------------------------

ggsave("./outputs/figures/figS10-trajectories-by-lab-fixed-cofactors.png",
       plot = fig.comb.standard,
       width = 9, 
       height = 5.8,
       dpi = 600)


