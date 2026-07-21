# This file: - Shows predictions for each article in the database
#            - Generates example trajectories across cofactors
#            - Computes which cofactors have the biggest effect on each metric
#            - Computes which cofactor has the biggest effect overall


# Load predictions & differences -----------------------------------------------

# Load the predictions
df.pred <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Load and compile all the pairwise differences
df.route.diffs <- readRDS("./outputs/predictions/pred-route-differences-1000.RDS")
df.dose.diffs <- readRDS("./outputs/predictions/pred-dose-differences-1000.RDS")
df.age.diffs <- readRDS("./outputs/predictions/pred-age-differences-1000.RDS")
df.sex.diffs <- readRDS("./outputs/predictions/pred-sex-differences-1000.RDS")
df.species.diffs <- readRDS("./outputs/predictions/pred-species-differences-1000.RDS")


# Figure 3 ---------------------------------------------------------------------

## Panels A-B ------------------------------------------------------------------

### Prep lab-specific predictions  ----------------------------------------------

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


df.pred.medians$assay_new <- df.pred.medians$assay_idx
df.pred.medians$assay_new[df.pred.medians$assay_idx < 4] <- "PCR"
df.pred.medians$assay_new[df.pred.medians$assay_idx  == -9999] <- "PCR"
df.pred.medians$assay_new[df.pred.medians$assay_idx == 4] <- "Culture"
df.pred.medians$assay_new <- factor(df.pred.medians$assay_new,
                                    levels = c("Culture", "PCR"))

## Get counts 
df.medians.counts <- df.pred.medians %>%
  group_by(assay_new, tissue_name, tissue_idx) %>%
  summarize(counts = n(), .groups = "drop")

df.medians.counts$yaxis[df.medians.counts$assay_new == "PCR"] <- 8.5
df.medians.counts$yaxis[df.medians.counts$assay_new != "PCR"] <- 8.5 - 1


### B: Timing --------------------------------------------------------------------

fig.time <- ggplot() + 
  geom_segment(data = subset(df.pred.medians, tissue_idx != 5 & assay_name == "PCR"),
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_segment(data = subset(df.pred.medians, tissue_idx != 5 & assay_name != "PCR"),
               aes(x = first_pos_median_median, 
                   xend = peak_median_median,
                   y = 0,
                   yend = titer_mean_median,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_segment(data = subset(df.pred.medians, tissue_idx == 5 & assay_name == "PCR"), 
               aes(x = first_pos_median_median, 
                   xend = first_pos_median_median + 2,
                   y = 0,
                   yend = 2,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_segment(data = subset(df.pred.medians, tissue_idx == 5 & assay_name != "PCR"), 
               aes(x = first_pos_median_median, 
                   xend = first_pos_median_median + 2,
                   y = 0,
                   yend = 2,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_segment(data = subset(df.pred.medians, tissue_idx != 5 & assay_name == "PCR"),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_segment(data = subset(df.pred.medians, tissue_idx != 5 & assay_name != "PCR"),
               aes(x = peak_median_median, 
                   xend = last_median_median,
                   y = titer_mean_median,
                   yend = 0,
                   color = assay_new,
                   group = group), 
               linewidth = 0.1,
               alpha = 1) +
  geom_text(data = df.medians.counts,
            aes(x = 22, y = yaxis, color = assay_new,
                label = counts),
            size = 2.5, fontface = "bold") +
  scale_color_manual(values = c("PCR" = "grey55",
                                "Culture" = "#00B0F6")) +
  labs(x = "Days post infection", y = "Viral titer (log10)", color = "Detection Assay") +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  scale_x_continuous(breaks = seq(0, 45, 6)) +
  coord_cartesian(clip = "off") +
  facet_wrap(.~ tissue_name, nrow = 1) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(text = element_text(size = 11),
        legend.position = "top",
        legend.key.size = unit(0.7, "line"), 
        legend.key = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "white", color = "black", size = 0.5),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_text()); fig.time



### A: Percent --------------------------------------------------------------------

df.pred.medians$tissue_assay <- df.pred.medians$tissue_idx
df.pred.medians$tissue_assay[df.pred.medians$assay_new == "Culture"] <- df.pred.medians$tissue_assay[df.pred.medians$assay_new == "Culture"] + 0.15
df.pred.medians$tissue_assay[df.pred.medians$assay_new == "PCR"] <- df.pred.medians$tissue_assay[df.pred.medians$assay_new == "PCR"] - 0.15

fig.percent <- ggplot(df.pred.medians) + 
  geom_point(aes(x = tissue_assay,
                 y = percent_positive_median * 100,
                 color = assay_new,
                 group = group), 
             alpha = 0.5) +
  scale_x_continuous(breaks = seq(1, 6, 1),
                     labels = c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI")) +
  labs(x = "Tissue Sampled", y = "Probability of\npositivity (%)") +
  scale_color_manual(values = c("PCR" = "grey55",
                                "Culture" = "#00B0F6")) +
  coord_cartesian(clip = "off", ylim = c(0, 100)) +
  facet_wrap("" ~ .) +
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



### Combine panels A & B -------------------------------------------------------

fig.top <- fig.percent + labs(tag = "A") +
  fig.time + labs(tag = "B") + 
  plot_layout(nrow = 1, widths = c(1, 6)); fig.top



## Panel C ---------------------------------------------------------------------

### Prep -----------------------------------------------------------------------

# Function to compare the mean differences in predictions
compute_mean_differences <- function(df.diffs){
  df.diffs.mean <- df.diffs[0, ]
  
  for (tissue.ii in unique(df.diffs$tissue_idx)) {
    for (assay.ii in unique(df.diffs$assay_idx)) {
      df.sub <- subset(df.diffs, 
                         tissue_idx == tissue.ii &
                         assay_idx == assay.ii)
      
      
      new_row <- df.sub[1, ]
      new_row$diff_auc <- mean(abs(df.sub$diff_auc))
      new_row$diff_percent <- mean(abs(df.sub$diff_percent))
      new_row$diff_first <- mean(abs(df.sub$diff_first))
      new_row$diff_titer <- mean(abs(df.sub$diff_titer))
      new_row$diff_peak <- mean(abs(df.sub$diff_peak))
      new_row$diff_last <- mean(abs(df.sub$diff_last))
      new_row$diff_duration <- mean(abs(df.sub$diff_duration))
      
      df.diffs.mean <- rbind(df.diffs.mean, new_row)
    }
  }
  return(df.diffs.mean)
}

# Apply function to all cofactors
df.route <- compute_mean_differences(df.route.diffs)
df.dose <- compute_mean_differences(df.dose.diffs)
df.age <- compute_mean_differences(df.age.diffs)
df.sex <- compute_mean_differences(df.sex.diffs)
df.species <- compute_mean_differences(df.species.diffs)

# Combine
df.diffs.mean <- rbind(df.route, df.dose, df.age, df.sex, df.species)

# Assign tissue names
df.rankings <- assign_tissue_names(df.diffs.mean)

# Convert to long data frame (one metric per row)
df.rankings <- df.rankings %>% 
  pivot_longer(cols = starts_with("diff"),
               names_to = "metric",
               values_to = "difference")

# Scale against maximum of the means, for each metric, location, and tissue
df.rankings$relative_effect <- NA
for (metric.ii in unique(df.rankings$metric)) {
  for (loc.ii in unique(df.rankings$tissue_name)) {
    for (assay.ii in unique(df.rankings$assay_idx)) {
      auc.sub <- subset(df.rankings, metric == metric.ii &
                          assay_idx == assay.ii &
                          tissue_name == loc.ii)
      
      max_diff <- max(auc.sub$difference)
      df.rankings$relative_effect[df.rankings$metric == metric.ii &
                                    df.rankings$tissue_name == loc.ii &
                                    df.rankings$assay_idx == assay.ii] <- 
        df.rankings$difference[df.rankings$metric == metric.ii & 
                                 df.rankings$tissue_name == loc.ii &
                                 df.rankings$assay_idx == assay.ii] / max_diff
    }
    
  }
}

# Assign metric names
df.rankings <- assign_metric_names_units(df.rankings)

# Upper GI only has the first two metrics because it require invasive sampling
df.rankings <- subset(df.rankings, !(tissue_idx == 5 & 
                                       metric_name %notin% c("Probability of positivity (%)",
                                                             "Time to detectability (days)")))

# Find which rows have the largest summative effect
cof.sums <- data.frame(cofactor = character(),
                       sum = numeric())
for (cof.ii in unique(df.rankings$cofactor)) {
  auc.sub <- subset(df.rankings, cofactor == cof.ii & assay_idx == 4)
  
  cof.sums <- rbind(cof.sums, 
                    data.frame(cofactor = cof.ii,
                               sum = sum(auc.sub$relative_effect)))
}
cof.sums <- cof.sums[order(cof.sums$sum, decreasing = TRUE), ]

# Factor them for plotting
df.rankings$cofactor <- factor(df.rankings$cofactor, 
                               levels = rev(cof.sums$cofactor))

# Assign assay names
df.rankings <- assign_assay_names(df.rankings)

# Get the rankings
df.rankings <- df.rankings %>%
  group_by(tissue_name, metric_name, assay_name) %>%
  mutate(rank = rank(relative_effect))


# Change probability units to range from 0 to 100 (not 0 to 1)
df.rankings$difference[df.rankings$metric_name == "Probability of positivity (%)"] <-
  100 * df.rankings$difference[df.rankings$metric_name == "Probability of positivity (%)"]
df.rankings$difference[df.rankings$metric_name == "Probability of positivity (%)"] <-
  round(df.rankings$difference[df.rankings$metric_name == "Probability of positivity (%)"], 1)


### Plot -----------------------------------------------------------------------

fig.culture <- 
  ggplot(subset(df.rankings, assay_idx == 4), 
         aes(y = factor(tissue_name, levels = rev(levels(tissue_name))), x = cofactor,
             fill = relative_effect)) +
  geom_tile(color = "black") +
  geom_tile(data = subset(df.rankings, assay_idx == 4 & rank == 5), 
            aes(x = cofactor, 
                y = factor(tissue_name, levels = rev(levels(tissue_name)))), 
            color = "black", linewidth = 0.8, fill = NA) +
  geom_text(aes(label = round(difference, 2)), size = 3 )+
  facet_wrap(. ~ factor(metric_name, levels = rev(levels(metric_name))), ncol = 4) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", 
                      breaks = seq(0, 1, 0.5),
                      limits = c(0, 1)) +
  coord_cartesian(clip = "off") +
  labs(fill = "Relative\nEffect",
       x = "Cofactor", y = "Tissue Sampled",
       title = "Effects of cofactors on model predictions") +
  guides(fill = guide_colorbar(
    barwidth = 4,
    barheight = 0.8,
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    direction = "horizontal")) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35,  hjust = 1, vjust = 1),
        text = element_text(size = 10),
        legend.position = c(0.88, 0.2),
        legend.title = element_text(size = 9, vjust = 1),
        legend.text = element_text(size =8),
        legend.key = element_rect(color = "black", size = 10),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, 
                                    colour = "black",
                                    size=0.5)); fig.culture



## Panel D ---------------------------------------------------------------------

### Prep -----------------------------------------------------------------------

# Find aggregate effects
cof.agg <- df.rankings %>%
  group_by(cofactor, assay_idx) %>%
  summarise(
    sum = sum(relative_effect))

# Compute relative effects 
cof.agg.relative <- cof.agg %>%
  group_by(assay_idx) %>%                     # Group by the metric_name column
  mutate(max_sum = max(sum),                  # Calculate the maximum sum for each metric_name
         scaled_sum = sum / max_sum) %>%      # Rescale the sum column
  ungroup() %>%                               # Ungroup the dataframe
  select(-max_sum)   


### Plot -----------------------------------------------------------------------

fig.agg <- ggplot(subset(cof.agg.relative, assay_idx == 4)) + 
  geom_bar(aes(y = scaled_sum, x = cofactor, 
               fill = as.character(assay_idx)), 
           stat = "identity", position = "dodge", color = "black") +
  labs(y = "Relative effect across\ntissues and metrics", x = "Cofactor") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("1" = "#F8766D",
                               "4" = "#00B0F6"),
                    limits = c("4", "1"),
                    labels = c("Culture",
                               "Total RNA")) +
  guides() +
  theme(
    text = element_text(size = 10),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size =8),
    legend.key = element_rect(color = "white"),
    legend.key.size = unit(0.3, "cm"),
    legend.box.margin = margin(-15, 0, -10, 5),
    axis.text.x = element_text(angle = 35,  hjust = 1, vjust = 1),
    plot.margin = margin(-2, 0, 0, 0),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "white",
                                    colour = "black",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                    colour = "light grey"), 
    panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                    colour = "light grey"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(),
    legend.background = element_rect(fill = "transparent"),
    panel.border = element_rect(fill = NA, 
                                colour = "black",
                                size=0.5)); fig.agg


### Combine Panels C & D --------------------------------------------------------------

fig.bottom <- (fig.culture + labs(tag = "C")) + 
  (fig.agg + labs(tag = "D") + plot_spacer() + plot_layout(nrow = 2, heights = c(1, 0.5))) +
  plot_layout(ncol = 2, widths = c(4.5, 1)); fig.bottom


## Combine all panels together -------------------------------------------------

fig3 <- fig.top / fig.bottom + plot_layout(nrow = 2, heights = c(0.8, 2)); fig2


## Save ------------------------------------------------------------------------

ggsave('./outputs/figures/fig3-variability-across-cofactors.pdf',
       plot = fig2,
       width = 9.7, 
       height = 7,
       dpi = 600)


# Figure S16 -------------------------------------------------------------------

## Plot  -----------------------------------------------------------------------

figS16 <- 
  ggplot(subset(df.rankings, assay_idx == 4), 
         aes(y = factor(tissue_name, 
                        levels = rev(levels(tissue_name))), 
             x = cofactor,
             fill = relative_effect)) +
  geom_tile(color = "black") +
  geom_tile(data = subset(df.rankings, assay_idx == 4 & rank == 5), 
            aes(x = cofactor, 
                y = factor(tissue_name, levels = rev(levels(tissue_name)))), 
            color = "black", linewidth = 0.8, fill = NA) +
  geom_text(aes(label = round(relative_effect, 2)), size = 3 )+
  facet_wrap(. ~ factor(metric_name, levels = rev(levels(metric_name))), ncol = 4) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", 
                      breaks = seq(0, 1, 0.5),
                      limits = c(0, 1)) +
  coord_cartesian(clip = "off") +
  labs(fill = "Relative\nEffect",
       x = "Cofactor", y = "Tissue Sampled") +
  guides(fill = guide_colorbar(
    barwidth = 4,
    barheight = 0.8,
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    direction = "horizontal")) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35,  hjust = 1, vjust = 1),
        text = element_text(size = 10),
        legend.position = c(0.88, 0.2),
        legend.title = element_text(size = 9, vjust = 1),
        legend.text = element_text(size =8),
        legend.key = element_rect(color = "black", size = 10),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, 
                                    colour = "black",
                                    size=0.5)); figS16


## Save  -----------------------------------------------------------------------

ggsave('./outputs/figures/figS16-relative-cofactor-effects.png',
       plot = figS16,
       width = 7, 
       height = 3.7,
       dpi = 600) 

