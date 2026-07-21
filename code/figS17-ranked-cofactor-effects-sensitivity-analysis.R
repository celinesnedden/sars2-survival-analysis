# This file: - Runs cofactor effect size comparisons, but using the median difference
#            - Recreates Figure 2 for this sensitivity analysis


# Load predictions & differences -----------------------------------------------

# Load and compile all the pairwise differences for the main doses
df.route.diffs <- readRDS("./outputs/predictions/pred-route-differences-1000.RDS")
df.dose.diffs <- readRDS("./outputs/predictions/pred-dose-differences-1000.RDS")
df.age.diffs <- readRDS("./outputs/predictions/pred-age-differences-1000.RDS")
df.sex.diffs <- readRDS("./outputs/predictions/pred-sex-differences-1000.RDS")
df.species.diffs <- readRDS("./outputs/predictions/pred-species-differences-1000.RDS")

# Load and compile all the pairwise differences for the lower dose range
df.route.diffs.low <- readRDS("./outputs/predictions/pred-route-differences-lowdoses-200.RDS")
df.dose.diffs.low <- readRDS("./outputs/predictions/pred-dose-differences-lowdoses-200.RDS")
df.age.diffs.low <- readRDS("./outputs/predictions/pred-age-differences-lowdoses-200.RDS")
df.sex.diffs.low <- readRDS("./outputs/predictions/pred-sex-differences-lowdoses-200.RDS")
df.species.diffs.low <- readRDS("./outputs/predictions/pred-species-differences-lowdoses-200.RDS")


# Panel A ---------------------------------------------------------------------

# Function to compare the mean differences in predictions
compute_median_differences <- function(df.diffs){
  df.diffs.mean <- df.diffs[0, ]
  
  for (tissue.ii in unique(df.diffs$tissue_idx)) {
    for (assay.ii in unique(df.diffs$assay_idx)) {
      df.sub <- subset(df.diffs, 
                       tissue_idx == tissue.ii &
                         assay_idx == assay.ii)
      
      
      new_row <- df.sub[1, ]
      new_row$diff_auc <- median(abs(df.sub$diff_auc))
      new_row$diff_percent <- median(abs(df.sub$diff_percent))
      new_row$diff_first <- median(abs(df.sub$diff_first))
      new_row$diff_titer <- median(abs(df.sub$diff_titer))
      new_row$diff_peak <- median(abs(df.sub$diff_peak))
      new_row$diff_last <- median(abs(df.sub$diff_last))
      new_row$diff_duration <- median(abs(df.sub$diff_duration))
      
      df.diffs.mean <- rbind(df.diffs.mean, new_row)
    }
  }
  return(df.diffs.mean)
}

# Apply function to all cofactors
df.route <- compute_median_differences(df.route.diffs)
df.dose <- compute_median_differences(df.dose.diffs)
df.age <- compute_median_differences(df.age.diffs)
df.sex <- compute_median_differences(df.sex.diffs)
df.species <- compute_median_differences(df.species.diffs)

# Combine
df.diffs.median <- rbind(df.route, df.dose, df.age, df.sex, df.species)


## Prep ---------------------------------------------------------------------

# Add tissue location names
df.rankings <- assign_tissue_names(df.diffs.median)

# Convert to long data frame (one metric per row)
df.rankings <- df.rankings %>% 
  pivot_longer(cols = starts_with("diff"),
               names_to = "metric",
               values_to = "difference")


# For a given metric, location, and tissue, rescale everything 
#     against the maximum effect
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

# Factor the cofactors for plotting
df.rankings$cofactor <- factor(df.rankings$cofactor, 
                               levels = rev(cof.sums$cofactor))

# Assing assay names
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


## Plot ---------------------------------------------------------------------

fig.culture.median <- 
  ggplot(subset(df.rankings, assay_idx == 4), 
         aes(y = factor(tissue_name, levels = rev(levels(tissue_name))), x = cofactor,
             fill = relative_effect)) +
  geom_tile(color = "black") +
  geom_tile(color = "black") +
  geom_text(aes(label = round(difference, 2)), size = 3 )+
  geom_tile(data = subset(df.rankings, assay_idx == 4 & rank == 5), 
            aes(x = cofactor, 
                y = factor(tissue_name, levels = rev(levels(tissue_name)))), 
            color = "black", linewidth = 0.8, fill = NA) +
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
                                    size=0.5)); fig.culture.median


# Panel B ---------------------------------------------------------------------

## Prep ---------------------------------------------------------------------

# Compute aggregate relative effects
cof.agg <- df.rankings %>%
  group_by(cofactor, assay_idx) %>%
  summarise(
    sum = sum(relative_effect))

cof.agg.relative <- cof.agg %>%
  group_by(assay_idx) %>%                     # Group by the metric column
  mutate(max_sum = max(sum),                  # Calculate the maximum sum for each metric
         scaled_sum = sum / max_sum) %>%      # Rescale the sum column
  ungroup() %>%                               # Ungroup the dataframe
  select(-max_sum)   


## Plot ---------------------------------------------------------------------

fig.agg.median <- ggplot(subset(cof.agg.relative, assay_idx %in% c(1, 4))) + 
  geom_bar(aes(y = scaled_sum, x = cofactor, 
               fill = as.character(assay_idx)), 
           stat = "identity", position = "dodge", color = "black") +
  labs(y = "Relative Effect", x = "Cofactor") +
  facet_wrap(.~ "All tissues & metrics") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("1" = "#F8766D",
                               "4" = "#00B0F6"),
                    limits = c("4", "1"),
                    labels = c("Culture",
                               "Total RNA")) +
  guides() +
  theme(
    text = element_text(size = 10),
    legend.position = c(0.28, 0.8),
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
                                size=0.5)); fig.agg.median


# Panel C ---------------------------------------------------------------------

## Prep ---------------------------------------------------------------------

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
df.route.low <- compute_mean_differences(df.route.diffs.low)
df.dose.low <- compute_mean_differences(df.dose.diffs.low)
df.age.low <- compute_mean_differences(df.age.diffs.low)
df.sex.low <- compute_mean_differences(df.sex.diffs.low)
df.species.low <- compute_mean_differences(df.species.diffs.low)

# Combine
df.diffs.mean <- rbind(df.route.low, df.dose.low, df.age.low, 
                       df.sex.low, df.species.low)

# Add tissue location names
df.rankings <- assign_tissue_names(df.diffs.mean)

# Convert to long data frame (one metric per row)
df.rankings <- df.rankings %>% 
  pivot_longer(cols = starts_with("diff"),
               names_to = "metric",
               values_to = "difference")


# Scale against maximum of the mean for each metric, location, and tissue
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
df.rankings<- assign_metric_names_units(df.rankings)

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



## Plot ---------------------------------------------------------------------

fig.culture.dose <- 
  ggplot(subset(df.rankings, assay_idx == 4), 
         aes(y = factor(tissue_name, levels = rev(levels(tissue_name))), x = cofactor,
             fill = relative_effect)) +
  geom_tile(color = "black") +
  geom_tile(color = "black") +
  geom_text(aes(label = round(difference, 2)), size = 3 ) +
  geom_tile(data = subset(df.rankings, assay_idx == 4 & rank == 5), 
            aes(x = cofactor, 
                y = factor(tissue_name, levels = rev(levels(tissue_name)))), 
            color = "black", linewidth = 0.8, fill = NA) +
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
                                    size=0.5)); fig.culture.dose


# Panel D --------------------------------------------------------------------

## Prep ------------------------------------------------------------------------

cof.agg <- df.rankings %>%
  group_by(cofactor, assay_idx) %>%
  summarise(
    sum = sum(relative_effect))

cof.agg.relative <- cof.agg %>%
  group_by(assay_idx) %>%                        # Group by the metric column
  mutate(max_sum = max(sum),                  # Calculate the maximum sum for each metric
         scaled_sum = sum / max_sum) %>%       # Rescale the sum column
  ungroup() %>%                               # Ungroup the dataframe
  select(-max_sum)   


## Plot ------------------------------------------------------------------------

fig.agg.dose <- ggplot(subset(cof.agg.relative, assay_idx %in% c(1, 4 ))) + 
  geom_bar(aes(y = scaled_sum, x = cofactor, 
               fill = as.character(assay_idx)), 
           stat = "identity", position = "dodge", color = "black") +
  #scale_fill_discrete(limits = c("route", "species", "dose", "age", "sex")) +
  labs(y = "Relative Effect", x = "Cofactor") +
  facet_wrap(.~ "All tissues & metrics") +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("1" = "#F8766D",
                               "4" = "#00B0F6"),
                    limits = c("4", "1"),
                    labels = c("Culture",
                               "Total RNA")) +
  guides() +
  theme(
    text = element_text(size = 10),
    legend.position = c(0.28, 0.8),
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
                                size=0.5)); fig.agg.dose


# Combine ----------------------------------------------------------------------

fig.supp.median <- fig.culture.median + labs(tag = "a") +
  ((fig.agg.median + labs(tag = "b")) / plot_spacer() / 
     plot_spacer() + plot_layout(heights = c(1, -0.156, 1))) + 
  plot_layout( widths = c(1, 0.25)); fig.supp.median

fig.supp.dose <- fig.culture.dose + labs(tag = "c") +
  ((fig.agg.dose + labs(tag = "d")) / plot_spacer() / 
     plot_spacer() + plot_layout(heights = c(1, -0.156, 1))) + 
  plot_layout( widths = c(1, 0.25)); fig.supp.dose

figS17 <- fig.supp.median  / fig.supp.dose; figS17



# Save ----------------------------------------------------------------------

ggsave('./outputs/figures/figS17-ranked-cofactor-sensitivity-analysis.png',
       plot = figS17,
       width = 10, 
       height = 8,
       dpi = 600) 
