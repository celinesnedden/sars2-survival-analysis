# This file: - Creates time-series stratified by demographic effects, integrating across all routes


# A. Average across routes -----------------------------------------------------

# Load predictions --------------------------------------------------------------

# Load full set of predictions
df <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Only include key tissues, culture assays, and mid-range dose
df <- subset(df, tissue_idx %in% c(1, 4, 6) & dose_total == 4 & assay_idx == 4)

# Assign all names
df <- assign_all_names(df)


# Compute median & percentiles -------------------------------------------------

## Species ---------------------------------------------------------------------

df.species <- df %>%
  group_by(sp_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 

## Age -------------------------------------------------------------------------

df.age <- df %>%
  group_by(age_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 


## Sex -------------------------------------------------------------------------

df.sex <- df %>%
  group_by(sex_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 


# Combine predictions -----------------------------------------------------------

# Add cofactor names & types for plotting
df.species$cofactor <- "Species"
df.age$cofactor <- "Age"
df.sex$cofactor <- "Sex"

# Combine, so they can all be plotted together in the same ggplot call
df.demog <- rbind(df.species, df.age, df.sex)

# Set acronyms for plotting ease
df.demog$cof_acronym[df.demog$sp_idx == 1] <- "RM"
df.demog$cof_acronym[df.demog$sp_idx == 2] <- "CM"
df.demog$cof_acronym[df.demog$sp_idx == 3] <- "AGM"
df.demog$cof_acronym[df.demog$age_idx == 1] <- "Juvenile"
df.demog$cof_acronym[df.demog$age_idx == 2] <- "Adult"
df.demog$cof_acronym[df.demog$age_idx == 3] <- "Geriatric"
df.demog$cof_acronym[df.demog$sex_idx == 0] <- "Female"
df.demog$cof_acronym[df.demog$sex_idx == 1] <- "Male"

df.demog$cof_acronym <- factor(df.demog$cof_acronym,
                               levels = c("RM", "CM", "AGM",
                                          "Juvenile", "Adult", "Geriatric",
                                          "Female", "Male"))

df.demog$group <- paste0(df.demog$tissue_idx, df.demog$cofactor)


# Plot  -----------------------------------------------------------------------

fig <- ggplot(df.demog, aes(color = cof_acronym, fill = cof_acronym)) + 
  geom_segment(aes(x = first_pos_median_median, xend = peak_median_median,
                   y = 0, yend = titer_mean_median),
               linewidth = 1, alpha = 0.9) +
  geom_segment(aes(x = peak_median_median, xend = last_median_median,
                   y = titer_mean_median, yend = 0),
               linewidth = 1, alpha = 0.9) +
  geom_col(aes(x = 12, y = 6 * percent_positive_median), 
           position = position_dodge(width = 2), width = 1.5,
           color = "black") +
  geom_text(aes(x = 9, y = 5, label = tissue_name), 
            color = "black", size = 3, alpha = 0.5) +
  facet_grid(factor(tissue_name, levels = rev(c("Lower GI", "Lung", "Nose"))) ~ 
               factor(cofactor, levels = c("Species", "Age", "Sex"))) +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",
       fill = "Cofactor") +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  scale_y_continuous(limits = c(0, 5.5)) +
  coord_cartesian(clip = "off", xlim = c(0, 13)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * 100/6, 
                                         breaks = seq(0, 100, 25),
                                         name = "Probability of positivity (%)"),
                     limits = c(0, 5.5)) +
  guides(fill = guide_legend(nrow = 3), color = "none") +
  scale_color_manual(values = c("RM" = "#F7766D",
                                "CM" = "#C77BFE",
                                "AGM" = "#00BE67", 
                                "Juvenile" = "#01BFC4",
                                "Adult" = "#FE61CB",
                                "Geriatric" = "#4062BB",
                                "Male" = "#00A9FF",
                                "Female" = "#DFAF2C")) +
  scale_fill_manual(values = c("RM" = "#F7766D",
                               "CM" = "#C77BFE",
                               "AGM" = "#00BE67", 
                               "Juvenile" = "#01BFC4",
                               "Adult" = "#FE61CB",
                               "Geriatric" = "#4062BB",
                               "Male" = "#00A9FF",
                               "Female" = "#DFAF2C")) +
  theme(text = element_text(size = 10),
        legend.position = c("bottom"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_blank(),
        strip.background = element_rect(fill = "white", color = "white")); fig


## Save  ------------------------------------------------------------------------
#
#ggsave("./outputs/figures/figS18-demographic-effects-all-routes.png",
#       plot = fig,
#       width = 6, 
#       height = 4.2,
#       dpi = 600)


# Load predictions --------------------------------------------------------------

# Load full set of predictions
df <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Only include aerosol individuals, key tissues, culture assays, and mid-range dose
df <- subset(df, tissue_idx %in% c(1, 4, 6) & dose_total == 4 & assay_idx == 4 & route_idx == 4)

# Assign all names
df <- assign_all_names(df)


# Compute median & percentiles -------------------------------------------------

## Species ---------------------------------------------------------------------

df.species <- df %>%
  group_by(sp_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 

## Age -------------------------------------------------------------------------

df.age <- df %>%
  group_by(age_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 


## Sex -------------------------------------------------------------------------

df.sex <- df %>%
  group_by(sex_idx, tissue_idx, tissue_name) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.025),
                     q95 = ~quantile(., probs = 0.975)))) 


# Combine predictions -----------------------------------------------------------

# Add cofactor names & types for plotting
df.species$cofactor <- "Species"
df.age$cofactor <- "Age"
df.sex$cofactor <- "Sex"

# Combine, so they can all be plotted together in the same ggplot call
df.demog <- rbind(df.species, df.age, df.sex)

# Set acronyms for plotting ease
df.demog$cof_acronym[df.demog$sp_idx == 1] <- "RM"
df.demog$cof_acronym[df.demog$sp_idx == 2] <- "CM"
df.demog$cof_acronym[df.demog$sp_idx == 3] <- "AGM"
df.demog$cof_acronym[df.demog$age_idx == 1] <- "Juvenile"
df.demog$cof_acronym[df.demog$age_idx == 2] <- "Adult"
df.demog$cof_acronym[df.demog$age_idx == 3] <- "Geriatric"
df.demog$cof_acronym[df.demog$sex_idx == 0] <- "Female"
df.demog$cof_acronym[df.demog$sex_idx == 1] <- "Male"

df.demog$cof_acronym <- factor(df.demog$cof_acronym,
                               levels = c("RM", "CM", "AGM",
                                          "Juvenile", "Adult", "Geriatric",
                                          "Female", "Male"))

df.demog$group <- paste0(df.demog$tissue_idx, df.demog$cofactor)

# Add age names
df.age <- df
df.age$cof_acronym <- df.age$age_idx
df.age$cof_acronym[df.age$cof_acronym == 1] <- "Juvenile"
df.age$cof_acronym[df.age$cof_acronym == 2] <- "Adult"
df.age$cof_acronym[df.age$cof_acronym == 3] <- "Geriatric"
df.age$cofactor <- "Age"

# Add species names
df.sp <- df
df.sp$cof_acronym <- df.sp$sp_idx
df.sp$cof_acronym[df.sp$cof_acronym == 1] <- "RM"
df.sp$cof_acronym[df.sp$cof_acronym == 2] <- "CM"
df.sp$cof_acronym[df.sp$cof_acronym == 3] <- "AGM"
df.sp$cofactor <- "Species"

# Add sex names
df.sex <- df
df.sex$cof_acronym <- df.sex$sex_idx
df.sex$cof_acronym[df.sex$cof_acronym == 0] <- "Female"
df.sex$cof_acronym[df.sex$cof_acronym == 1] <- "Male"
df.sex$cofactor <- "Sex"

# Combine together for plotting
df.draws <- rbind(df.age, df.sex, df.sp)

# Subset to only some trajectories to show
set.seed(555)
samps <- sample(unique(df.draws$sample_num), 100, replace = TRUE)
df.draws <- subset(df.draws, sample_num %in% samps)


# Plot  -----------------------------------------------------------------------

fig.ae <- ggplot(df.demog, aes(color = cof_acronym, fill = cof_acronym)) + 
  geom_segment(aes(x = first_pos_median_median, xend = peak_median_median,
                   y = 0, yend = titer_mean_median),
               linewidth = 1, alpha = 0.9) +
  geom_segment(aes(x = peak_median_median, xend = last_median_median,
                   y = titer_mean_median, yend = 0),
               linewidth = 1, alpha = 0.9) +
  geom_segment(data = df.draws, 
               aes(x = peak_median, xend = last_median,
                   y = titer_mean, yend = 0),
               linewidth = 0.05, alpha = 0.05) +
  geom_segment(data = df.draws, 
               aes(x = first_pos_median, xend = peak_median,
                   y = 0, yend = titer_mean),
               linewidth = 0.05, alpha = 0.05) +
  geom_col(aes(x = 12, y = 6 * percent_positive_median), 
           position = position_dodge(width = 2), width = 1.5,
           color = "black") +
  geom_text(aes(x = 9, y = 5.9, label = tissue_name), 
            color = "black", size = 3, alpha = 0.5) +
  facet_grid(factor(tissue_name, levels = rev(c("Lower GI", "Lung", "Nose"))) ~ 
               factor(cofactor, levels = c("Species", "Age", "Sex"))) +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",
       fill = "Cofactor") +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(clip = "off", xlim = c(0, 13)) +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . * 100/6, 
                                         breaks = seq(0, 100, 25),
                                         name = "Probability of positivity (%)")) +
  guides(fill = guide_legend(nrow = 3), color = "none") +
  scale_color_manual(values = c("RM" = "#F7766D",
                                "CM" = "#C77BFE",
                                "AGM" = "#00BE67", 
                                "Juvenile" = "#01BFC4",
                                "Adult" = "#FE61CB",
                                "Geriatric" = "#4062BB",
                                "Male" = "#00A9FF",
                                "Female" = "#DFAF2C")) +
  scale_fill_manual(values = c("RM" = "#F7766D",
                               "CM" = "#C77BFE",
                               "AGM" = "#00BE67", 
                               "Juvenile" = "#01BFC4",
                               "Adult" = "#FE61CB",
                               "Geriatric" = "#4062BB",
                               "Male" = "#00A9FF",
                               "Female" = "#DFAF2C")) +
  theme(text = element_text(size = 10),
        legend.position = c("bottom"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_blank(),
        strip.background = element_rect(fill = "white", color = "white")); fig.ae

# Combine ---------------------------------------------------------------------

fig.comb <- fig + labs(tag = "a") + 
  fig.ae + labs(tag = "b") + plot_layout(nrow = 1); fig.comb



# Save  ------------------------------------------------------------------------

ggsave("./outputs/figures/figS18-demographic-effects-all-routes.png",
       plot = fig.comb,
       width = 9, 
       height = 4.2,
       dpi = 600)
