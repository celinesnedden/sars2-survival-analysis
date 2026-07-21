# This file: - Determines whether demographic differences are significant
#            - Plots the metric-specific diffs & signifs between demographics


# Load predictions --------------------------------------------------------------

# Load full set of predictions
df <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Only include key tissues, culture assays, and mid-range dose
df <- subset(df, tissue_idx %in% c(1, 4, 6) & dose_total == 4 & assay_idx == 4)

# Assign all names
df <- assign_all_names(df)


# Compute differences among demographics  --------------------------------------

# Change duration column name
colnames(df)[colnames(df) == "duration_median"] <- "duration"

# Compute differences
df.sex.diffs <- get_pairwise_differences(df, "sex")
df.age.diffs <- get_pairwise_differences(df, "age")
df.species.diffs <- get_pairwise_differences(df, "species")
df.demog.sig.diff <- rbind(df.sex.diffs, df.age.diffs, df.species.diffs)

# Calculate quantiles 
df.demog.sig.diff.quantiles <- df.demog.sig.diff %>%
  group_by(cofactor, cofactor1, cofactor2, tissue_idx) %>%
  summarise(across(c(diff_percent, diff_first, 
                     diff_peak, diff_titer, diff_last,
                     diff_duration, diff_auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95))))


# Check for significance -------------------------------------------------------

# Data frame to fill in 
df.sigs <- data.frame(tissue_idx = numeric(),
                      cofactor = character(),
                      cofactor1 = numeric(),
                      cofactor2 = numeric(),
                      metric = character(),
                      signif = numeric(),
                      sign = numeric())

# Determine significance (Does 90% interval include 0?)
for (cof in unique(df.demog.sig.diff.quantiles$cofactor)) {
  df.sub1 <- subset(df.demog.sig.diff.quantiles, cofactor == cof)
  
  for (cof1 in unique(df.sub1$cofactor1)) {
    df.sub2 <- subset(df.sub1, cofactor1 == cof1)
    
    for (cof2 in unique(df.sub2$cofactor2)) {
      df.sub3 <- subset(df.sub2, cofactor2 == cof2)
      
      for (tissue.ii in unique(df.sub3$tissue_idx)){
        df.row <- subset(df.sub3, tissue_idx == tissue.ii)
        
        percent_signif <- sign(df.row$diff_percent_q5) * sign(df.row$diff_percent_q95)
        percent_sign <- sign(df.row$diff_percent_q5)
        
        first_signif <- sign(df.row$diff_first_q5) * sign(df.row$diff_first_q95)
        first_sign <- sign(df.row$diff_first_q5)
        
        peak_signif <- sign(df.row$diff_peak_q5) * sign(df.row$diff_peak_q95)
        peak_sign <- sign(df.row$diff_peak_q5)
        
        titer_signif <- sign(df.row$diff_titer_q5) * sign(df.row$diff_titer_q95)
        titer_sign <- sign(df.row$diff_titer_q5)
        
        last_signif <- sign(df.row$diff_last_q5) * sign(df.row$diff_last_q95)
        last_sign <- sign(df.row$diff_last_q5)
        
        duration_signif <- sign(df.row$diff_duration_q5) * sign(df.row$diff_duration_q95)
        duration_sign <- sign(df.row$diff_duration_q5)
        
        auc_signif <- sign(df.row$diff_auc_q5) * sign(df.row$diff_auc_q95)
        auc_sign <- sign(df.row$diff_auc_q5)
        
        df.add <- data.frame(tissue_idx = rep(tissue.ii, 7) ,
                             cofactor = rep(cof, 7),
                             cofactor1 = rep(cof1, 7),
                             cofactor2 = rep(cof2, 7),
                             metric = c("Probability of positivity", "Time to detectability",
                                        "Time to peak titer", "Peak titer",
                                        "Time to undetectability", "Duration of infection",
                                        "AUC"),
                             signif = c(as.numeric(percent_signif), 
                                        as.numeric(first_signif), 
                                        as.numeric(peak_signif),
                                        as.numeric(titer_signif), 
                                        as.numeric(last_signif), 
                                        as.numeric(duration_signif),
                                        as.numeric(auc_signif)),
                             sign = c(as.numeric(percent_sign), 
                                      as.numeric(first_sign), 
                                      as.numeric(peak_sign),
                                      as.numeric(titer_sign), 
                                      as.numeric(last_sign), 
                                      as.numeric(duration_sign),
                                      as.numeric(auc_sign)),
                             size = c(as.numeric(df.row$diff_percent_median),
                                      as.numeric(df.row$diff_first_median),
                                      as.numeric(df.row$diff_peak_median),
                                      as.numeric(df.row$diff_titer_median),
                                      as.numeric(df.row$diff_last_median),
                                      as.numeric(df.row$diff_duration_median),
                                      as.numeric(df.row$diff_auc_median)))
        
        df.sigs <- rbind(df.sigs, df.add)
        
      }
    }
  }
}

# Change key for significance
df.sigs$signif[df.sigs$signif == -1] <- 0
df.sigs$signif[df.sigs$signif == 1 & df.sigs$sign == -1] <- -1

df.sigs$symbol[df.sigs$signif == "1"] <- "+"
df.sigs$symbol[df.sigs$signif == "-1"] <- "-"
df.sigs$symbol[df.sigs$signif == 0] <- ""

# Add tissue ID numbers (for consistency when combining dfs)
df.sigs$tissue_idx[df.sigs$tissue_idx == "Nose"] <- 1
df.sigs$tissue_idx[df.sigs$tissue_idx == "Lung"] <- 4
df.sigs$tissue_idx[df.sigs$tissue_idx == "Lower GI"] <- 6

# Add cofactor numbers (for consistency when combining dfs)
df.demog$cofactor_num[df.demog$cofactor == "Species"] <- df.demog$sp_idx[df.demog$cofactor == "Species"] 
df.demog$cofactor_num[df.demog$cofactor == "Age"] <- df.demog$age_idx[df.demog$cofactor == "Age"] 
df.demog$cofactor_num[df.demog$cofactor == "Sex"] <- df.demog$sex_idx[df.demog$cofactor == "Sex"] 

# Combine into one database (for easier plotting purposes)
df.sigs$cofactor1_value <- NA
df.sigs$cofactor2_value <- NA
for (row_num in 1:nrow(df.sigs)) {
  
  df.demog.sub1 <- subset(df.demog, tissue_idx == df.sigs$tissue_idx[row_num] &
                            cofactor == str_to_title(df.sigs$cofactor[row_num]) &
                            cofactor_num == df.sigs$cofactor1[row_num])
  df.demog.sub2 <- subset(df.demog, tissue_idx == df.sigs$tissue_idx[row_num] &
                            cofactor == str_to_title(df.sigs$cofactor[row_num]) &
                            cofactor_num == df.sigs$cofactor2[row_num])
  
  
  if (df.sigs$metric[row_num] == "Probability of positivity") {
    val_cof1 <- as.numeric(df.demog.sub1$percent_positive_median)
    val_cof2 <- as.numeric(df.demog.sub2$percent_positive_median)
  }
  else if (df.sigs$metric[row_num] == "Time to detectability") {
    val_cof1 <- as.numeric(df.demog.sub1$first_pos_median_median)
    val_cof2 <- as.numeric(df.demog.sub2$first_pos_median_median)
  }
  else if (df.sigs$metric[row_num] == "Time to peak titer") {
    val_cof1 <- as.numeric(df.demog.sub1$peak_median_median)
    val_cof2 <- as.numeric(df.demog.sub2$peak_median_median)
  }
  else if (df.sigs$metric[row_num] == "Peak titer") {
    val_cof1 <- as.numeric(df.demog.sub1$titer_mean_median)
    val_cof2 <- as.numeric(df.demog.sub2$titer_mean_median)
  }
  else if (df.sigs$metric[row_num] == "Time to undetectability") {
    val_cof1 <- as.numeric(df.demog.sub1$last_median_median)
    val_cof2 <- as.numeric(df.demog.sub2$last_median_median)
  }
  else if (df.sigs$metric[row_num] == "Duration of infection") {
    val_cof1 <- as.numeric(df.demog.sub1$duration_median_median)
    val_cof2 <- as.numeric(df.demog.sub2$duration_median_median)
  }
  else if (df.sigs$metric[row_num] == "AUC") {
    val_cof1 <- as.numeric(df.demog.sub1$auc_median)
    val_cof2 <- as.numeric(df.demog.sub2$auc_median)
  }
  
  df.sigs$cofactor1_value[row_num] <- val_cof1
  df.sigs$cofactor2_value[row_num] <- val_cof2
}

# Set consistent acronym names for plotting
df.sigs$cof1_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor1 == 1] <- "RM"
df.sigs$cof1_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor1 == 2] <- "CM"
df.sigs$cof1_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor1 == 3] <- "AGM"
df.sigs$cof2_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor2 == 1] <- "RM"
df.sigs$cof2_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor2 == 2] <- "CM"
df.sigs$cof2_acronym[df.sigs$cofactor == "species" & df.sigs$cofactor2 == 3] <- "AGM"
df.sigs$cof1_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor1 == 1] <- "Juvenile"
df.sigs$cof1_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor1 == 2] <- "Adult"
df.sigs$cof1_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor1 == 3] <- "Geriatric"
df.sigs$cof2_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor2 == 1] <- "Juvenile"
df.sigs$cof2_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor2 == 2] <- "Adult"
df.sigs$cof2_acronym[df.sigs$cofactor == "age" & df.sigs$cofactor2 == 3] <- "Geriatric"
df.sigs$cof1_acronym[df.sigs$cofactor == "sex" & df.sigs$cofactor1 == 0] <- "Female"
df.sigs$cof1_acronym[df.sigs$cofactor == "sex" & df.sigs$cofactor1 == 1] <- "Male"
df.sigs$cof2_acronym[df.sigs$cofactor == "sex" & df.sigs$cofactor2 == 0] <- "Female"
df.sigs$cof2_acronym[df.sigs$cofactor == "sex" & df.sigs$cofactor2 == 1] <- "Male"

# Factor the acronyms for appropriate plotting order
df.sigs$cof1_acronym <- factor(df.sigs$cof1_acronym,
                               levels = c("RM", "CM", "AGM",
                                          "Juvenile", "Adult", "Geriatric",
                                          "Female", "Male"))
df.sigs$cof2_acronym <- factor(df.sigs$cof2_acronym,
                               levels = c("RM", "CM", "AGM",
                                          "Juvenile", "Adult", "Geriatric",
                                          "Female", "Male"))


# Plot each metric ------------------------------------------------------------

df.color <- assign_colors()

# Line type by significance
df.sigs$line[df.sigs$signif == 0] <- 0.25
df.sigs$line[df.sigs$signif == -1] <- 0.5
df.sigs$line[df.sigs$signif == 1] <- 0.5


## Probability of positivity  --------------------------------------------------

fig.demog.percent <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = percent_positive_median * 100,
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Probability of positivity"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = cofactor1_value * 100, yend = cofactor2_value * 100,
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  labs(y = "Probability (%)") +
  coord_cartesian(clip = "off") +
  facet_wrap(.~ "Probability of positivity") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.percent


## Time to detectability -------------------------------------------------------

fig.demog.first <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = log2(first_pos_median_median),
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Time to detectability"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = log2(cofactor1_value), yend = log2(cofactor2_value),
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  scale_y_continuous(breaks = seq(-1, 2, 1),
                     labels = 2^seq(-1, 2, 1)) +
  labs(y = "Days since inoculation") +
  facet_wrap(.~ "Time to detectability") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.first


## Time to peak titer ----------------------------------------------------------

fig.demog.peak <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = log2(peak_median_median),
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Time to peak titer"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = log2(cofactor1_value), yend = log2(cofactor2_value),
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  scale_y_continuous(breaks = seq(-1, 2, 1),
                     labels = 2^seq(-1, 2, 1)) +
  labs(y = "Days since inoculation") +
  facet_wrap(.~ "Time to peak titer") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.peak


## Time to undetectability -----------------------------------------------------

fig.demog.last <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = log2(last_median_median),
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Time to undetectability"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = log2(cofactor1_value), yend = log2(cofactor2_value),
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  scale_y_continuous(breaks = seq(-1, 7, 0.5),
                     labels = round(2^seq(-1, 7, 0.5), 1)) +
  labs(y = "Days since inoculation") +
  facet_wrap(.~ "Time to undetectability") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_text(size = 9),
        #strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.last


## Duration -------------------------------------------------------------------

fig.demog.duration <- ggplot(df.demog) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  geom_point(aes(x = as.numeric(cof_acronym), y = log2(duration_median_median),
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Duration of infection"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = log2(cofactor1_value), yend = log2(cofactor2_value),
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  scale_y_continuous(breaks = c(log2(5), log2(6), log2(7), log2(8), log2(10)),
                     labels = c(5, 6, 7, 8, 10)) +
  labs(y = "Days") +
  facet_wrap(.~ "Duration of infection") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.duration


## Peak titer -------------------------------------------------------------------

fig.demog.titer <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = titer_mean_median,
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  geom_segment(data = subset(df.sigs, metric == "Peak titer"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = cofactor1_value, yend = cofactor2_value,
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  labs(y = "Viral titer (log10 pfu)") +
  facet_wrap(.~ "Peak titer") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        #strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.titer


## AUC -------------------------------------------------------------------------

fig.demog.auc <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = auc_median,
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_vline(aes(xintercept = 3.5), linewidth = 0.6, color = "grey45") +
  geom_vline(aes(xintercept = 6.5), linewidth = 0.6, color = "grey45") +
  geom_segment(data = subset(df.sigs, metric == "AUC"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = cofactor1_value, yend = cofactor2_value,
                   linewidth = as.character(signif), color = as.character(tissue_idx),
                   linetype = as.character(signif))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8, 
                                    "-1" = 0.8)) +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid", 
                                   "-1" = "solid")) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI)) +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  labs(y = "AUC") +
  facet_wrap(.~ "AUC") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.demog.auc


## Legend -------------------------------------------------------------------

fig.for.legend <- ggplot(df.demog) +
  geom_point(aes(x = as.numeric(cof_acronym), y = titer_mean_median,
                 fill = as.character(tissue_idx)), shape = 21,
             size = 2.2, alpha = 0.8) +
  geom_segment(data = subset(df.sigs, metric == "Peak titer"),
               aes(x = as.numeric(cof1_acronym), xend = as.numeric(cof2_acronym),
                   y = cofactor1_value, yend = cofactor2_value,
                   linewidth = as.character(abs(signif)), color = as.character(tissue_idx),
                   linetype = as.character(abs(signif)))) +
  scale_linewidth_manual(values = c("0" = 0.5,
                                    "1" = 0.8), guide = "none") +
  scale_linetype_manual(values = c("0" = "dotted",
                                   "1" = "solid"),
                        labels = c("Insignificant", "Significant")) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI),
                    labels = c("Nose", "Lung", "Lower GI")) +
  scale_color_manual(values = c("1" = df.color$Nose,
                                "4" = df.color$Lung,
                                "6" = df.color$Lower.GI),
                     guide = "none") +
  scale_x_continuous(breaks = seq(1, max(as.numeric(df.demog$cof_acronym)), 1),
                     labels = levels(df.demog$cof_acronym)) +
  labs(y = "Viral titer (log10 pfu)", fill = "Tissue", linetype = "Difference") +
  facet_wrap(.~ "Peak titer") +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 10),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.for.legend

fig.legend <- as_ggplot(get_legend(fig.for.legend))


# Plot total counts ------------------------------------------------------------

# Calculate the number of significant differences among cofactor pairs
df.sigs.sum <- df.sigs %>%
  group_by(tissue_idx, cofactor, cofactor1, cofactor2) %>%
  summarise(sum_abs = sum(abs(signif))) 
df.sigs.sum.overall <- df.sigs %>%
  group_by(cofactor, cofactor1, cofactor2) %>%
  summarise(sum_abs = sum(abs(signif))) 

# Set acronyms for consistency
df.sigs.sum$cof1_acronym[df.sigs.sum$cofactor == "species" & df.sigs.sum$cofactor1 == "1"] <- "RM"
df.sigs.sum$cof1_acronym[df.sigs.sum$cofactor == "species" & df.sigs.sum$cofactor1 == "2"] <- "CM"
df.sigs.sum$cof2_acronym[df.sigs.sum$cofactor == "species" & df.sigs.sum$cofactor2 == "2"] <- "CM"
df.sigs.sum$cof2_acronym[df.sigs.sum$cofactor == "species" & df.sigs.sum$cofactor2 == "3"] <- "AGM"
df.sigs.sum$cof1_acronym[df.sigs.sum$cofactor == "age" & df.sigs.sum$cofactor1 == "1"] <- "Juvenile"
df.sigs.sum$cof1_acronym[df.sigs.sum$cofactor == "age" & df.sigs.sum$cofactor1 == "2"] <- "Adult"
df.sigs.sum$cof2_acronym[df.sigs.sum$cofactor == "age" & df.sigs.sum$cofactor2 == "2"] <- "Adult"
df.sigs.sum$cof2_acronym[df.sigs.sum$cofactor == "age" & df.sigs.sum$cofactor2 == "3"] <- "Geriatric"
df.sigs.sum$cof1_acronym[df.sigs.sum$cofactor == "sex" & df.sigs.sum$cofactor1 == "0"] <- "Female"
df.sigs.sum$cof2_acronym[df.sigs.sum$cofactor == "sex" & df.sigs.sum$cofactor2 == "1"] <- "Male"

# Labels for cofactor comparisons
df.sigs.sum$cof_comp <- paste0(df.sigs.sum$cof1_acronym, " vs. ", df.sigs.sum$cof2_acronym)
df.sigs.sum$cof_comp <- factor(df.sigs.sum$cof_comp,
                               levels = rev(c("RM vs. CM",
                                              "CM vs. AGM",
                                              "RM vs. AGM", 
                                              "Juvenile vs. Adult",
                                              "Adult vs. Geriatric",
                                              "Juvenile vs. Geriatric",
                                              "Female vs. Male")))

# Cofactor labels
df.sigs.sum.overall$cof1_acronym[df.sigs.sum.overall$cofactor == "species" & df.sigs.sum.overall$cofactor1 == "1"] <- "RM"
df.sigs.sum.overall$cof1_acronym[df.sigs.sum.overall$cofactor == "species" & df.sigs.sum.overall$cofactor1 == "2"] <- "CM"
df.sigs.sum.overall$cof2_acronym[df.sigs.sum.overall$cofactor == "species" & df.sigs.sum.overall$cofactor2 == "2"] <- "CM"
df.sigs.sum.overall$cof2_acronym[df.sigs.sum.overall$cofactor == "species" & df.sigs.sum.overall$cofactor2 == "3"] <- "AGM"

df.sigs.sum.overall$cof1_acronym[df.sigs.sum.overall$cofactor == "age" & df.sigs.sum.overall$cofactor1 == "1"] <- "Juvenile"
df.sigs.sum.overall$cof1_acronym[df.sigs.sum.overall$cofactor == "age" & df.sigs.sum.overall$cofactor1 == "2"] <- "Adult"
df.sigs.sum.overall$cof2_acronym[df.sigs.sum.overall$cofactor == "age" & df.sigs.sum.overall$cofactor2 == "2"] <- "Adult"
df.sigs.sum.overall$cof2_acronym[df.sigs.sum.overall$cofactor == "age" & df.sigs.sum.overall$cofactor2 == "3"] <- "Geriatric"

df.sigs.sum.overall$cof1_acronym[df.sigs.sum.overall$cofactor == "sex" & df.sigs.sum.overall$cofactor1 == "0"] <- "Female"
df.sigs.sum.overall$cof2_acronym[df.sigs.sum.overall$cofactor == "sex" & df.sigs.sum.overall$cofactor2 == "1"] <- "Male"

df.sigs.sum.overall$cof_comp <- paste0(df.sigs.sum.overall$cof1_acronym, " vs. ", df.sigs.sum.overall$cof2_acronym)

df.sigs.sum.overall$cof_comp <- factor(df.sigs.sum.overall$cof_comp,
                                       levels = rev(c("RM vs. CM",
                                                      "CM vs. AGM",
                                                      "RM vs. AGM", 
                                                      "Juvenile vs. Adult",
                                                      "Adult vs. Geriatric",
                                                      "Juvenile vs. Geriatric",
                                                      "Female vs. Male"
                                       )))

fig.sigs <- ggplot(df.sigs.sum) +
  geom_tile(aes(x = as.character(tissue_idx),
                y = cof_comp, fill = as.character(tissue_idx),
                alpha = sum_abs/max(df.sigs.sum$sum_abs)),
            color = "black") +
  geom_text(aes(x = as.character(tissue_idx),
                y = cof_comp, 
                label = sum_abs), size = 3.5) +
  geom_tile(data = df.sigs.sum.overall,
            aes(x = "7",
                y = cof_comp, fill = "grey",
                alpha = sum_abs/9),
            color = "black") +
  geom_text(data = df.sigs.sum.overall,
            aes(x = "7",
                y = cof_comp, 
                label = sum_abs), size = 3.5) +
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0), labels = c("Nose", "Lung", "Lower GI", "Overall")) +
  coord_cartesian(clip = "off") +
  facet_wrap(.~ "# significant differences") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.placement = "outside",
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); fig.sigs


# Combine ----------------------------------------------------------------------

# Combine panels in 2 rows
fig.demog <- fig.sigs + labs(tag = "a") +
  fig.demog.percent + labs(tag = "b") + 
  fig.demog.first + labs(tag = "c") + 
  fig.demog.peak +  labs(tag = "d") +
  fig.demog.titer + labs(tag = "e") + 
  fig.demog.last + labs(tag = "f") + 
  fig.demog.duration + labs(tag = "g") +
  fig.demog.auc + labs(tag = "h") +
  plot_layout(nrow = 2) ; fig.demog

# Add legend along the top
fig.combined <- (fig.legend / fig.demog ) + 
  plot_layout(heights = c(0.05, 1)); fig.combined 


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS19-demographic-effects-signifs.png",
       plot = fig.combined,
       width = 10, 
       height = 5,
       dpi = 600)

