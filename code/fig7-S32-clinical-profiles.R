# This file:  - Compares AUC values across routes and 2 key doses


# Prep -------------------------------------------------------------------------

# Load the data
df.all <- fread( "./outputs/predictions/pred-across-cofactors-1000.csv")

# Subset to tissues of interest & demographic factors
df <- subset(df.all, assay_idx %in% c(1, 4) & tissue_idx %in% c(1, 4, 6) & 
                     sex_idx == 0 & age_idx == 2 & sp_idx == 1 &
                     dose_total %in% c(4, 7))

# Assign route & location names
df <- assign_route_names(df)
df <- assign_location_names(df)
df$route_name <- factor(df$route_name,
                        levels = rev(levels(df$route_name)))

# Remove columns that aren't needed
df.auc <- subset(df, select = -c(organ_group,
                                    organ_idx,
                                    dose_nose,
                                    dose_throat,
                                    dose_trachea,
                                    dose_lung,
                                    dose_gi,
                                    location_idx,
                                    lab_idx ,
                                    tissue_idx,
                                    percent_positive,
                                    first_pos_median,
                                    peak_median,
                                    titer_mean,
                                    last_median,
                                    duration_median))

# Expand to include one column per tissue location, conditional on all other cofactors
df.auc <- df.auc  %>% 
  pivot_wider(names_from = location_name,
              values_from = auc)

# Sample numbers that exist multiple times in the dataframe will create nested
#    entries in the corresponding tissue column, so we need to unnest them
df.auc <- unnest(df.auc, c('Nose', 'Lung', 'Lower GI'))

# Assign the route names
df.auc <- assign_route_names(df.auc)

# Calculate the medians per tissue, route, dose, and assay
df.auc.median <- df.auc %>%
  group_by(route_name, dose_total, assay_idx) %>%
  summarise(across(c(Nose, Lung, 'Lower GI'), ~median(.x, na.rm = TRUE)))

# Get the ratios with AE exposure at 10^7 dose for the NOSE
ae_nose <- df.auc.median$Nose[df.auc.median$route_name == "AE" & 
                                df.auc.median$assay_idx == 4 & 
                                df.auc.median$dose_total == 7]
df.auc.median$Nose_ratio <- df.auc.median$Nose / ae_nose

# Get the ratios with AE exposure at 10^7 dose for the GI
ae_gi <- df.auc.median$`Lower GI`[df.auc.median$route_name == "AE" & 
                                    df.auc.median$assay_idx == 4 & 
                                    df.auc.median$dose_total == 7]
df.auc.median$GI_ratio <- df.auc.median$`Lower GI` / ae_gi

# Get the ratios with AE exposure at 10^7 dose for the LUNG
ae_lung <- df.auc.median$Lung[df.auc.median$route_name == "AE" & 
                                df.auc.median$assay_idx == 4 & 
                                df.auc.median$dose_total == 7]
df.auc.median$Lung_ratio <- df.auc.median$`Lung` / ae_lung

# Combine into just a sample-by-sample dataframe
df.auc.samp <- df.auc %>%
  group_by(across(-c(Nose, Lung, `Lower GI`))) %>%
  summarise(
    Nose = na.omit(Nose)[1],
    Lung = na.omit(Lung)[1],
    GI   = na.omit(`Lower GI`)[1],
    .groups = "drop"
  )


# Set route names
df.auc.samp <- assign_route_names(df.auc.samp)

# Calculate the ratios based on the median AE
df.auc.samp$Nose_ratio <- df.auc.samp$Nose / ae_nose
df.auc.samp$GI_ratio <- df.auc.samp$GI / ae_gi
df.auc.samp$Lung_ratio <- df.auc.samp$Lung / ae_lung

# Significances ----------------------------------------------------------

# Subset to only extreme ends of the dose range
df.dose.full.sig <- subset(df.all, dose_total %in% c(1.2, 7.4) &
                             #dose_total %in% c(4, 7) & 
                             assay_idx == 4 &
                             tissue_idx %in% c(1, 4, 6))

# Set group 
df.dose.full.sig$group <- paste0(df.dose.full.sig$sample_num, "-", 
                                 df.dose.full.sig$tissue_idx)


# Quickly calculate differences for each metric
max_dose <- max(df.dose.full.sig$dose_total)
min_dose <- min(df.dose.full.sig$dose_total)

dt <- as.data.table(df.dose.full.sig)

dt[, `:=`(
  auc_diff      = ifelse(dose_total == max_dose,
                         auc -
                           auc[dose_total ==  min_dose],
                         NA_real_)
  
), by = .(sample_num, tissue_idx, route_idx, sex_idx, age_idx, sp_idx)]

df.dose.full.sig.diff <- as.data.frame(dt)

# Calculate quantiles of the differences
df.dose.diff.sig.quantiles <- df.dose.full.sig.diff %>%
  group_by(tissue_idx, route_idx) %>%
  summarise(across(c(auc_diff), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     q5 = ~quantile(., probs = 0.05, na.rm = TRUE),
                     q95 = ~quantile(., probs = 0.95, na.rm = TRUE))))


# Determine whether the difference is signficant (based on 90% CI)
df.sigs <- data.frame(tissue_idx = numeric(),
                      route_idx = numeric(),
                      metric = character(),
                      quantiles = character(),
                      signif = numeric(),
                      sign = numeric())
for (route.ii in 1:5) {
  for (tissue.ii in c(1, 4, 6)) {
    df.row <- subset(df.dose.diff.sig.quantiles, tissue_idx == tissue.ii & route_idx == route.ii)
    
    auc_signif <- sign(df.row$auc_diff_q5) * sign(df.row$auc_diff_q95)
    auc_sign <- sign(df.row$auc_diff_q5)
    
    
    
    df.add <- data.frame(tissue_idx = rep(tissue.ii, 1) ,
                         route_idx = rep(route.ii, 1) ,
                         metric = c("AUC"),
                         quantiles = c(paste0(df.row$auc_diff_q5, ", ", df.row$auc_diff_q95)),
                         signif = c(auc_signif),
                         sign = c(auc_sign),
                         size = c(as.numeric(df.row$auc_diff_median)))
    
    df.sigs <- rbind(df.sigs, df.add)
  }
}

df.sigs$signif[df.sigs$signif == -1] <- 0
df.sigs$signif[df.sigs$signif == 1 & df.sigs$sign == -1] <- -1

df.sigs$tissue_name[df.sigs$tissue_idx == 1] <- "Nose"
df.sigs$tissue_name[df.sigs$tissue_idx == 4] <- "Lung"
df.sigs$tissue_name[df.sigs$tissue_idx == 6] <- "Lower GI"

df.sigs$route_name[df.sigs$route_idx == 1] <- "IN"
df.sigs$route_name[df.sigs$route_idx == 2] <- "IT"
df.sigs$route_name[df.sigs$route_idx == 3] <- "IN+IT"
df.sigs$route_name[df.sigs$route_idx == 4] <- "AE"
df.sigs$route_name[df.sigs$route_idx == 5] <- "IG"
df.sigs$route_name <- factor(df.sigs$route_name,
                             levels = c("IN", "IT", "IN+IT",
                                        "AE", "IG"))


df.sigs$symbol[df.sigs$signif == "1"] <- "+"
df.sigs$symbol[df.sigs$signif == "-1"] <- "-"
df.sigs$symbol[df.sigs$signif == 0] <- ""


fig.sigs <- ggplot(df.sigs) + 
  geom_tile(aes(y = route_name, x = tissue_name), fill = "white", color = "black") + 
  geom_text(aes(y = route_name, x = tissue_name, label = symbol), 
            size = 3.1) +
  scale_x_discrete(limits = c("Nose", "Lung", "Lower GI"),
                   labels = c("Nasal  \nShedding", "Lung  \nSeverity", "Lower GI\nShedding"),
                   expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0),
                   limits = rev(c("IN", "IT", "IN+IT", "AE", "IG"))) +
  facet_wrap(.~ "Significance\nof Dose Effect") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()); fig.sigs



# Nose vs GI Shedding ------------------------------------------------------------

# Set coordinates of route labels
df.auc.median$Nose_label[df.auc.median$route_name == "IN"] <- 0.6
df.auc.median$GI_label[df.auc.median$route_name == "IN"] <- 0.5
df.auc.median$Nose_label[df.auc.median$route_name == "IN + IT"] <- 0.70
df.auc.median$GI_label[df.auc.median$route_name == "IN + IT"] <- 0.4
df.auc.median$Nose_label[df.auc.median$route_name == "AE"] <- 0.75
df.auc.median$GI_label[df.auc.median$route_name == "AE"] <- 0.8
df.auc.median$Nose_label[df.auc.median$route_name == "IT"] <- 0.4
df.auc.median$GI_label[df.auc.median$route_name == "IT"] <- 0.32
df.auc.median$Nose_label[df.auc.median$route_name == "IG"] <- 0.31
df.auc.median$GI_label[df.auc.median$route_name == "IG"] <- 0.37


figA <- ggplot(subset(df.auc.median, assay_idx == 4),
               aes(y = log2(Nose_ratio), x = log2(GI_ratio))) +
  geom_path(aes(color = route_name), alpha = 1) +
  geom_point(aes(fill = route_name, 
                 shape = as.character(dose_total)), 
             size = 3, alpha = 1) +
  geom_text(data = subset(df.auc.median, dose_total == 7 & assay_idx == 4),
            aes(y = log2(Nose_label), x = log2(GI_label), 
                label = route_name, color = route_name), 
            size = 3, fontface = "bold") +
  scale_shape_manual(values = c("4" = 21,
                                "7" = 24),
                     labels = c(expression(paste(10^4, " pfu")),
                                expression(paste(10^7, " pfu")))) +
  scale_y_continuous(#limits = c(-3, 1), 
                     breaks = seq(-5, 2, 0.5),
                     labels = round(2^seq(-5, 2, 0.5), digits = 2)) +
  scale_x_continuous(#limits = c(-2, 0), 
                     breaks = seq(-5, 2, 0.5),
                     labels = round(2^seq(-5, 2, 0.5), digits = 2)) +
  coord_cartesian(clip = 'off') +
  guides(color = "none", fill = "none") +
  labs(x = "GI Shedding\n(Relative AUC)", y = "Nasal Shedding\n(Relative AUC)",
       fill = "Exposure\nRoute", shape = "Exposure Dose") +
  theme(legend.position = c(0.75, 0.15),
        text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", 
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); figA


# Nose Shedding vs. Lung Severity --------------------------------------------

# Set coordinates of route labels
df.auc.median$Nose_label[df.auc.median$route_name == "IN"] <- 0.66
df.auc.median$Lung_label[df.auc.median$route_name == "IN"] <- 1.35
df.auc.median$Nose_label[df.auc.median$route_name == "IN + IT"] <- 0.70
df.auc.median$Lung_label[df.auc.median$route_name == "IN + IT"] <- 1.19
df.auc.median$Nose_label[df.auc.median$route_name == "AE"] <- 0.92
df.auc.median$Lung_label[df.auc.median$route_name == "AE"] <- 0.92
df.auc.median$Nose_label[df.auc.median$route_name == "IT"] <- 0.4
df.auc.median$Lung_label[df.auc.median$route_name == "IT"] <- 1.23
df.auc.median$Nose_label[df.auc.median$route_name == "IG"] <- 0.28
df.auc.median$Lung_label[df.auc.median$route_name == "IG"] <- 1.17


figB <- ggplot(subset(df.auc.median, assay_idx == 4),
               aes(y = log2(Nose_ratio), x = log2(Lung_ratio))) +
  geom_path(aes(color = route_name), alpha = 1) +
  geom_point(aes(fill = route_name, 
                 shape = as.character(dose_total)), 
             size = 3, alpha = 1) +
  geom_text(data = subset(df.auc.median, dose_total == 7 & assay_idx == 4),
            aes(x = log2(Lung_label), y = log2(Nose_label), label = route_name,
                color = route_name), size = 3, fontface = "bold") +
  scale_shape_manual(values = c("4" = 21,
                                "7" = 24),
                     labels = c(expression(paste(10^4, " pfu")),
                                expression(paste(10^7, " pfu")))) +
  scale_y_continuous(#limits = c(-3, 1), 
                     breaks = seq(-5, 2, 0.5),
                     labels = round(2^seq(-5, 2, 0.5), digits = 2)) +
  scale_x_continuous(#limits = c(-0.25, 0.52), 
                     breaks = seq(-2.25, 2.5, 0.25),
                     labels = round(2^seq(-2.25, 2.5, 0.25), digits = 2)) +
  coord_cartesian(clip = 'off') +
  guides(color = "none", fill = "none") +
  labs(x = "Lung Severity\n(Relative AUC)", 
       y = "Nasal Shedding\n(Relative AUC)",
       fill = "Exposure\nRoute", shape = "Exposure Dose") +
  theme(legend.position = "none",
        text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); figB


# Nose Density -----------------------------------------------------------------

df.auc.samp$density_color[df.auc.samp$dose_total == 7] <- "black"
df.auc.samp$density_color[df.auc.samp$dose_total == 4 &
                            df.auc.samp$route_idx == 1] <- "red"
df.auc.samp$density_color[df.auc.samp$dose_total == 4 &
                            df.auc.samp$route_idx == 2] <- "yellow"
df.auc.samp$density_color[df.auc.samp$dose_total == 4 &
                            df.auc.samp$route_idx == 3] <- "green"
df.auc.samp$density_color[df.auc.samp$dose_total == 4 &
                            df.auc.samp$route_idx == 4] <- "blue"
df.auc.samp$density_color[df.auc.samp$dose_total == 4 &
                            df.auc.samp$route_idx == 5] <- "purple"


df.auc.samp <- df.auc.samp %>% arrange(route_name)


fig.nose <- ggplot() +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 4),
                      aes(x = log2(Nose_ratio), 
                          y = route_name, fill = route_name,
                          color = route_name,
                          alpha = as.character(dose_total)),
                      rel_min_height = 0.01, scale = 1) +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 7),
                      aes(x = log2(Nose_ratio), 
                          y = route_name, fill = route_name,
                          alpha = as.character(dose_total)),
                      rel_min_height = 0.005, scale = 1) +
  scale_x_continuous(limits = c(-5, 1), 
                     breaks = seq(-5, 1, 1),
                     labels = c(".03", ".06", ".12", ".25", ".5", "1", "2")
                     #labels = round(2^seq(-5, 2, 1), digits = 2)
                     ) +
  scale_y_discrete(limits = rev(c("", "IN", "IT", "IN + IT", "AE", "IG")),
                   expand = c(0, 0.2)) +
  scale_alpha_manual(values = c("4" = 0.5, "7" = 1),
                     labels = c(expression(paste(10^4, " pfu")),
                                expression(paste(10^7, " pfu")))) +
  labs(x = "Nasal Shedding\n(Relative AUC)", alpha = "Dose") +
  guides(fill = "none", color = "none",
         alpha = guide_legend(override.aes = list(color = c("grey33", "black")))) +
  theme(legend.position = c(0.2, 0.8),
        text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "line"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.nose


# Lung Density --------------------------------------------------------

fig.lung <- ggplot() +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 7),
                      aes(x = log2(Lung_ratio), 
                          y = route_name, fill = route_name),
                      alpha = 1, rel_min_height = 0.01, scale = 1) +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 4),
                      aes(x = log2(Lung_ratio), 
                          y = route_name, fill = route_name,
                          color = route_name),
                      alpha = 0.5, scale = 1, 
                      rel_min_height = 0.01) +
  scale_x_continuous(limits = c(-3, 3), 
                     breaks = seq(-3, 3, 1),
                     labels = c(".125", ".25", ".5", "1", "2", "4", "8")) +
  scale_y_discrete(limits = rev(c("", "IN", "IT", "IN + IT", "AE", "IG")),
                   expand = c(0, 0.2)) +
  labs(x = "Lung Severity\n(Relative AUC)") +
  theme(legend.position = "none",
        text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.lung

# GI Density --------------------------------------------------------

fig.gi <- ggplot() +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 7),
                      aes(x = log2(GI_ratio), 
                          y = route_name, fill = route_name),
                      alpha = 0.9, rel_min_height = 0.01, scale = 1) +
  geom_density_ridges(data = subset(df.auc.samp, assay_idx == 4 & dose_total == 4),
                      aes(x = log2(GI_ratio), 
                          y = route_name, fill = route_name,
                          color = route_name),
                      alpha = 0.5, scale = 1,
                      rel_min_height = 0.01) +
  scale_x_continuous(#limits = c(-3, 3), 
                     breaks = seq(-3, 3, 1),
                     labels = c(".125", ".25", ".5", "1", "2", "4", "8")) +
  scale_y_discrete(limits = rev(c("", "IN", "IT", "IN + IT", "AE", "IG")),
                   expand = c(0, 0.2)) +
  labs(x = "GI Shedding\n(Relative AUC)") +
  theme(legend.position = "none",
        text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.gi


# Fig 6: Combine -------------------------------------------------------------------

fig7 <- figA + labs(tag = "A") + 
  figB + labs(tag = "B") +
  fig.sigs + labs(tag = "C") + 
  plot_layout(nrow = 1, widths = c(1, 1, 0.7)); fig7


ggsave("./outputs/figures/fig7-clinical-profiles.pdf",
       plot = fig7,
       width = 6.5, 
       height = 3.25,
       dpi = 600)



# Fig 32: Combine -------------------------------------------------------------------

fig.draws <- fig.nose + labs(tag = "a") +
  fig.lung + labs(tag = "b") + 
  fig.gi + labs(tag = "c"); fig.draws

ggsave("./outputs/figures/figS32-clinical-profiles-variability.png",
       plot = fig.draws,
       width = 8, 
       height = 2.2,
       dpi = 600)

