# This file: - Plots the main effects of dose


# Prep -------------------------------------------------------------------------

# Predictions across full dose suite
df.dose.full <- fread("./outputs/predictions/pred-across-doses-450.csv")
df.dose.fixed <- subset(df.dose.full, age_idx == 2 & sp_idx == 1 & sex_idx == 0)

# Data passed to stan for plotting sample sizes
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS")

# Set groups by sample number, tissue, and assay, for plotting
df.dose.fixed$group <- paste0(df.dose.fixed$sample_num, "-", 
                              df.dose.fixed$tissue_idx, 
                              "-", df.dose.fixed$assay_idx)

# Calculate medians & quantiles for plotting
df.dose.fixed.quantiles <- df.dose.fixed %>%
  group_by(dose_total, tissue_idx, route_idx, assay_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 

df.dose.full.quantiles <- df.dose.full %>%
  group_by(dose_total, tissue_idx, route_idx, assay_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
                     duration_median, auc), 
                   list(
                     median = ~median(.),
                     q5 = ~quantile(., probs = 0.05),
                     q95 = ~quantile(., probs = 0.95)))) 


# Plot -------------------------------------------------------------------------

## Panel A: ID50 per routes -----------------------------------------------------

### Compute ID50 ---------------------------------------------------------------

# Identify at which dose, we reach ~50% positivity 
df.id50 <- data.frame(route_idx = numeric(),
                      assay_idx = numeric(),
                      tissue_idx = numeric(),
                      id50 = numeric(),
                      diff50 = numeric())

for (route.ii in 1:5) {
  for (assay.ii in c(1, 4)) {
    for (tissue.ii in c(1, 4, 6)) {
      df.sub <- subset(df.dose.full.quantiles, route_idx == route.ii &
                         assay_idx == assay.ii & tissue_idx == tissue.ii)
      
      df.sub$diff50 <- 0.5 - df.sub$percent_positive_median
      
      row_closest <- which(abs(df.sub$diff50) == min(abs(df.sub$diff50)))
      
      df.id50 <- rbind(df.id50,
                       data.frame(route_idx = route.ii,
                                  assay_idx = assay.ii,
                                  tissue_idx = tissue.ii,
                                  closest_percent = df.sub$percent_positive_median[row_closest],
                                  id50 = df.sub$dose_total[row_closest],
                                  diff50 = min(abs(df.sub$diff50))))
      
    }
  }
}

# Set assay names
df.id50 <- assign_assay_names(df.id50)
df.id50$assay_name <- factor(df.id50$assay_name, 
                             levels = c("Culture", "Total RNA"))

# Set route names
df.id50$route_name[df.id50$route_idx == 1] <- "IN"
df.id50$route_name[df.id50$route_idx == 2] <- "IT"
df.id50$route_name[df.id50$route_idx == 3] <- "IN+IT"
df.id50$route_name[df.id50$route_idx == 4] <- "AE"
df.id50$route_name[df.id50$route_idx == 5] <- "IG"
df.id50$route_name <- factor(df.id50$route_name,
                             levels = c("IN", "IT", "IN+IT",
                                        "AE", "IG"))

# Flag exposed tissues
df.id50$exposed <- "Not exposed"
df.id50$exposed[df.id50$route_idx %in% c(1, 3, 4) & df.id50$tissue_idx == 1] <- "Exposed"
df.id50$exposed[df.id50$route_idx %in% c(2, 3, 4) & df.id50$tissue_idx == 4] <- "Exposed"
df.id50$exposed[df.id50$route_idx %in% c(5) & df.id50$tissue_idx == 6] <- "Exposed"

# Determine if any are upper or lower bounds
df.id50$id50_exp <- df.id50$id50
df.id50$id50_exp[df.id50$id50 == 1.2 & df.id50$diff50 >= 0.1] <- 
  paste0('<', df.id50$id50_exp[df.id50$id50 == 1.2 & df.id50$diff50 >= 0.1])
df.id50$id50_exp[df.id50$id50 == 7.4 & df.id50$diff50 >= 0.1] <- 
  paste0(">", df.id50$id50_exp[df.id50$id50 == 7.4 & df.id50$diff50 >= 0.1])


### Plot ---------------------------------------------------------------

fig.id50 <- ggplot(df.id50) + 
  geom_tile(aes(y = as.character(tissue_idx), x = route_name,
                fill = as.character(tissue_idx), alpha = exposed),
            linewidth = 0.4, color = "black") +
  geom_text(aes(y = as.character(tissue_idx), x = route_name,
                label = id50_exp), parse = F, size = 3.1) +
  scale_alpha_manual(values = c("Exposed" = 0.85, "Not exposed" = 0.15)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev(c("1", "4", "6"))) + 
  scale_fill_manual(values = c("1" = df.color$Nose,
                               "4" = df.color$Lung,
                               "6" = df.color$Lower.GI)) +
  facet_grid(assay_name ~ "Median ID50 (log10)") +
  labs(y = "", x = "Exposure Route") +
  theme(text = element_text(size = 10),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.key.size = unit(0.7, "line"),
        axis.title = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5)); fig.id50


## Panel B: Culture probability ------------------------------------------------

# Subset Stan data to just AE
dat.stan.ae <- subset(as.data.frame(dat.stan), route == 4)

# Set groups for plotting consistency
df.dose.fixed.quantiles$group <- paste0(df.dose.fixed.quantiles$tissue_idx,
                                        df.dose.fixed.quantiles$route_idx,
                                        df.dose.fixed.quantiles$assay_idx)

# Plot
fig.dose.percent.ci <- ggplot(subset(df.dose.fixed.quantiles, route_idx == 4 & assay_idx == 4)) +
  geom_histogram(data = dat.stan.ae,
                 aes(x = log10(dose_AE)), color = "grey44", fill = "grey77",
                 breaks = seq(1.2, 5, 0.2)) + 
  geom_ribbon(aes(x = dose_total, ymin = percent_positive_q5 * 100, 
                  ymax = percent_positive_q95 * 100,
                  fill = as.character(tissue_idx)), 
              alpha = 0.2) +
  geom_line(aes(x = dose_total, y = percent_positive_median * 100,
                color = as.character(tissue_idx), group = group),
            alpha = 1, linewidth = 1) +
  geom_hline(aes(yintercept = 50), linetype = "dotted") +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c(expression(paste(10^1)), 
                                expression(paste(10^2)),
                                expression(paste(10^3)),
                                expression(paste(10^4)), 
                                expression(paste(10^5)),
                                expression(paste(10^6)),
                                expression(paste(10^7))),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 25),
                     limits = c(0, 100)) + #,
  coord_cartesian(clip = "off") +
  labs(x = "Exposure dose (pfu)", y = "Probability (%)") +
  facet_wrap(.~ "Probability of culture positivity") +
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
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        axis.title = element_text(size = 9),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5)); fig.dose.percent.ci



## Panel C: Significances -------------------------------------------------------

# Subset to only extreme ends of the dose range
df.dose.full.sig <- subset(df.dose.full, dose_total %in% c(1.2, 7.4) & assay_idx == 4)

# Set group 
df.dose.full.sig$group <- paste0(df.dose.full.sig$sample_num, "-", 
                                 df.dose.full.sig$tissue_idx)


# Quickly calculate differences for each metric
dt <- as.data.table(df.dose.full.sig)

dt[, `:=`(
  percent_diff  = ifelse(dose_total == 7.4,
                         percent_positive * 100 -
                           percent_positive[dose_total == 1.2] * 100,
                         NA_real_),
  
  first_diff    = ifelse(dose_total == 7.4,
                         first_pos_median -
                           first_pos_median[dose_total == 1.2],
                         NA_real_),
  
  peak_diff     = ifelse(dose_total == 7.4,
                         peak_median -
                           peak_median[dose_total == 1.2],
                         NA_real_),
  
  titer_diff    = ifelse(dose_total == 7.4,
                         titer_mean -
                           titer_mean[dose_total == 1.2],
                         NA_real_),
  
  last_diff     = ifelse(dose_total == 7.4,
                         last_median -
                           last_median[dose_total == 1.2],
                         NA_real_),
  
  duration_diff = ifelse(dose_total == 7.4,
                         duration_median -
                           duration_median[dose_total == 1.2],
                         NA_real_),
  
  auc_diff      = ifelse(dose_total == 7.4,
                         auc -
                           auc[dose_total == 1.2],
                         NA_real_)
  
), by = .(sample_num, tissue_idx, route_idx, sex_idx, age_idx, sp_idx)]

df.dose.full.sig.diff <- as.data.frame(dt)


# Convert to a long dataframe
df.dose.full.sig.diff.long <- df.dose.full.sig.diff %>%
  pivot_longer(percent_diff:last_diff,
               names_to = "metric",
               values_to = "diff") 
df.dose.full.sig.diff.long <- subset(df.dose.full.sig.diff.long, !is.na(diff))

# Calculate quantiles of the differences
df.dose.diff.sig.quantiles <- df.dose.full.sig.diff %>%
  group_by(tissue_idx, route_idx) %>%
  summarise(across(c(percent_diff, first_diff, 
                     peak_diff, titer_diff, last_diff,
                     duration_diff, auc_diff), 
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
    
    percent_signif <- sign(df.row$percent_diff_q5) * sign(df.row$percent_diff_q95)
    percent_sign <- sign(df.row$percent_diff_q5)
    
    first_signif <- sign(df.row$first_diff_q5) * sign(df.row$first_diff_q95)
    first_sign <- sign(df.row$first_diff_q5)
    
    peak_signif <- sign(df.row$peak_diff_q5) * sign(df.row$peak_diff_q95)
    peak_sign <- sign(df.row$peak_diff_q5)
    
    titer_signif <- sign(df.row$titer_diff_q5) * sign(df.row$titer_diff_q95)
    titer_sign <- sign(df.row$titer_diff_q5)
    
    last_signif <- sign(df.row$last_diff_q5) * sign(df.row$last_diff_q95)
    last_sign <- sign(df.row$last_diff_q5)
    
    duration_signif <- sign(df.row$duration_diff_q5) * sign(df.row$duration_diff_q95)
    duration_sign <- sign(df.row$duration_diff_q5)
    
    auc_signif <- sign(df.row$auc_diff_q5) * sign(df.row$auc_diff_q95)
    auc_sign <- sign(df.row$auc_diff_q5)
    
    
    
    df.add <- data.frame(tissue_idx = rep(tissue.ii, 7) ,
                         route_idx = rep(route.ii, 7) ,
                         metric = c("Probability of positivity (%)", "Time to detectability (days)",
                                    "Time to peak titer (days)", "Peak titer (log10 pfu)",
                                    "Time to undetectability (days)", "Duration of infection (days)",
                                    "AUC"),
                         quantiles = c(paste0(df.row$percent_diff_q5, ", ", df.row$percent_diff_q95),
                                       paste0(df.row$first_diff_q5, ", ", df.row$first_diff_q95),
                                       paste0(df.row$peak_diff_q5, ", ", df.row$peak_diff_q95),
                                       paste0(df.row$titer_diff_q5, ", ", df.row$titer_diff_q95),
                                       paste0(df.row$last_diff_q5, ", ", df.row$last_diff_q95),
                                       paste0(df.row$duration_diff_q5, ", ", df.row$duration_diff_q95),
                                       paste0(df.row$auc_diff_q5, ", ", df.row$auc_diff_q95)),
                         signif = c(percent_signif, first_signif, peak_signif,
                                    titer_signif, last_signif, duration_signif,
                                    auc_signif),
                         sign = c(percent_sign, first_sign, peak_sign,
                                  titer_sign, last_sign, duration_sign,
                                  auc_sign),
                         size = c(as.numeric(df.row$percent_diff_median),
                                  as.numeric(df.row$first_diff_median),
                                  as.numeric(df.row$peak_diff_median),
                                  as.numeric(df.row$titer_diff_median),
                                  as.numeric(df.row$last_diff_median),
                                  as.numeric(df.row$duration_diff_median),
                                  as.numeric(df.row$auc_diff_median)))
    
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

for (metric.ii in unique(df.sigs$metric)){
  for (tissue.ii in unique(df.sigs$tissue_idx)) {
    df.sigs.sub <- subset(df.sigs, metric == metric.ii &
                            tissue_idx == tissue.ii)
    # Max effect should only consider the statistically significant effects
    max_effect <- max(abs(df.sigs.sub$size[df.sigs.sub$symbol != ""]))
    
    df.sigs$size_relative[df.sigs$metric == metric.ii &
                            df.sigs$tissue_idx == tissue.ii] <- 
      df.sigs$size[df.sigs$metric == metric.ii &
                     df.sigs$tissue_idx == tissue.ii] / max_effect
  }
}

# Note the warning that it retunrs -Inf does not affect plotting or computation
#   because these sets don't have a significant effect & will plot as white

# Make the size differences more obvious
df.sigs$size_adj <- df.sigs$size
df.sigs$size_adj[df.sigs$metric != "Probability of positivity (%)"] <- 
  round(df.sigs$size_adj[df.sigs$metric != "Probability of positivity (%)"], 1)
df.sigs$size_adj[df.sigs$metric == "Probability of positivity (%)"] <- 
  round(df.sigs$size_adj[df.sigs$metric == "Probability of positivity (%)"])

#df.sigs$size_adj <- gsub("-", "↓ ", df.sigs$size_adj)
#df.sigs$size_adj <- ifelse(grepl("↓ ", df.sigs$size_adj), df.sigs$size_adj, paste0("↑ ", df.sigs$size_adj))


df.sigs$size_adj <- ifelse(grepl("-", df.sigs$size_adj), df.sigs$size_adj, paste0("+", df.sigs$size_adj))
#df.sigs$size_adj <- gsub("-(?=\\d)", "- ", df.sigs$size_adj, perl = TRUE)


### Plot -------------------------------------------------------------------------

fig.dose.sigs <- ggplot() + 
  geom_tile(data = subset(df.sigs, symbol == ""), 
            aes(y = route_name, x = metric),
            fill = "white", color = "black", linewidth = 0.2) +
  geom_tile(data = subset(df.sigs, symbol != ""), 
            aes(y = route_name, x = metric, 
                fill = tissue_name, alpha = abs(size_relative)),
            color = "black", linewidth = 0.2) +
  geom_text(data = subset(df.sigs, symbol != ""),
            aes(y = route_name, x = metric, label = size_adj), 
            size = 3.1) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 4.5, ymax = 5.5,
           color = "black", fill = NA, linewidth = 1) +
  annotate("rect", xmin = 0.5, xmax = 6.5, ymin = 1.5, ymax = 3.5,
           color = "black", fill = NA, linewidth = 1) +
  scale_fill_manual(values = c("Nose" = df.color$Nose,
                               "Lung" = df.color$Lung,
                               "Lower GI" = df.color$Lower.GI)) +
  facet_grid(factor(tissue_name, levels = c("Nose", "Lung", "Lower GI")) ~ "Diff. between Max & Min Dose") +
  labs(fill = "") +
  scale_y_discrete(expand = c(0, 0),
                   limits = rev(c("IN", "IT", "IN+IT", "AE", "IG"))) +
  scale_x_discrete(expand = c(0, 0),
                   limits = c("Probability of positivity (%)",
                              "Time to detectability (days)",
                              "Time to peak titer (days)",
                              "Peak titer (log10 pfu)",
                              "Time to undetectability (days)",
                              "Duration of infection (days)")) +
  coord_cartesian(clip = "off") +
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
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5)); fig.dose.sigs



## Panel D: Trajectory, AE ------------------------------------------------------

# Subset to AE data with culture
df.dose.fixed.traj.ae <- subset(df.dose.fixed, assay_idx == 4 & route_idx == 4)

# Get median & quantiles for each dose, tissue, and assay
df.dose.fixed.quantiles.ae <- df.dose.fixed.traj.ae %>%
  group_by(dose_total, tissue_idx, route_idx, assay_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
  ), 
  list(
    median = ~median(.),
    q5 = ~quantile(., probs = 0.05),
    q95 = ~quantile(., probs = 0.95)))) 

# Plot for the nose
fig.traj.nose <- 
  ggplot(subset(df.dose.fixed.quantiles.ae, assay_idx == 4 & tissue_idx == 1)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(.~ "Aerosol exposure", ncol = 1) +
  labs(x = "days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#D8CDDF", high = df.color$Nose) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(2, 2, 0, 2),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                         linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5)); fig.traj.nose


# Plot for the lung
fig.traj.lung <- 
  ggplot(subset(df.dose.fixed.quantiles.ae, assay_idx == 4 & tissue_idx == 4)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FBEFC6", high = df.color$Lung) +
  theme(text = element_text(size = 11),
        legend.position =  "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 2, 0, 2),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        strip.text = element_blank()); fig.traj.lung


# Plot for the GI
fig.traj.gi <- 
  ggplot(subset(df.dose.fixed.quantiles.ae, assay_idx == 4 & tissue_idx == 6)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FFDAD6", high = df.color$Lower.GI) +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title.y = element_blank(),
        axis.title = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        strip.text = element_blank()); fig.traj.gi


fig.traj.legend <- 
  ggplot(subset(df.dose.fixed.quantiles, assay_idx == 4 & tissue_idx == 6)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14)) +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "Days post infection", y = "Viral titer (log10)", 
       color = "Dose\n(log10 pfu)") +
  scale_color_gradient(low = "grey88", high = "grey22") +
  theme(text = element_text(size = 9),
        legend.position = c(0.3, 0.8),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        strip.text = element_blank()); fig.traj.legend

# Extract the legend
g <- ggplotGrob(fig.traj.legend)
legend <- g$grobs[[which(g$layout$name == "guide-box-inside")]] # legend at the bottom
fig.legend <- ggdraw(legend); fig.legend # wrap as ggplot compatible item


# Add legend to nose trajectory panel
fig.traj.nose.legend <- 
  fig.traj.nose +
  annotation_custom(
    grob = ggplotGrob(fig.legend),
    xmin = 10.2, xmax = 11.5,
    ymin = -1.5, ymax = 4.2
  ); fig.traj.nose.legend


fig.traj.ae <- fig.traj.nose.legend + labs(tag = "D") + 
  fig.traj.lung + fig.traj.gi + plot_layout(ncol = 1); fig.traj.ae


## Panel E: Trajectory, IT ------------------------------------------------------------

df.dose.fixed.it.traj <- subset(df.dose.fixed, assay_idx == 4 & route_idx == 2)


df.dose.fixed.quantiles.it.traj <- df.dose.fixed.it.traj %>%
  group_by(dose_total, tissue_idx, route_idx, assay_idx) %>%
  summarise(across(c(percent_positive, first_pos_median, 
                     peak_median, titer_mean, last_median,
  ), 
  list(
    median = ~median(.),
    q5 = ~quantile(., probs = 0.05),
    q95 = ~quantile(., probs = 0.95)))) 



fig.traj.nose.it <- 
  ggplot(subset(df.dose.fixed.quantiles.it.traj, assay_idx == 4 & tissue_idx == 1)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(.~ "Intratracheal exposure", ncol = 1) +
  labs(x = "days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#D8CDDF", high = df.color$Nose) +
  #guides(color = "none") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(2, 2, 0, 2),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        linetype = "solid"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_text()); fig.traj.nose.it


fig.traj.lung.it <- 
  ggplot(subset(df.dose.fixed.quantiles.it.traj, assay_idx == 4 & tissue_idx == 4)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FBEFC6", high = df.color$Lung) +
  #guides(color = "none") +
  theme(text = element_text(size = 11),
        legend.position =  "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(0, 2, 0, 2),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_blank()); fig.traj.lung.it



fig.traj.gi.it <- 
  ggplot(subset(df.dose.fixed.quantiles.it.traj, assay_idx == 4 & tissue_idx == 6)) +
  geom_segment(aes(x = first_pos_median_median,
                   xend = peak_median_median,
                   y = 0, yend = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  geom_segment(aes(x = peak_median_median,
                   xend = last_median_median,
                   yend = 0, y = titer_mean_median,
                   color = dose_total),
               linewidth = 0.5, alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14), clip = "on") +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "Days post infection", y = "Viral titer (log10 pfu)",linetype = "Dose (pfu)") +
  scale_color_gradient(low = "#FFDAD6", high = df.color$Lower.GI) +
  #guides(color = "none") +
  theme(text = element_text(size = 11),
        legend.position = "none",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        axis.title = element_text(size = 9),
        axis.title.y = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.5, linetype = "solid"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.text = element_blank()); fig.traj.gi.it


fig.traj.it <- fig.traj.nose.it + labs(tag = "E") + 
  fig.traj.lung.it + fig.traj.gi.it + plot_layout(ncol = 1); fig.traj.it


## Legend ----------------------------------------------------------------

fig.color.legend <- 
  ggplot(subset(df.dose.fixed.quantiles, assay_idx == 4)) +
  geom_point(aes(x = first_pos_median_median,
                 y = 0, 
                 color = as.character(tissue_idx)),
             alpha = 1) +
  scale_x_continuous(breaks = seq(0, 100, 4)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(0, 14)) +
  facet_wrap(tissue_idx ~., ncol = 1, scales = "free_y") +
  labs(x = "Days post infection", y = "Viral titer (log10)", 
       color = "") +
  scale_color_manual(values = c("1" = "#85649B",
                                "4" = "#E3B710",
                                "6" = "#F11B00"),
                     labels = c("Nose", "Throat", "Lower GI")) +
  theme(text = element_text(size = 9),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.5, "line"),
        legend.key = element_rect(fill = "white", color = "white"),
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(0, 2, 2, 2),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 1, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.5),
        strip.text = element_blank()); fig.color.legend


# Extract the legend
g <- ggplotGrob(fig.color.legend)
legend <- g$grobs[[which(g$layout$name == "guide-box-top")]] # legend at the bottom
fig.legend.color <- ggdraw(legend); fig.legend.color # wrap as ggplot compatible item



# Combine  ------------------------------------------------------------

fig.dose.comb <- 
  ((fig.id50 + labs(tag = "A")) / (fig.dose.percent.ci +labs(tag = "B"))) |
  (fig.dose.sigs + labs(tag = "C")) |
  fig.traj.ae | fig.traj.it; fig.dose.comb

fig.dose.comb <- (fig.legend.color / fig.dose.comb) + plot_layout(heights = c(0.1, 1))


# Save  ------------------------------------------------------------------------

ggsave("./outputs/figures/fig6-dose-effects-culture.pdf",
       plot = fig.dose.comb,
       width = 8.8, 
       height = 5.8,
       dpi = 600)
