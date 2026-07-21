# This file: - Computes the significance of the effects of dose on all metrics, based on total RNA


# Prep -------------------------------------------------------------------------

df.dose.full <- fread("./outputs/predictions/pred-across-doses-450.csv")

# Subset to only extreme ends of the dose range
df.dose.full.sig <- subset(df.dose.full, dose_total %in% c(1.2, 7.4) & assay_idx == 1)

# Calculate differences for each metric (note this takes a long time to run)
df.dose.full.sig.diff <- df.dose.full.sig
for (route.ii in unique(df.dose.full.sig.diff$route_idx)) {
  for (sample.ii in unique(df.dose.full.sig.diff$sample_num)) {
    for (tissue.ii in unique(df.dose.full.sig.diff$tissue_idx)) {
      for (sex.ii in unique(df.dose.full.sig.diff$sex_idx)) {
        for (age.ii in unique(df.dose.full.sig.diff$age_idx)) {
          for (sp.ii in unique(df.dose.full.sig.diff$sp_idx)) {
            
            df.sub <- subset(df.dose.full.sig.diff, sample_num == sample.ii &
                               tissue_idx == tissue.ii & route_idx == route.ii & 
                               sex_idx == sex.ii & age_idx == age.ii & sp_idx == sp.ii)
            
            percent_diff <- (df.sub$percent_positive[df.sub$dose_total == 7.4] * 100 -#* 10 -
                               df.sub$percent_positive[df.sub$dose_total == 1.2] * 100) #* 10)
            first_diff <- (df.sub$first_pos_median[df.sub$dose_total == 7.4] -
                             df.sub$first_pos_median[df.sub$dose_total == 1.2]) 
            peak_diff <- (df.sub$peak_median[df.sub$dose_total == 7.4] -
                            df.sub$peak_median[df.sub$dose_total == 1.2]) 
            titer_diff <- (df.sub$titer_mean[df.sub$dose_total == 7.4] -
                             df.sub$titer_mean[df.sub$dose_total == 1.2]) 
            last_diff <- (df.sub$last_median[df.sub$dose_total == 7.4] -
                            df.sub$last_median[df.sub$dose_total == 1.2]) 
            duration_diff <- (df.sub$duration_median[df.sub$dose_total == 7.4] -
                                df.sub$duration_median[df.sub$dose_total == 1.2]) 
            
            df.dose.full.sig.diff$percent_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                                 df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                                 df.dose.full.sig.diff$route_idx == route.ii &
                                                 df.dose.full.sig.diff$sex_idx == sex.ii &
                                                 df.dose.full.sig.diff$age_idx == age.ii &
                                                 df.dose.full.sig.diff$sp_idx == sp.ii] <- c(percent_diff, NA)
            
            df.dose.full.sig.diff$first_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                               df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                               df.dose.full.sig.diff$route_idx == route.ii &
                                               df.dose.full.sig.diff$sex_idx == sex.ii &
                                               df.dose.full.sig.diff$age_idx == age.ii &
                                               df.dose.full.sig.diff$sp_idx == sp.ii] <- c(first_diff, NA)
            
            df.dose.full.sig.diff$peak_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                              df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                              df.dose.full.sig.diff$route_idx == route.ii &
                                              df.dose.full.sig.diff$sex_idx == sex.ii &
                                              df.dose.full.sig.diff$age_idx == age.ii &
                                              df.dose.full.sig.diff$sp_idx == sp.ii] <- c(peak_diff, NA)
            
            df.dose.full.sig.diff$titer_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                               df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                               df.dose.full.sig.diff$route_idx == route.ii &
                                               df.dose.full.sig.diff$sex_idx == sex.ii &
                                               df.dose.full.sig.diff$age_idx == age.ii &
                                               df.dose.full.sig.diff$sp_idx == sp.ii] <- c(titer_diff, NA)
            
            df.dose.full.sig.diff$last_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                              df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                              df.dose.full.sig.diff$route_idx == route.ii &
                                              df.dose.full.sig.diff$sex_idx == sex.ii &
                                              df.dose.full.sig.diff$age_idx == age.ii &
                                              df.dose.full.sig.diff$sp_idx == sp.ii] <- c(last_diff, NA)
            
            df.dose.full.sig.diff$duration_diff[df.dose.full.sig.diff$sample_num == sample.ii &
                                                  df.dose.full.sig.diff$tissue_idx == tissue.ii &
                                                  df.dose.full.sig.diff$route_idx == route.ii &
                                                  df.dose.full.sig.diff$sex_idx == sex.ii &
                                                  df.dose.full.sig.diff$age_idx == age.ii &
                                                  df.dose.full.sig.diff$sp_idx == sp.ii] <- c(duration_diff, NA)
            
            #print(df.dose.full.sig.diff)
            
          }
        }
      }
      
    }
  }
}

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
                     duration_diff), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     q5 = ~quantile(., probs = 0.05, na.rm = TRUE),
                     q95 = ~quantile(., probs = 0.95, na.rm = TRUE))))


# Determine whether the difference is signficant (based on 90% CI)
df.sigs <- data.frame(tissue_idx = numeric(),
                      route_idx = numeric(),
                      metric = character(),
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
    
    
    df.add <- data.frame(tissue_idx = rep(tissue.ii, 6) ,
                         route_idx = rep(route.ii, 6) ,
                         metric = c("Probability of positivity (%)", "Time to detectability (days)",
                                    "Time to peak titer (days)", "Peak titer (log10 pfu)",
                                    "Time to undetectability (days)", "Duration of infection (days)"),
                         signif = c(percent_signif, first_signif, peak_signif,
                                    titer_signif, last_signif, duration_signif),
                         sign = c(percent_sign, first_sign, peak_sign,
                                  titer_sign, last_sign, duration_sign),
                         size = c(as.numeric(df.row$percent_diff_median),
                                  as.numeric(df.row$first_diff_median),
                                  as.numeric(df.row$peak_diff_median),
                                  as.numeric(df.row$titer_diff_median),
                                  as.numeric(df.row$last_diff_median),
                                  as.numeric(df.row$duration_diff_median)))
    
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


# Plot -------------------------------------------------------------------------

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
  facet_grid(factor(tissue_name, levels = c("Nose", "Lung", "Lower GI")) ~ "Dose Effect") +
  labs(fill = "Dose Effect") +
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
        strip.text.x = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 1)); fig.dose.sigs


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS30-dose-signifs-totalRNA.png",
       plot = fig.dose.sigs,
       width = 2.8, 
       height = 5.8,
       dpi = 600)

