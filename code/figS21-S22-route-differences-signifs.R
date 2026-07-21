

# Generate pairwise differences ------------------------------------------------

# Load predictions
df.route.diffs <- readRDS("./outputs/predictions/pred-route-differences-1000.RDS")

# Restrict to mid-level dose & culture assays
df.diffs <- subset(df.route.diffs, assay_idx == 4 & dose_total == 4)
df.diffs.totalRNA <- subset(df.route.diffs, assay_idx == 1 & dose_total == 4)

# Assign tissue names
df.diffs <- assign_tissue_names(df.diffs)
df.diffs.totalRNA <- assign_tissue_names(df.diffs.totalRNA)

# Set upper GI differences to zero for all metrics but probability & detectability
df.diffs$diff_peak[df.diffs$tissue_idx == 5] <- 0
df.diffs$diff_titer[df.diffs$tissue_idx == 5] <- 0
df.diffs$diff_last[df.diffs$tissue_idx == 5] <- 0
df.diffs$diff_duration[df.diffs$tissue_idx == 5] <- 0

df.diffs.totRNA$diff_peak[df.diffs.totRNA$tissue_idx == 5] <- 0
df.diffs.totRNA$diff_titer[df.diffs.totRNA$tissue_idx == 5] <- 0
df.diffs.totRNA$diff_last[df.diffs.totRNA$tissue_idx == 5] <- 0
df.diffs.totRNA$diff_duration[df.diffs.totRNA$tissue_idx == 5] <- 0


# Check for significance for culture -------------------------------------------

# Get quantiles for differences
df.diffs.median <- df.diffs %>%
  group_by(tissue_name, cofactor1, cofactor2, assay_idx) %>%
  summarise(across(c(diff_percent, diff_first, 
                     diff_peak, diff_titer, diff_last,
                     diff_duration, diff_auc), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     q5 = ~quantile(., probs = 0.05, na.rm = TRUE),
                     q95 = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


# Create dataframe to store significances
df.sigs <- data.frame(tissue_name = numeric(),
                      route1 = numeric(),
                      route2 = numeric(),
                      metric = character(),
                      signif = numeric(),
                      sign = numeric())

# Loop over each metric, tissue, and route comparison to check significance
for (tissue.ii in unique(df.diffs.median$tissue_name)) {
  
  for (route1 in 1:5) {
    
    comp_routes <- unique(df.diffs.median$cofactor2[df.diffs.median$cofactor1 == route1])
    
    for (route2 in comp_routes) {
      
      df.row <- subset(df.diffs.median, cofactor1 == route1 & 
                         cofactor2 == route2 &
                         tissue_name == tissue.ii)
      
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
      
      
      df.add <- data.frame(tissue_name = rep(tissue.ii, 7) ,
                           route1 = rep(route1, 7) ,
                           route2 = rep(route2, 7),
                           metric = c("Probability of positivity", "Time to detectability",
                                      "Time to peak titer", "Peak titer",
                                      "Time to undetectability", "Duration of infection",
                                      "AUC"),
                           signif = c(percent_signif, first_signif, peak_signif,
                                      titer_signif, last_signif, duration_signif,
                                      auc_signif),
                           sign = c(percent_sign, first_sign, peak_sign,
                                    titer_sign, last_sign, duration_sign,
                                    auc_sign),
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

# Update significanc categoricals
df.sigs$signif[df.sigs$signif == -1] <- 0
df.sigs$signif[df.sigs$signif == 1 & df.sigs$sign == -1] <- -1

# Add symbol for effect direction
df.sigs$symbol[df.sigs$signif == "1"] <- "+"
df.sigs$symbol[df.sigs$signif == "-1"] <- "-"
df.sigs$symbol[df.sigs$signif == 0] <- ""

# Add route names
df.sigs$route1_name[df.sigs$route1 == 1] <- "IN"
df.sigs$route1_name[df.sigs$route1 == 2] <- "IT"
df.sigs$route1_name[df.sigs$route1 == 3] <- "IN+IT"
df.sigs$route1_name[df.sigs$route1 == 4] <- "AE"
df.sigs$route1_name[df.sigs$route1 == 5] <- "IG"
df.sigs$route2_name[df.sigs$route2 == 1] <- "IN"
df.sigs$route2_name[df.sigs$route2 == 2] <- "IT"
df.sigs$route2_name[df.sigs$route2 == 3] <- "IN+IT"
df.sigs$route2_name[df.sigs$route2 == 4] <- "AE"
df.sigs$route2_name[df.sigs$route2 == 5] <- "IG"

# Specify route comparison
df.sigs$route_comp <- paste0(df.sigs$route1_name, 
                             " vs. ", 
                             df.sigs$route2_name)

# Factor metrics for plotting order
df.sigs$metric <- factor(df.sigs$metric,
                         levels = c("Probability of positivity",
                                    "Time to detectability",
                                    "Time to peak titer",
                                    "Peak titer",
                                    "Time to undetectability",
                                    "Duration of infection",
                                    "AUC"))

# Factor tissue names for plotting order
df.sigs$tissue_name <- factor(df.sigs$tissue_name,
                              levels = c("Nose", "Throat",
                                         "Trachea", "Lung",
                                         "Upper GI", "Lower GI"))



# Check for significance for total RNA -----------------------------------------

# Get quantiles for differences
df.diffs.totalRNA.median <- df.diffs.totalRNA %>%
  group_by(tissue_name, cofactor1, cofactor2, assay_idx) %>%
  summarise(across(c(diff_percent, diff_first, 
                     diff_peak, diff_titer, diff_last,
                     diff_duration, diff_auc), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     q5 = ~quantile(., probs = 0.05, na.rm = TRUE),
                     q95 = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


# Create dataframe to store significances
df.sigs.totalRNA <- data.frame(tissue_name = numeric(),
                               route1 = numeric(),
                               route2 = numeric(),
                               metric = character(),
                               signif = numeric(),
                               sign = numeric())

# Loop over each metric, tissue, and route comparison to check significance
for (tissue.ii in unique(df.diffs.totalRNA.median$tissue_name)) {
  
  for (route1 in 1:5) {
    
    comp_routes <- unique(df.diffs.totalRNA.median$cofactor2[df.diffs.totalRNA.median$cofactor1 == route1])
    
    for (route2 in comp_routes) {
      
      df.row <- subset(df.diffs.totalRNA.median, cofactor1 == route1 & 
                                                 cofactor2 == route2 &
                                                 tissue_name == tissue.ii)
      
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
      
      
      df.add <- data.frame(tissue_name = rep(tissue.ii, 7) ,
                           route1 = rep(route1, 7) ,
                           route2 = rep(route2, 7),
                           metric = c("Probability of positivity", "Time to detectability",
                                      "Time to peak titer", "Peak titer",
                                      "Time to undetectability", "Duration of infection",
                                      "AUC"),
                           signif = c(percent_signif, first_signif, peak_signif,
                                      titer_signif, last_signif, duration_signif,
                                      auc_signif),
                           sign = c(percent_sign, first_sign, peak_sign,
                                    titer_sign, last_sign, duration_sign,
                                    auc_sign),
                           size = c(as.numeric(df.row$diff_percent_median),
                                    as.numeric(df.row$diff_first_median),
                                    as.numeric(df.row$diff_peak_median),
                                    as.numeric(df.row$diff_titer_median),
                                    as.numeric(df.row$diff_last_median),
                                    as.numeric(df.row$diff_duration_median),
                                    as.numeric(df.row$diff_auc_median)))
      
      df.sigs.totalRNA <- rbind(df.sigs.totalRNA, df.add)
      
    }
  }
}

# Update significanc categoricals
df.sigs.totalRNA$signif[df.sigs.totalRNA$signif == -1] <- 0
df.sigs.totalRNA$signif[df.sigs.totalRNA$signif == 1 & df.sigs.totalRNA$sign == -1] <- -1

# Add symbol for effect direction
df.sigs.totalRNA$symbol[df.sigs.totalRNA$signif == "1"] <- "+"
df.sigs.totalRNA$symbol[df.sigs.totalRNA$signif == "-1"] <- "-"
df.sigs.totalRNA$symbol[df.sigs.totalRNA$signif == 0] <- ""

# Add route names
df.sigs.totalRNA$route1_name[df.sigs.totalRNA$route1 == 1] <- "IN"
df.sigs.totalRNA$route1_name[df.sigs.totalRNA$route1 == 2] <- "IT"
df.sigs.totalRNA$route1_name[df.sigs.totalRNA$route1 == 3] <- "IN+IT"
df.sigs.totalRNA$route1_name[df.sigs.totalRNA$route1 == 4] <- "AE"
df.sigs.totalRNA$route1_name[df.sigs.totalRNA$route1 == 5] <- "IG"
df.sigs.totalRNA$route2_name[df.sigs.totalRNA$route2 == 1] <- "IN"
df.sigs.totalRNA$route2_name[df.sigs.totalRNA$route2 == 2] <- "IT"
df.sigs.totalRNA$route2_name[df.sigs.totalRNA$route2 == 3] <- "IN+IT"
df.sigs.totalRNA$route2_name[df.sigs.totalRNA$route2 == 4] <- "AE"
df.sigs.totalRNA$route2_name[df.sigs.totalRNA$route2 == 5] <- "IG"

# Specify route comparison
df.sigs.totalRNA$route_comp <- paste0(df.sigs.totalRNA$route1_name, 
                             " vs. ", 
                             df.sigs.totalRNA$route2_name)

# Factor metrics for plotting order
df.sigs.totalRNA$metric <- factor(df.sigs.totalRNA$metric,
                         levels = c("Probability of positivity",
                                    "Time to detectability",
                                    "Time to peak titer",
                                    "Peak titer",
                                    "Time to undetectability",
                                    "Duration of infection",
                                    "AUC"))

# Factor tissue names for plotting order
df.sigs.totalRNA$tissue_name <- factor(df.sigs.totalRNA$tissue_name,
                              levels = c("Nose", "Throat",
                                         "Trachea", "Lung",
                                         "Upper GI", "Lower GI"))


# Figure 21 --------------------------------------------------------------------

## Top Panels: Culture --------------------------------------------------------------------

df.sigs$route_comp_factor <- factor(df.sigs$route_comp,
                                    levels = rev(c("IN vs. IT", "IN vs. IN+IT",
                                                   "IN vs. AE", "IN vs. IG",
                                                   "IT vs. IN+IT", "IT vs. AE",
                                                   "IT vs. IG",
                                                   "IN+IT vs. AE", "IN+IT vs. IG",
                                                   "AE vs. IG")))

df.sigs.upgi <- subset(df.sigs, tissue_name == "Upper GI" & metric %notin% c("Probability of positivity","Time to detectability"))
df.sigs <- subset(df.sigs, !(tissue_name == "Upper GI" & metric %notin% c("Probability of positivity",
                                                                          "Time to detectability")))


### Plot ------------------------------------------------------------------------

fig.culture.sigs <- ggplot(df.sigs) + 
  geom_tile(aes(x = tissue_name, y = route_comp),
            fill = "white", color = "black", linewidth = 0.5) +
  geom_tile(data = df.sigs.upgi,
              aes(x = tissue_name, y = route_comp),
            fill = "grey55", color = "black", linewidth = 0.5) +
  geom_tile(data = subset(df.sigs, symbol != ""),
            aes(x = tissue_name, y = route_comp, fill = symbol),
            color = "black", linewidth = 0.5) +
  geom_text(aes(x = tissue_name, y = route_comp_factor, label = symbol),
            size = 4, fontface = "bold") +
  facet_wrap(.~ metric, ncol = 4) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("-" = "#457b9d",
                               "+" = "#e63946"),
                    labels = c("X smaller than Y", 
                               "X larger than Y")) +
  labs(fill = "Difference (X vs. Y)") +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.culture.sigs



## Bottom Panels: Total RNA --------------------------------------------------------------------

df.sigs.totalRNA$route_comp_factor <- factor(df.sigs.totalRNA$route_comp,
                                    levels = rev(c("IN vs. IT", "IN vs. IN+IT",
                                                   "IN vs. AE", "IN vs. IG",
                                                   "IT vs. IN+IT", "IT vs. AE",
                                                   "IT vs. IG",
                                                   "IN+IT vs. AE", "IN+IT vs. IG",
                                                   "AE vs. IG")))

df.sigs.totalRNA.upgi <- subset(df.sigs.totalRNA, tissue_name == "Upper GI" & metric %notin% c("Probability of positivity","Time to detectability"))
df.sigs.totalRNA <- subset(df.sigs.totalRNA, !(tissue_name == "Upper GI" & metric %notin% c("Probability of positivity",
                                                                          "Time to detectability")))


### Plot ------------------------------------------------------------------------

fig.sigs.totalRNA <- ggplot(df.sigs.totalRNA) + 
  geom_tile(aes(x = tissue_name, y = route_comp),
            fill = "white", color = "black", linewidth = 0.5) +
  geom_tile(data = df.sigs.totalRNA.upgi,
            aes(x = tissue_name, y = route_comp),
            fill = "grey55", color = "black", linewidth = 0.5) +
  geom_tile(data = subset(df.sigs.totalRNA, symbol != ""),
            aes(x = tissue_name, y = route_comp, fill = symbol),
            color = "black", linewidth = 0.5) +
  geom_text(aes(x = tissue_name, y = route_comp_factor, label = symbol),
            size = 4, fontface = "bold") +
  facet_wrap(.~ metric, ncol = 4) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c("-" = "#457b9d",
                               "+" = "#e63946"),
                    labels = c("X smaller than Y", 
                               "X larger than Y")) +
  labs(fill = "Difference (X vs. Y)") +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "line"),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(fill = "transparent", color = "transparent"),
        legend.box.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.sigs.totalRNA


## Combine ---------------------------------------------------------------------

fig.comb.sigs <- 
  (fig.culture.sigs + labs(tag = "Culture Results") + theme(legend.position = "none")) / 
  fig.sigs.totalRNA + labs(tag = "Total RNA Results"); fig.comb.sigs


## Save ------------------------------------------------------------------------

ggsave('./outputs/figures/figS21-route-differences-signifs.png',
       plot = fig.comb.sigs,
       width = 9.5, 
       height = 10,
       dpi = 600)



# Figure 22 --------------------------------------------------------------------

## Panel A: Culture Comparison Counts ------------------------------------------

# Compute total number of differences per exposure route comparison
df.sigs.counts <- df.sigs %>%
  group_by(route_comp_factor, route1_name, route2_name) %>%
  summarize(num = sum(abs(signif)))

fig.culture.counts <- ggplot(df.sigs.counts) + 
  geom_col(aes(x = reorder(route_comp_factor, -num), y = num)) +
  geom_text(aes(x = reorder(route_comp_factor, -num), y = num - 2,
                label = num), color = "white") +
  scale_x_discrete() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 21), breaks = seq(0, 20, 5)) +
  labs (y = "Total number of\nsignificant differences\nfor culture") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x  = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.culture.counts


## Panel B: Culture Counts Per Route -------------------------------------------

df.sigs.counts.byroute <- df.sigs.counts %>%
  pivot_longer(cols = c(route1_name, route2_name),
               names_to = "route_position",
               values_to = "route") %>%
  group_by(route) %>%
  summarise(total_num = sum(num, na.rm = TRUE))

fig.culture.counts.byroute <- ggplot(df.sigs.counts.byroute) + 
  geom_col(aes(x = reorder(route, -total_num), y = total_num)) +
  geom_text(aes(x = reorder(route, -total_num), y = total_num - 5,
                label = total_num), color = "white") +
  scale_x_discrete() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) +
  labs (y = "Total number of\nsignificant differences\nwith all other routes\nfor culture") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x  = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.culture.counts.byroute

## Panel C: Total RNA Comparison Counts ----------------------------------------

# Compute total number of differences per exposure route comparison
df.sigs.totalRNA.counts <- df.sigs.totalRNA %>%
  group_by(route_comp_factor, route1_name, route2_name) %>%
  summarize(num = sum(abs(signif)))

fig.totalRNA.counts <- ggplot(df.sigs.totalRNA.counts) + 
  geom_col(aes(x = reorder(route_comp_factor, -num), y = num)) +
  geom_text(aes(x = reorder(route_comp_factor, -num), y = num - 2,
                label = num), color = "white") +
  scale_x_discrete() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22), breaks = seq(0, 20, 5)) +
  labs (y = "Total number of\nsignificant differences\nfor total RNA") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x  = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.totalRNA.counts


## Panel D: Total RNA Counts Per Route -------------------------------------------

df.sigs.totalRNA.counts.byroute <- df.sigs.totalRNA.counts %>%
  pivot_longer(cols = c(route1_name, route2_name),
               names_to = "route_position",
               values_to = "route") %>%
  group_by(route) %>%
  summarise(total_num = sum(num, na.rm = TRUE))

fig.totalRNA.counts.byroute <- ggplot(df.sigs.totalRNA.counts.byroute) + 
  geom_col(aes(x = reorder(route, -total_num), y = total_num)) +
  geom_text(aes(x = reorder(route, -total_num), y = total_num - 5,
                label = total_num), color = "white") +
  scale_x_discrete() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 70)) +
  labs (y = "Total number of\nsignificant differences\nwith all other routes\nfor total RNA") +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x  = element_blank(),
        plot.background = element_rect(fill = "white", colour = NA_character_),
        panel.background = element_rect(fill = "white",
                                        colour = "black", linewidth = 0.5,
                                        size = 0.8, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black",
                                        linewidth = 0.8)); fig.totalRNA.counts.byroute



## Combine -------------------------------------------------------------------

fig.counts <- fig.culture.counts + labs(tag = "a") + 
  fig.culture.counts.byroute + labs(tag = "b") +
  fig.totalRNA.counts + labs(tag = "c") + 
  fig.totalRNA.counts.byroute + labs(tag = "d") +
  plot_layout(nrow = 2, widths = c(10, 5), byrow =  TRUE); fig.counts


## Plot -----------------------------------------------------------------------

ggsave('./outputs/figures/figS22-route-differences-signifs-counts.png',
       plot = fig.counts,
       width = 6, 
       height = 5,
       dpi = 600)

