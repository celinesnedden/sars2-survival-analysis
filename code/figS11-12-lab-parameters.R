# This file: - Plots the median & 90% CrI for all of the lab parameters

# Prep -------------------------------------------------------------------------

library("ggh4x")

# Load the model fit
fit <- readRDS('./outputs/fits/fit-main.RDS')

# Extract the draws from all parameters
df.draws <- fit$draws(format = "df")

# Get only the lab effect parameters
df.draws <- df.draws[, grepl("lab", names(df.draws))]

# Extract the medians & quantiles
df.summary <- summarise_draws(df.draws, 
                              median, 
                              ~quantile(.x, probs = c(0.05, 0.95)))

# Round the median & quantiles to 2 digits
df.summary$median_rounded <- round(df.summary$median, digits = 2)
df.summary$q5_rounded <- round(df.summary$`5%`, digits = 2)
df.summary$q95_rounded <- round(df.summary$`95%`, digits = 2)
df.summary$quantiles_combined <- paste0(df.summary$median_rounded,
                                        " [", df.summary$q5_rounded,
                                        ", ", df.summary$q95_rounded,
                                        "]")

# Remove text from variable names to just have the metric
df.summary$variable <- str_remove_all(df.summary$variable, "median")
df.summary$variable <- str_remove_all(df.summary$variable, "_")
df.summary$variable <- str_remove_all(df.summary$variable, "lab")

# Get the organ group 
df.summary$location_id <- sub(".*,\\s*([0-9])[^0-9]?.*", "\\1", df.summary$variable)
df.summary$lab_id <- sub(".*\\[([0-9]+),.*", "\\1", df.summary$variable)
df.summary$metric <- sub("\\[.*", "", df.summary$variable)

# Load the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS")

# Get the unique assay-lab pairs
df.pairs <- data.frame(lab = dat.stan$lab,
                       organ_group = dat.stan$organ_location,
                       assay = dat.stan$assay) %>%
  distinct(lab, assay, organ_group)

labs_with_unknown <- unique(df.pairs$lab[df.pairs$assay == -9999])
df.pairs$assay[df.pairs$assay == -9999] <- 1

# Get the labs that correspond with each assay
labs_assay1 <- df.pairs$lab[df.pairs$assay == 1]
labs_assay2 <- df.pairs$lab[df.pairs$assay == 2]
labs_assay3 <- df.pairs$lab[df.pairs$assay == 3]
labs_assay4 <- df.pairs$lab[df.pairs$assay == 4]

# Adding the assay type to each lab
df.summary$assay_idx[df.summary$lab_id %in% labs_assay1] <- 1
df.summary$assay_idx[df.summary$lab_id %in% labs_assay2] <- 2
df.summary$assay_idx[df.summary$lab_id %in% labs_assay3] <- 3
df.summary$assay_idx[df.summary$lab_id %in% labs_assay4] <- 4

# Add assay names
df.summary$assay_name[df.summary$lab_id %in% labs_assay1] <- "Total RNA"
df.summary$assay_name[df.summary$lab_id %in% labs_assay2] <- "gRNA"
df.summary$assay_name[df.summary$lab_id %in% labs_assay3] <- "sgRNA"
df.summary$assay_name[df.summary$lab_id %in% labs_assay4] <- "Culture"
df.summary$assay_name[is.na(df.summary$assay_name)] <- "Unknown"
df.summary$assay_name <- factor(df.summary$assay_name, 
                                levels = c("Total RNA", "gRNA", 
                                           "sgRNA", "Culture", "Unknown"))

# Add tissue group names
df.summary$organ_group[df.summary$location_id == 1] <- "URT"
df.summary$organ_group[df.summary$location_id == 2] <- "LRT"
df.summary$organ_group[df.summary$location_id == 3] <- "GI"
df.summary$organ_group <- factor(df.summary$organ_group, 
                                 levels = c("URT", "LRT", "GI"))

# Add metric names
df.summary$metric[df.summary$metric == "first"] <- "Time to detectability"
df.summary$metric[df.summary$metric == "peak"] <- "Time to peak titer"
df.summary$metric[df.summary$metric == "percent"] <- "Probability of positivity"
df.summary$metric[df.summary$metric == "titer"] <- "Peak titer"
df.summary$metric[df.summary$metric == "last"] <- "Time to undetectability"
df.summary$metric[df.summary$metric == "titerobssd"] <- "Peak titer SD"
df.summary$metric <- factor(df.summary$metric,
                            levels = c("Probability of positivity",
                                       "Time to detectability",
                                       "Time to peak titer",
                                       "Peak titer", "Peak titer SD",
                                       "Time to undetectability"))

# Flag params for tissues without fitting info for each lab
df.summary.data <- df.summary[0, ]
for (row_num in 1:nrow(df.pairs)) {
  df.next <- subset(df.summary, lab_id == df.pairs$lab[row_num] & 
                                assay_idx == df.pairs$assay[row_num] & 
                                location_id == df.pairs$organ_group[row_num])
  
  df.summary.data <- rbind(df.summary.data, df.next)
}
df.summary.extrap <- anti_join(df.summary, df.summary.data)

# Combine them for plotting
df.summary.data$type <- "With Data"
df.summary.extrap$type <- "Without Data"
df.summary.comb <- rbind(df.summary.data, df.summary.extrap)

# Reset labs with unknown assay
df.summary.comb$assay_name[df.summary.comb$lab_id %in% labs_with_unknown] <- "Unknown"


# Fig S11 -------------------------------------------------------------------------

## Plot -------------------------------------------------------------------------

fig <- ggplot(subset(df.summary.comb, metric %in% c("Probability of positivity", 
                                                     "Time to detectability", 
                                                     "Time to peak titer"))) +
  geom_vline(aes(xintercept = 0), color = "black") +
  geom_segment(aes(x = q5_rounded, xend = q95_rounded,
                   y = as.character(lab_id), yend = as.character(lab_id),
                   color = type),
               linewidth = 0.1) +
  geom_point(aes(x = median_rounded, y = as.character(lab_id),
                 fill = type, color = type), 
              color = "black", size = 1, shape = 21) + 
  facet_nested(assay_name ~ metric + organ_group, 
               scales = "free", space = "free_y") +
  scale_y_discrete(breaks = as.character(seq(0, 170, 1))) +
  labs(x = "Parameter estimate", y = "Lab") + 
  theme(text = element_text(size = 9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust = 1),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.spacing = unit(2, "pt"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
        strip.text.x = element_text(face = "bold")
        ); fig

## Save -------------------------------------------------------------------------

ggsave('./outputs/figures/figS11-lab-parameters.png',
       plot = fig,
       width = 7, 
       height = 7,
       dpi = 600)


# Fig S12 ----------------------------------------------------------------------

## Plot ------------------------------------------------------------------------

fig <- ggplot(subset(df.summary.comb, metric %notin% c("Probability of positivity", 
                                                       "Time to detectability", 
                                                       "Time to peak titer"))) +
  geom_vline(aes(xintercept = 0), color = "black") +
  geom_segment(aes(x = q5_rounded, xend = q95_rounded,
                   y = as.character(lab_id), yend = as.character(lab_id),
                   color = type),
               linewidth = 0.1) +
  geom_point(aes(x = median_rounded, y = as.character(lab_id),
                 fill = type, color = type), 
             color = "black", size = 1, shape = 21) + 
  facet_nested(assay_name ~ metric + organ_group, 
               scales = "free", space = "free_y") +
  scale_y_discrete(breaks = as.character(seq(0, 170, 1))) +
  labs(x = "Parameter estimate", y = "Lab") + 
  theme(text = element_text(size = 9),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust = 1),
        legend.position = "top",
        legend.key = element_blank(),
        legend.title = element_blank(),
        panel.spacing = unit(2, "pt"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.y = element_blank(),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.y = element_text(face = "bold", angle = 0, hjust = 0),
        strip.text.x = element_text(face = "bold")
  ); fig


## Save -------------------------------------------------------------------------

ggsave('./outputs/figures/figS12-lab-parameters.png',
       plot = fig,
       width = 7, 
       height = 7,
       dpi = 600)


