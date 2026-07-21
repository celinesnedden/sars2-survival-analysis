# This file: - Compares the predicted and observed peak titers across individuals & sample locations


# Prep -------------------------------------------------------------------------

# Load the model fit
fit <- readRDS("./outputs/fits/fit-main.RDS")

# Load & prep the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS") 

# Get observed peak titers from those with that information
observed_peaks <- dat.stan$peak_observed_titer[dat.stan$has_titer == 1]

# Get predicted peak titers from the corresponding individuals 
rows <- which(dat.stan$has_titer == 1) # corresponding rows
true_peaks <- fit$summary(c(paste0("true_peak_titer[", rows, "]")))[2] # model predictions

# Get their assay types
assay_types <- dat.stan$assay[dat.stan$has_titer == 1]
assay_full <- dat.stan$assay_full[dat.stan$has_titer == 1]


# Combine into a dataframe for plotting purposes
df <- data.frame(observed = observed_peaks,
                 true = true_peaks,
                 assay_type = assay_types,
                 assay_full = assay_full)

# Add assay nammes
df$assay_names[df$assay_type == 1] <- "Total RNA"
df$assay_names[df$assay_type == 2] <- "gRNA"
df$assay_names[df$assay_type == 3] <- "sgRNA"
df$assay_names[df$assay_type == 4] <- "Culture (PFU)"
df$assay_names[df$assay_type == -9999] <- "Unknown"
df$assay_names <- factor(df$assay_names, levels = c("Total RNA", 
                                                    "gRNA", "sgRNA", 
                                                    "Culture (PFU)",
                                                    "Unknown"))

# Get correlation, to report in MS
cor(df$mean[df$assay_names == "Culture (PFU)"], df$observed[df$assay_names == "Culture (PFU)"])
mean(df$mean[df$assay_names == "Culture (PFU)"] - df$observed[df$assay_names == "Culture (PFU)"])

cor(df$mean[df$assay_names != "Culture (PFU)"], df$observed[df$assay_names != "Culture (PFU)"])
mean(df$mean[df$assay_names != "Culture (PFU)"] - df$observed[df$assay_names != "Culture (PFU)"])


# Plot -------------------------------------------------------------------------

fig.titer <- ggplot(data = df) +
  geom_point(aes(x = observed, y = mean), #fill = assay_names),
             fill = "#C87EA4",
             alpha = 0.2, stroke = 0.2, 
             size = 2, shape = 24) + 
  scale_x_continuous(breaks = seq(0, 14, 4),
                     labels = c(expression(paste(10^0)),
                                expression(paste(10^4)),  
                                expression(paste(10^8)), 
                                expression(paste(10^12))
                     )) +
  scale_y_continuous(breaks = seq(0, 14, 4),
                     labels = c(expression(paste(10^0)),
                                expression(paste(10^4)), 
                                expression(paste(10^8)), 
                                expression(paste(10^12))
                     )) +
  facet_wrap(.~ assay_names, nrow = 1) +
  labs(x = "Largest observed titer",
       y = "Predicted peak titer", fill = "Assay") + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(text = element_text(size = 12),
        legend.position = c(0.8, 0.27),
        legend.key = element_rect(color = "transparent"),
        legend.key.size = unit(0.9, "line"),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(face = "bold"),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.titer



# Save -------------------------------------------------------------------------

ggsave('./outputs/figures/figS8-predicted-vs-observed-peak-titers.png',
       plot = fig.titer,
       width = 8, 
       height = 2.2,
       dpi = 600)
