# This file: - Creates Figure 2

## Prep -------------------------------------------------------------------------

# Load the model fit
fit <- readRDS("./outputs/fits/fit-main.RDS")

# Load event times used for fitting
dat <- read.csv("./data/df-event-times.csv")
dat <- assign_all_names(dat)

# Load & prep the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS") 

# Load the predictions
df.estimates <- fread("./outputs/predictions/pred-for-survival-curves-500.csv")

# Subset to total RNA data from IN-exposed individuals sampled in the nose
df <- subset(df.estimates, assay_idx == 1 & route_idx == 1 & tissue_idx == 1)

# Load & prep the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS") 
dat.surv <- as.data.frame(dat.stan)
dat.surv$tissue_idx <- dat.surv$tissue_location
dat.surv <- subset(dat.surv, assay == 1 & route == 1 & tissue_location == 1)
dat.surv <- dat.surv %>%                                          # Arrange for visualization
  arrange(first_lower_bounds, first_upper_bounds)
dat.surv$indiv <- (1:nrow(dat.surv))/nrow(dat.surv)               # Each row is an individual
dat.surv$right_censored[dat.surv$first_upper_bounds == 45] <- "Yes"
dat.surv.right <- subset(dat.surv, right_censored == "Yes")       # Subset by censored type
dat.surv.int <- subset(dat.surv, is.na(right_censored))           # Subset by censored type



## Panel 1: Percent Positive ----------------------------------------------------

# Generate survival curves based on parameter samples from above
df.surv <- get_weibull_curves(df, include_cure = TRUE, 
                              event = "first positive", 
                              x_vals = seq(0, 22, 0.25))

# Plot
fig.percent <- ggplot(subset(df.surv, tissue_idx == 1 & route_idx == 1)) +
  geom_density_ridges(aes(x = percent_positive, y = "BLAH"), 
                      bandwidth = 0.02,
                      fill = "#4A8396", rel_min_height = 0.005) +
  scale_x_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25),
                     limits = c(0, 1)) +
  labs(y = "", x = "% ever detectable") +
  coord_cartesian() +
  coord_flip() +
  theme(text = element_text(size = 13),
        legend.key.size = unit(0.9, "line"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.percent


## Panel 2: First positive ------------------------------------------------------

# Get Median & 90% CI for the curves
df.surv.median <- df.surv %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

# Plot
fig.first <- ggplot() +
  geom_segment(data = dat.surv.int,
               aes(x = first_lower_bounds, 
                   xend = first_upper_bounds,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               color = "#4F609C") +
  geom_point(data = dat.surv.int,
             aes(x = first_upper_bounds, y = indiv),
             alpha = 1, shape = 23,
             fill = "#4F609C") +
  geom_point(data = dat.surv.int,
             aes(x = first_lower_bounds, y = indiv),
             alpha = 1, shape = 21, size = 0.5,
             fill = "black") +
  geom_point(data = dat.surv.right,
             aes(x = first_lower_bounds, y = indiv),
             alpha = 1, shape = 23, size = 1.2,
             color = "#4F609C", fill = "white", stroke = 0.8) +
  geom_segment(data = dat.surv.right,
               aes(x = first_lower_bounds, 
                   xend = 20,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               linetype = "dashed",
               color = "#4F609C") +
  geom_line(data = subset(df.surv.median, tissue_idx == 1),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#4F609C",
            linewidth = 1) +
  geom_ribbon(data = subset(df.surv.median, tissue_idx == 1),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax =cdf_qhigh),
              alpha = 0.4, fill = "#4F609C") +
  coord_cartesian(xlim = c(0, 20)) +
  labs(x = "Days since inoculation", y = "% already detectable") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  theme(text = element_text(size = 13),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.first


## Panel 3: Peak Time -----------------------------------------------------------

# Generate survival curves based on parameter samples from above
df.surv.peak <- get_weibull_curves(df, include_cure = FALSE, 
                                   event = "peak time",
                                   x_vals = seq(0, 22, 0.25))

# Get Median & 90% CI for the curves
df.surv.peak.median <- df.surv.peak %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

# Only include data from individuals with peak information
dat.surv.peak <- subset(dat.surv, has_peak == 1) 

# Adjust times to be relative to first positives times, as we do inside 
#    the Stan model before fitting
dat.surv.peak$peak_lower_bounds_adj <- dat.surv.peak$peak_lower_bounds - 
  dat.surv.peak$first_upper_bounds   
dat.surv.peak$peak_upper_bounds_adj <- dat.surv.peak$peak_upper_bounds - 
  dat.surv.peak$first_lower_bounds  
dat.surv.peak$peak_obs_adj[dat.surv.peak$peak_lower_bounds_adj <= 0] <- 0
dat.surv.peak$peak_obs_adj[dat.surv.peak$peak_lower_bounds_adj > 0] <- 
  dat.surv.peak$peak_upper_bounds_adj[dat.surv.peak$peak_lower_bounds_adj > 0] - (
    dat.surv.peak$peak_upper_bounds[dat.surv.peak$peak_lower_bounds_adj > 0] - 
      dat.surv.peak$peak_observed_time[dat.surv.peak$peak_lower_bounds_adj > 0])
dat.surv.peak$peak_lower_bounds_adj[dat.surv.peak$peak_lower_bounds_adj < 0] <- 0

# Remove leftover individuals without information
dat.surv.peak <- subset(dat.surv.peak, !(peak_lower_bounds_adj == 0 & peak_upper_bounds_adj == 50))

# Arrange by event times for visualization
dat.surv.peak <- dat.surv.peak %>%
  arrange(peak_lower_bounds_adj, peak_obs_adj, peak_upper_bounds_adj)
dat.surv.peak$indiv <- (1:nrow(dat.surv.peak))/nrow(dat.surv.peak)

# Subset by censoring type
dat.surv.peak$right_censored[dat.surv.peak$peak_upper_bounds_adj == 50] <- "Yes"
dat.surv.peak.right <- subset(dat.surv.peak, right_censored == "Yes")
dat.surv.peak.int <- subset(dat.surv.peak, is.na(right_censored))

# Plot
fig.peak <- ggplot() +
  geom_segment(data = dat.surv.peak.int,
               aes(x = peak_lower_bounds_adj, 
                   xend = peak_upper_bounds_adj,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               color = "#9E4472") +
  geom_point(data = dat.surv.peak.int,
             aes(x = peak_upper_bounds_adj, y = indiv),
             alpha = 1, shape = 21, size = 0.5,
             fill = "black") +
  geom_point(data = dat.surv.peak.int,
             aes(x = peak_lower_bounds_adj, y = indiv),
             alpha = 1, shape = 21, size = 0.5,
             fill = "black") +
  geom_point(data = dat.surv.peak.int,
             aes(x = peak_obs_adj, y = indiv),
             alpha = 1, shape = 24, 
             fill = "#9E4472") +
  geom_point(data = dat.surv.peak.right,
             aes(x = peak_lower_bounds_adj, y = indiv),
             alpha = 1, shape = 24, size = 1.2,
             color = "#9E4472", fill = "white", stroke = 0.8) +
  geom_segment(data = dat.surv.peak.right,
               aes(x =  peak_lower_bounds_adj, 
                   xend = 20,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               linetype = "dashed",
               color = "#9E4472") +
  geom_line(data = subset(df.surv.peak.median, tissue_idx == 1),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#9E4472",
            linewidth = 1) +
  geom_ribbon(data = subset(df.surv.peak.median, tissue_idx == 1),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax =cdf_qhigh),
              alpha = 0.4, fill = "#9E4472") +
  coord_cartesian(xlim = c(0, 20)) +
  labs(x = "Days since detectability", y = "% already peaked") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  theme(text = element_text(size = 13),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.peak


## Panel 4: Peak Titer ----------------------------------------------------------

# Get observed peak titers from those with that information
observed_peaks <- dat.stan$peak_observed_titer[dat.stan$has_titer == 1]

# Get predicted peak titers from the corresponding individuals 
rows <- which(dat.stan$has_titer == 1) # corresponding rows
true_peaks <- fit$summary(c(paste0("true_peak_titer[", rows, "]")))[2] # model predictions

# Combine into a dataframe for plotting purposes
df.peak.titer <- data.frame(observed = observed_peaks,
                            true = true_peaks,
                            assay = dat.stan$assay[dat.stan$has_titer == 1])

# Get correlation, to report in MS, for PCR and culture individually
cor(df.peak.titer$mean, df.peak.titer$observed) # overall
cor(df.peak.titer$mean[df.peak.titer$assay < 4], df.peak.titer$observed[df.peak.titer$assay < 4])
cor(df.peak.titer$mean[df.peak.titer$assay == 4], df.peak.titer$observed[df.peak.titer$assay == 4])

# Get difference, to report in MS, for PCR and culture individually
mean(df.peak.titer$mean - df.peak.titer$observed) # overall
mean(df.peak.titer$mean[df.peak.titer$assay < 4] - df.peak.titer$observed[df.peak.titer$assay < 4])
mean(df.peak.titer$mean[df.peak.titer$assay == 4] - df.peak.titer$observed[df.peak.titer$assay == 4])

# Plot
fig.titer <- ggplot(data = subset(df.peak.titer, assay < 4)) +
  geom_point(aes(x = observed, y = mean),
             alpha = 0.3, stroke = 0.01, size = 1.5,
             shape = 24, fill = "#C87EA4") + 
  scale_x_continuous(breaks = seq(0, 14, 4),
                     labels = c(expression(paste(10^0)),
                                #expression(paste(10^2)), 
                                expression(paste(10^4)), 
                                #expression(paste(10^6)), 
                                expression(paste(10^8)), 
                                #expression(paste(10^10)),
                                expression(paste(10^12))#, 
                                #expression(paste(10^14))
                     ),
                     limits = c(-1, 13)) +
  scale_y_continuous(breaks = seq(0, 14, 4),
                     labels = c(expression(paste(10^0)),
                                #expression(paste(10^2)), 
                                expression(paste(10^4)), 
                                #expression(paste(10^6)), 
                                expression(paste(10^8)), 
                                #expression(paste(10^10)),
                                expression(paste(10^12))#, 
                                #expression(paste(10^14))
                     ),
                     limits = c(-1, 13)) +
  labs(x = "Largest observed titer",
       y = "Predicted peak titer") + 
  theme(text = element_text(size = 13),
        legend.position = "right",
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 11, face = "bold"),
        strip.text.y = element_blank(),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.titer

inset <- ggplot(data = subset(df.peak.titer, assay < 4))+
  geom_density(aes(x = mean - observed), color = "black",
               fill = "#C87EA4") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  labs(x = "Difference", y = "Density") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.9),
                     breaks = seq(0, 0.8, 0.4)) +
  theme(text = element_text(size = 9),
        legend.position = "right",
        legend.key.size = unit(0.9, "line"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); inset


fig.titer.combined <- 
  fig.titer +
  annotation_custom(
    grob = ggplotGrob(inset),
    xmin = 6.2, xmax = 13.5,
    ymin = -1.5, ymax = 5.2
  ); fig.titer.combined


## Panel 5: Last positive ------------------------------------------------------------------

# Generate survival curves based on parameter samples from above
df.surv.last <- get_weibull_curves(df, include_cure = FALSE, 
                                   event = "last positive",
                                   x_vals = seq(0, 36, 0.25))

# Get Median & 90% CI for the curves
df.surv.last.median <- df.surv.last %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

# Only include data from individuals with last positive information
dat.surv.last <- subset(dat.surv, has_last_positive == 1)

# Adjust times to be relative to peak times, as we do inside 
#    the Stan model before fitting
dat.surv.last$last_lower_bounds_adj <- dat.surv.last$last_lower_bounds - dat.surv.last$peak_upper_bounds
dat.surv.last$last_upper_bounds_adj <- dat.surv.last$last_upper_bounds - dat.surv.last$peak_lower_bounds
dat.surv.last$last_lower_bounds_adj[dat.surv.last$last_lower_bounds_adj < 0] <- 0

# Remove leftover individuals without information
dat.surv.last <- subset(dat.surv.last, !(last_lower_bounds_adj == 0 & last_upper_bounds_adj > 50))

# Arrange by event times for visualization
dat.surv.last <- dat.surv.last %>%
  arrange(last_lower_bounds_adj, last_upper_bounds_adj)
dat.surv.last$indiv <- (1:nrow(dat.surv.last))/nrow(dat.surv.last)

# Stratify by censoring type
dat.surv.last$right_censored[dat.surv.last$last_upper_bounds_adj > 50] <- "Yes"
dat.surv.last.right <- subset(dat.surv.last, right_censored == "Yes")
dat.surv.last.int <- subset(dat.surv.last, is.na(right_censored))


fig.last <- ggplot() +
  geom_segment(data = dat.surv.last.int,
               aes(x = last_lower_bounds_adj, 
                   xend = last_upper_bounds_adj,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               color = "#628D56") +
  geom_point(data = dat.surv.last.int,
             aes(x = last_upper_bounds_adj, y = indiv),
             alpha = 1, shape = 21, size = 0.5,
             fill = "black") +
  geom_point(data = dat.surv.last.int,
             aes(x = last_lower_bounds_adj, y = indiv),
             alpha = 1, shape = 22, 
             fill = "#628D56") +
  geom_point(data = dat.surv.last.right,
             aes(x = last_lower_bounds_adj, y = indiv),
             alpha = 1, shape = 22, size = 1.2,
             color = "#628D56", fill = "white", stroke = 0.8) +
  geom_segment(data = dat.surv.last.right,
               aes(x =  last_lower_bounds_adj, 
                   xend = 38,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               linetype = "dashed",
               color = "#628D56") +
  geom_line(data = subset(df.surv.last.median, tissue_idx == 1),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#628D56",
            linewidth = 1) +
  geom_ribbon(data = subset(df.surv.last.median, tissue_idx == 1),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax =cdf_qhigh),
              alpha = 0.4, fill = "#628D56") +
  coord_cartesian(xlim = c(0, 34)) +
  labs(x = "Days since peak", y = "% already undetectabile") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 40, 8)) +
  theme(text = element_text(size = 13),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.last


## Combine --------------------------------------------------------------------

fig.comb <- 
  fig.percent + labs(tag = "A") + 
  fig.first + labs(tag = "B") +
  fig.peak + labs(tag = "C") + 
  fig.titer.combined  + labs(tag = "D") + 
  fig.last + labs(tag = "E") +
  plot_layout(nrow = 1, 
              widths = c(1, 6, 6, 6, 6)); fig.comb

## Save --------------------------------------------------------------------

ggsave('./outputs/figures/fig2-example-model-fits.pdf',
       plot = fig.comb,
       width = 13, 
       height = 3,
       dpi = 600)


