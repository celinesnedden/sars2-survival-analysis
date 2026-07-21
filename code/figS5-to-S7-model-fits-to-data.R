# This file creates: model fits by route, tissue, and assay
# Dependencies: 

# Prep -------------------------------------------------------------------------

# Load event times used for fitting
dat <- read.csv("./data/df-event-times.csv")
dat <- assign_all_names(dat)

# Load & prep the data passed to Stan
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS") 

# Load the predictions
df.estimates <- fread("./outputs/predictions/pred-for-survival-curves-500.csv")


# Figure S5: Time to detectability ---------------------------------------------

## Prep -----------------------------------------------------------------------

# Get the curves for each route
df.surv.in <- get_weibull_curves(subset(df.estimates, route_idx == 1), 
                                 include_cure = TRUE, 
                                 event = "first positive",
                                 x_vals = seq(0, 20, 1))

df.surv.it <- get_weibull_curves(subset(df.estimates, route_idx == 2), 
                                 include_cure = TRUE, 
                                 event = "first positive",
                                 x_vals = seq(0, 20, 1))

df.surv.init <- get_weibull_curves(subset(df.estimates, route_idx == 3), 
                                   include_cure = TRUE, 
                                   event = "first positive",
                                   x_vals = seq(0, 20, 1))

df.surv.ae <- get_weibull_curves(subset(df.estimates, route_idx == 4), 
                                 include_cure = TRUE, 
                                 event = "first positive", 
                                 x_vals = seq(0, 20, 1))

df.surv.ig <- get_weibull_curves(subset(df.estimates, route_idx == 5), 
                                 include_cure = TRUE, 
                                 event = "first positive",
                                 x_vals = seq(0, 20, 1))


# Get the medians for each route
df.surv.median.in <- df.surv.in %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.median.it <- df.surv.it %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.median.init <- df.surv.init %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.median.ae <- df.surv.ae %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.median.ig <- df.surv.ig %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

# Combine them together
df.surv.median <- rbind(df.surv.median.in, df.surv.median.it,
                        df.surv.median.init, df.surv.median.ae,
                        df.surv.median.ig)


# Assign cofactor names
df.surv.median <- assign_route_names(df.surv.median)
df.surv.median <- assign_tissue_names(df.surv.median)
df.surv.median$assay_name[df.surv.median$assay_idx == 1] <- "PCR"
df.surv.median$assay_name[df.surv.median$assay_idx == 4] <- "Culture"

# Load dat.stan from model file to get data in correct form
dat.surv <- as.data.frame(dat.stan)
dat.surv$tissue_idx <- dat.surv$tissue_location
dat.surv <- subset(dat.surv, tissue_location %in% c(1, 2, 3, 4, 6))
dat.surv <- assign_tissue_names(dat.surv)
dat.surv$route_idx <- dat.surv$route
dat.surv <- assign_route_names(dat.surv)
dat.surv$assay_name[dat.surv$assay %in% c(1, 2, 3, -9999)] <- "PCR"
dat.surv$assay_name[dat.surv$assay == 4] <- "Culture"

# Order them based on bounds and censor type
dat.surv$indiv <- NA
dat.surv.ordered <- dat.surv[0, ]

for (route.ii in unique(dat.surv$route_idx)) {
  dat.surv.route <- subset(dat.surv, route_idx == route.ii)
  for (assay.ii in unique(dat.surv.route$assay_name)) {
    dat.surv.assay <- subset(dat.surv.route, assay_name == assay.ii)
    for (tissue.ii in unique(dat.surv.assay$tissue_idx)) {
      dat.surv.tissue <- subset(dat.surv.assay, tissue_idx == tissue.ii)
      
      dat.surv.tissue <- dat.surv.tissue %>%
        arrange(first_lower_bounds, first_upper_bounds)
      
      dat.surv.tissue$indiv <- (1:nrow(dat.surv.tissue))/nrow(dat.surv.tissue)
      
      dat.surv.ordered <- rbind(dat.surv.ordered, dat.surv.tissue)
    }
  }
}
dat.surv <- dat.surv.ordered
dat.surv$right_censored[dat.surv$first_upper_bounds == 45] <- "Yes"
dat.surv.right <- subset(dat.surv, right_censored == "Yes")
dat.surv.int <- subset(dat.surv, is.na(right_censored))


## Plot -----------------------------------------------------------------------

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
  geom_point(data = subset(dat.surv.right, first_lower_bounds <= 20),
             aes(x = first_lower_bounds, y = indiv),
             alpha = 1, shape = 23, size = 1.2,
             color = "#4F609C", fill = "white", stroke = 0.8) +
  geom_segment(data = subset(dat.surv.right, first_lower_bounds <= 20),
               aes(x = first_lower_bounds, 
                   xend = 20,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               linetype = "dashed",
               color = "#4F609C") +
  geom_point(data = subset(dat.surv.right, first_lower_bounds > 20),
             aes(x = 21, y = indiv),
             alpha = 1, shape = 23, size = 1.2,
             color = "#4F609C", fill = "white", stroke = 0.8) +
  geom_line(data = subset(df.surv.median),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#4F609C",
            linewidth = 1) +
  geom_ribbon(data = subset(df.surv.median),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax =cdf_qhigh),
              alpha = 0.4, fill = "#4F609C") +
  facet_grid(tissue_name+assay_name ~ route_name) +
  coord_cartesian(xlim = c(0, 20), clip = "off") +
  labs(x = "days since inoculation", y = "% already detectable") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.first

ggsave('./outputs/figures/figS5-model-fits-percent-and-detectability.png',
       plot = fig.first,
       width = 8.5, 
       height = 11,
       dpi = 600)



# Figure S6: Time to peak titer ------------------------------------------------

## Prep -----------------------------------------------------------------------

# Get the curves for each route
df.surv.peak.in <- get_weibull_curves(subset(df.estimates, route_idx == 1), 
                                      include_cure = FALSE, 
                                      event = "peak time", 
                                      x_vals = seq(0, 20, 1))

df.surv.peak.it <- get_weibull_curves(subset(df.estimates, route_idx == 2), 
                                      include_cure = FALSE, 
                                      event = "peak time", 
                                      x_vals = seq(0, 20, 1))

df.surv.peak.init <- get_weibull_curves(subset(df.estimates, route_idx == 3), 
                                      include_cure = FALSE, 
                                      event = "peak time", 
                                      x_vals = seq(0, 20, 1))

df.surv.peak.ae <- get_weibull_curves(subset(df.estimates, route_idx == 4), 
                                      include_cure = FALSE, 
                                      event = "peak time", 
                                      x_vals = seq(0, 20, 1))

df.surv.peak.ig <- get_weibull_curves(subset(df.estimates, route_idx == 5), 
                                      include_cure = FALSE, 
                                      event = "peak time", 
                                      x_vals = seq(0, 20, 1))


# Get the medians for each route
df.surv.peak.median.in <- df.surv.peak.in %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

df.surv.peak.median.it <- df.surv.peak.it %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

df.surv.peak.median.init <- df.surv.peak.init %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

df.surv.peak.median.ae <- df.surv.peak.ae %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

df.surv.peak.median.ig <- df.surv.peak.ig %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


# Combine them together
df.surv.peak.median <- rbind(df.surv.peak.median.in, df.surv.peak.median.it,
                             df.surv.peak.median.init, df.surv.peak.median.ae,
                             df.surv.peak.median.ig)

# Assign cofactor names
df.surv.peak.median <- assign_route_names(df.surv.peak.median)
df.surv.peak.median <- assign_tissue_names(df.surv.peak.median)
df.surv.peak.median$assay_name[df.surv.peak.median$assay_idx == 1] <- "PCR"
df.surv.peak.median$assay_name[df.surv.peak.median$assay_idx == 4] <- "Culture"

# Get data for peak into basic correct form
dat.surv.peak <- as.data.frame(dat.stan)
dat.surv.peak <- subset(dat.surv.peak, has_peak == 1)
dat.surv.peak$tissue_idx <- dat.surv.peak$tissue_location
dat.surv.peak <- subset(dat.surv.peak, tissue_location %in% c(1, 2, 3, 4, 6))
dat.surv.peak <- assign_tissue_names(dat.surv.peak)
dat.surv.peak$route_idx <- dat.surv.peak$route
dat.surv.peak <- assign_route_names(dat.surv.peak)
dat.surv.peak$assay_name[dat.surv.peak$assay %in% c(1, 2, 3, -9999)] <- "PCR"
dat.surv.peak$assay_name[dat.surv.peak$assay == 4] <- "Culture"

# Adjust the lower & upper bounds based on other time to detectability
dat.surv.peak$peak_lower_bounds_adj <- dat.surv.peak$peak_lower_bounds - dat.surv.peak$first_upper_bounds
dat.surv.peak$peak_upper_bounds_adj <- dat.surv.peak$peak_upper_bounds - dat.surv.peak$first_lower_bounds
dat.surv.peak$peak_obs_adj[dat.surv.peak$peak_lower_bounds_adj <= 0] <- 0
dat.surv.peak$peak_obs_adj[dat.surv.peak$peak_lower_bounds_adj > 0] <- 
  dat.surv.peak$peak_upper_bounds_adj[dat.surv.peak$peak_lower_bounds_adj > 0] - (
    dat.surv.peak$peak_upper_bounds[dat.surv.peak$peak_lower_bounds_adj > 0] - 
      dat.surv.peak$peak_observed_time[dat.surv.peak$peak_lower_bounds_adj > 0]
  )
dat.surv.peak$peak_lower_bounds_adj[dat.surv.peak$peak_lower_bounds_adj < 0] <- 0
dat.surv.peak <- subset(dat.surv.peak, !(peak_lower_bounds_adj == 0 & peak_upper_bounds_adj == 50))

# Order them by bounds and censor type
dat.surv.peak$indiv <- NA
dat.surv.peak.ordered <- dat.surv.peak[0, ]
for (route.ii in unique(dat.surv.peak$route_idx)) {
  dat.surv.peak.route <- subset(dat.surv.peak, route_idx == route.ii)
  for (assay.ii in unique(dat.surv.peak.route$assay_name)) {
    dat.surv.peak.assay <- subset(dat.surv.peak.route, assay_name == assay.ii)
    for (tissue.ii in unique(dat.surv.peak.assay$tissue_idx)) {
      dat.surv.peak.tissue <- subset(dat.surv.peak.assay, tissue_idx == tissue.ii)
      
      dat.surv.peak.tissue <- dat.surv.peak.tissue %>%
        arrange(peak_lower_bounds_adj, peak_obs_adj, peak_upper_bounds_adj)
      
      dat.surv.peak.tissue$indiv <- (1:nrow(dat.surv.peak.tissue))/nrow(dat.surv.peak.tissue)
      
      dat.surv.peak.ordered <- rbind(dat.surv.peak.ordered, dat.surv.peak.tissue)
    }
  }
}


dat.surv.peak <- dat.surv.peak.ordered
dat.surv.peak$right_censored[dat.surv.peak$peak_upper_bounds_adj == 50] <- "Yes"
dat.surv.peak$peak_upper_bounds_adj[dat.surv.peak$peak_upper_bounds_adj>20] <- 21

# Subset by censor type
dat.surv.peak.right <- subset(dat.surv.peak, right_censored == "Yes")
dat.surv.peak.int <- subset(dat.surv.peak, is.na(right_censored))


## Plot -----------------------------------------------------------------------

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
  geom_line(data = subset(df.surv.peak.median),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#9E4472",
            linewidth = 1) +s
  geom_ribbon(data = subset(df.surv.peak.median),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax =cdf_qhigh),
              alpha = 0.4, fill = "#9E4472") +
  facet_grid(tissue_name+assay_name ~ route_name) +
  coord_cartesian(xlim = c(0, 20), clip = "off") +
  labs(x = "days since detectability", y = "% already peaked") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 20, 4)) +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.peak



ggsave('./outputs/figures/figS6-model-fits-peak-time.png',
       plot = fig.peak,
       width = 8.5, 
       height = 11,
       dpi = 600)



# Figure S7: Time to undetectability -------------------------------------------

## Prep -----------------------------------------------------------------------

# Get the curves for each route
df.surv.last.in <- get_weibull_curves(subset(df.estimates, route_idx == 1), 
                                      include_cure = FALSE, 
                                      event = "last positive",
                                      x_vals = seq(0, 25, 1))

df.surv.last.it <- get_weibull_curves(subset(df.estimates, route_idx == 2), 
                                      include_cure = FALSE, 
                                      event = "last positive",
                                      x_vals = seq(0, 25, 1))

df.surv.last.init <- get_weibull_curves(subset(df.estimates, route_idx == 3), 
                                        include_cure = FALSE, 
                                        event = "last positive",
                                        x_vals = seq(0, 25, 1))

df.surv.last.ae <- get_weibull_curves(subset(df.estimates, route_idx == 4), 
                                      include_cure = FALSE, 
                                      event = "last positive",
                                      x_vals = seq(0, 25, 1))

df.surv.last.ig <- get_weibull_curves(subset(df.estimates, route_idx == 5), 
                                      include_cure = FALSE, 
                                      event = "last positive",
                                      x_vals = seq(0, 25, 1))

# Get the medians
df.surv.last.median.in <- df.surv.last.in %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.last.median.it <- df.surv.last.it %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.last.median.init <- df.surv.last.init %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.last.median.ae <- df.surv.last.ae %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 


df.surv.last.median.ig <- df.surv.last.ig %>%
  group_by(route_idx, tissue_idx, assay_idx, day_post_infection) %>%
  summarise(across(c(cdf), 
                   list(
                     median = ~median(., na.rm = TRUE),
                     qlow = ~quantile(., probs = 0.05, na.rm = TRUE),
                     qhigh = ~quantile(., probs = 0.95, na.rm = TRUE)))) 

# Bind them together
df.surv.last.median <- rbind(df.surv.last.median.in, df.surv.last.median.it,
                             df.surv.last.median.init, df.surv.last.median.ae,
                             df.surv.last.median.ig)

# Assign cofactor names
df.surv.last.median <- assign_route_names(df.surv.last.median)
df.surv.last.median <- assign_tissue_names(df.surv.last.median)
df.surv.last.median$assay_name[df.surv.last.median$assay_idx == 1] <- "PCR"
df.surv.last.median$assay_name[df.surv.last.median$assay_idx == 4] <- "Culture"

# Set up data for plotting
dat.surv.last <- as.data.frame(dat.stan)
dat.surv.last <- subset(dat.surv.last, has_last_positive == 1)
dat.surv.last$tissue_idx <- dat.surv.last$tissue_location
dat.surv.last <- subset(dat.surv.last, tissue_location %in% c(1, 2, 3, 4, 6))
dat.surv.last <- assign_tissue_names(dat.surv.last)
dat.surv.last$route_idx <- dat.surv.last$route
dat.surv.last <- assign_route_names(dat.surv.last)
dat.surv.last$assay_name[dat.surv.last$assay %in% c(1, 2, 3, -9999)] <- "PCR"
dat.surv.last$assay_name[dat.surv.last$assay == 4] <- "Culture"

dat.surv.last$last_lower_bounds_adj <- dat.surv.last$last_lower_bounds - dat.surv.last$peak_upper_bounds
dat.surv.last$last_upper_bounds_adj <- dat.surv.last$last_upper_bounds - dat.surv.last$peak_lower_bounds
dat.surv.last$last_lower_bounds_adj[dat.surv.last$last_lower_bounds_adj < 0] <- 0

dat.surv.last <- subset(dat.surv.last, !(last_lower_bounds_adj == 0 & last_upper_bounds_adj > 40))

dat.surv.last$indiv <- NA
dat.surv.last.ordered <- dat.surv.last[0, ]

for (route.ii in unique(dat.surv.last$route_idx)) {
  dat.surv.last.route <- subset(dat.surv.last, route_idx == route.ii)
  for (assay.ii in unique(dat.surv.last.route$assay_name)) {
    dat.surv.last.assay <- subset(dat.surv.last.route, assay_name == assay.ii)
    for (tissue.ii in unique(dat.surv.last.assay$tissue_idx)) {
      dat.surv.last.tissue <- subset(dat.surv.last.assay, tissue_idx == tissue.ii)
      
      dat.surv.last.tissue <- dat.surv.last.tissue %>%
        arrange(last_lower_bounds_adj, last_upper_bounds_adj)
      
      dat.surv.last.tissue$indiv <- (1:nrow(dat.surv.last.tissue))/nrow(dat.surv.last.tissue)
      
      dat.surv.last.ordered <- rbind(dat.surv.last.ordered, dat.surv.last.tissue)
    }
  }
}
dat.surv.last <- dat.surv.last.ordered

dat.surv.last$last_upper_bounds_adj[dat.surv.last$last_upper_bounds_adj > 24] <- 24.5
dat.surv.last$last_lower_bounds_adj[dat.surv.last$last_lower_bounds_adj > 24] <- 24.5

dat.surv.last$right_censored[dat.surv.last$last_upper_bounds_adj > 50] <- "Yes"

dat.surv.last.right <- subset(dat.surv.last, right_censored == "Yes")
dat.surv.last.int <- subset(dat.surv.last, is.na(right_censored))


## Plot -----------------------------------------------------------------------

fig.last <- ggplot() +
  coord_cartesian(xlim = c(0, 24), clip = "off") +
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
                   xend = 25,
                   y = indiv, yend = indiv),
               alpha = 0.5, linewidth = 0.2,
               linetype = "dashed",
               color = "#628D56") +
  geom_line(data = subset(df.surv.last.median),
            aes(x = day_post_infection, y = cdf_median),
            alpha = 1, color = "#628D56",
            linewidth = 1) +
  geom_ribbon(data = subset(df.surv.last.median),
              aes(x = day_post_infection, ymin = cdf_qlow, 
                  ymax = cdf_qhigh),
              alpha = 0.4, fill = "#628D56") +
  facet_grid(tissue_name+assay_name ~ route_name) +
  
  labs(x = "days since peak", y = "% already undetectabile") +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = seq(0, 100, 25)) +
  scale_x_continuous(breaks = seq(0, 40, 8)) +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "light grey"), 
        panel.grid.minor = element_line(size = 0.08, linetype = 'solid',
                                        colour = "light grey")); fig.last


ggsave('./outputs/figures/figS7-model-fits-undetectability.png',
       plot = fig.last,
       width = 8.5, 
       height = 11,
       dpi = 600)






