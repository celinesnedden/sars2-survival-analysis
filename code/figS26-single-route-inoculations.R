# This file: - Identifies the individuals exposed via single-route inoculation only
#            - Computes the inferred time to detectability & percent positive for these individuals


# Prep -------------------------------------------------------------------------

# Load the event time database
dat <- read.csv("./data/df-event-times.csv")

# Flag tissues as exposed / inoculated or not
dat$inoc_idx <- 0
dat$inoc_idx[dat$dose_nose > 0 & dat$location_grp == "URT" & dat$location_idx == 0] <- 1
dat$inoc_idx[dat$dose_throat > 0 & dat$location_grp == "URT" & dat$location_idx == 1] <- 1
dat$inoc_idx[dat$dose_trachea > 0 & dat$location_grp == "LRT" & dat$location_idx == 0] <- 1
dat$inoc_idx[dat$dose_lung > 0 & dat$location_grp == "LRT" & dat$location_idx == 1] <- 1
dat$inoc_idx[dat$dose_gi > 0 & dat$location_grp == "GI" & dat$location_idx == 0] <- 1

# Consolidate similar assay categories
dat$assay_full <- dat$assay_idx # keep fully stratified assay types for later plotting / analysis
dat$assay_full[is.na(dat$assay_full)] <- -9999 # unknown assays
dat$assay_idx[dat$assay_idx %in% c(3, 4)] <- 3
dat$assay_idx[dat$assay_idx %in% c(5, 6)] <- 4
dat$assay_idx[dat$assay_type == "totRNA" & dat$pcr_target_gene == "Unknown"] <- 1 
dat$assay_idx[dat$assay_type == "sgRNA" & dat$pcr_target_gene == "Unknown"] <- 3

# Distinguish between lab groups for lab / study effects
dat$lab_group <- dat$study_location # groups articles conducted at the same facility (e.g., Tulane Primate Center)
dat$lab_group[str_detect(dat$lab_group, "et al.")] <- dat$article_doi[str_detect(dat$lab_group, "et al.")] # one-off locations or those without sufficient detail for grouping
dat$lab_group <- paste0(dat$article_doi, "-", dat$assay_idx) # lab effects are assay-specific

# Convert TCID50 peak titers to PFU peak titers using standard conversion
dat$peak_observed_titer[!is.na(dat$peak_observed_titer) & 
                          dat$assay_idx == 5 & 
                          !is.na(dat$assay_idx) &
                          dat$peak_observed_titer < 100] <- 
  log10((10^(dat$peak_observed_titer[dat$assay_idx == 5 & 
                                       !is.na(dat$assay_idx) &
                                       dat$peak_observed_titer != 9999 &
                                       !is.na(dat$peak_observed_titer)])) * 0.69)

# Add individual ID numbers for plotting later
dat$indiv_idx <- as.numeric(factor(dat$indiv))

# Data for Stan
dat.stan <- prepare_stan_data_joint(dat)

# Get list of individuals we eventually want to analyze ------------------------

# Only want individuals with ID names & single-route exposures
#dat <- subset(dat, indiv_rep != "No ID")
dat <- subset(dat, inoc_route_subgrp %in% c("IN", "IB", "IT", "IG", "OC"))

## Find individuals with at least two sampled tissues
#multi_tissue_indivs <- c()
#for (indiv.ii in unique(dat$indiv)) {
#  num_tissues <- length(unique(subset(dat, indiv == indiv.ii)$location_subgrp))
#  if (num_tissues > 1) {
#    multi_tissue_indivs <- c(multi_tissue_indivs, indiv.ii)
#  }
#}

# Subset the data to this
#dat <- subset(dat, indiv %in% multi_tissue_indivs)

# Store the indiv id nums
indiv_ids <- unique(dat$indiv_idx)


# Extract true event times ------------------------------------------------------

# Load the model
fit <- readRDS("./outputs/fits/fit-main.RDS")

# Extract the draws from all parameters
df.draws <- fit$draws(format = "df")

# Get only the lab effect parameters
df.draws <- df.draws[, grepl("true_first", names(df.draws))]

# Extract the medians & quantiles
df.summary <- summarise_draws(df.draws, 
                              median, 
                              ~quantile(.x, probs = c(0.05, 0.95)))

# Add the associated assay, tissue, and individual for each of these true times
df.summary$assay_idx <- dat.stan$assay
df.summary$tissue_idx <- dat.stan$tissue_location
df.summary$indiv_idx <- dat.stan$indiv_idx
df.summary$ever_pos <- dat.stan$ever_positive

# Add the upper and lower bounds
df.summary$first_upb <- dat.stan$first_upper_bounds
df.summary$first_lob <- dat.stan$first_lower_bounds

# Subset to the individuals we will consider
df.summary <- subset(df.summary, indiv_idx %in% indiv_ids)

# Get their inoculation route
df.summary$route_name <- NA
for (row_num in 1:nrow(df.summary)) {
  indiv_route <- unique(dat$inoc_route_subgrp[dat$indiv_idx == df.summary$indiv_idx[row_num]])
  if (length(indiv_route) > 1) {print("WARNING")}
  else {df.summary$route_name[row_num] <- indiv_route}
}

# Remove ones with unknown tissue location
df.summary <- subset(df.summary, tissue_idx != -9999)

# Set tissue names
df.summary <- assign_tissue_names(df.summary)

# Factor routes
df.summary$route_name <- factor(df.summary$route_name,
                                levels = c("OC", "IN", "IT", "IB", "IG"))

# Factor the tissues (opposite so they order correctly along y axis)
df.summary$tissue_name <- factor(df.summary$tissue_name,
                                 levels = rev(levels(df.summary$tissue_name)))


# Plot -------------------------------------------------------------------------

## Panel A ---------------------------------------------------------------------

df.counts <- df.summary %>%
  count(tissue_name, route_name)

fig.counts <- ggplot(df.counts) + 
  geom_tile(aes(x = route_name, y = tissue_name, fill = log2(n)), color = "black",
            linewidth = 0.5) +
  geom_text(aes(x = route_name, y = tissue_name, label = n), color = "black") + 
  scale_fill_gradient(low  = "white", high = df.color$Throat) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Exposure Route", y = "Sampled Tissue", fill = "log2(N)") +
  facet_wrap(~ "Sample Size") +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.9, "line"),
        legend.key = element_rect(color = "white"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black")); fig.counts
  

## Panel B ---------------------------------------------------------------------

# Exclude inoculated tissues for clarity
df.summary$median_ref <- df.summary$median

df.summary$median[df.summary$route_name == "OC" & df.summary$tissue_name == "Nose"] <- 0
df.summary$median[df.summary$route_name == "IN" & df.summary$tissue_name == "Nose"] <- 0
df.summary$median[df.summary$route_name == "IT" & df.summary$tissue_name == "Trachea"] <- 0
df.summary$median[df.summary$route_name == "IB" & df.summary$tissue_name == "Lung"] <- 0
df.summary$median[df.summary$route_name == "IG" & df.summary$tissue_name == "Upper GI"] <- 0

df.median <- subset(df.summary, ever_pos == 1) %>%
  group_by(tissue_name, route_name) %>%
  summarise(
    median_value = median(median, na.rm = TRUE),
    median_ref = median(median_ref, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(route_name) %>% 
  mutate(tissue_order = dense_rank(median_value))
  

fig.TD <- ggplot(df.median) + 
  geom_tile(aes(x = route_name, y = tissue_name, fill = tissue_order), color = "black",
            linewidth = 0.5) +
  geom_text(data = subset(df.median, tissue_order != 1),
            aes(x = route_name, y = tissue_name, label = round(median_value, 1)), 
            color = "black") + 
  geom_text(data = subset(df.median, tissue_order == 1),
            aes(x = route_name, y = tissue_name, label = round(median_ref, 1)), 
            color = "black") + 
  scale_fill_gradient(low  = "white", high = df.color$Throat, breaks = 1:6, 
                      labels = c("Inoc.", 1:5),
                      guide = guide_legend(
                        override.aes = list(
                          shape = 14,
                          size  = 2,
                          colour = "black"))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Exposure Route", y = "Sampled Tissue", fill = "Order of\nPositivity") +
  facet_wrap(~ "Time to Detectability") +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.9, "line"),
        legend.key = element_rect(color = "white"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black")); fig.TD


## Panel C ---------------------------------------------------------------------

df.pos <- df.summary %>%
  count(tissue_name, route_name, ever_pos) %>%
  tidyr::complete(tissue_name, route_name, ever_pos, fill = list(n = 0)) %>%
  group_by(tissue_name, route_name) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

fig.percent <- ggplot(subset(df.pos, ever_pos == 1 & !(route_name == "IB" & tissue_name %in% c("Throat", "Trachea", "Upper GI")))) + 
  geom_tile(aes(x = route_name, y = tissue_name, fill = percent), color = "black",
            linewidth = 0.5) +
  geom_text(aes(x = route_name, y = tissue_name, label = round(percent, 1)), 
            color = "black") + 
  scale_fill_gradient(low  = "white", high = df.color$Throat) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Exposure Route", y = "Sampled Tissue", fill = "Percent") +
  facet_wrap(~ "Percent that test positive") +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.9, "line"),
        legend.key = element_rect(color = "white"),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color = "black"),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "transparent", color = "black")); fig.percent


# Combine ---------------------------------------------------------------------

fig <- fig.counts + labs(tag = "a") + fig.TD + labs(tag = "b") + 
  fig.percent + labs(tag = "c") +
  plot_layout(nrow = 1); fig


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS26-single-route-inoculations.png",
       plot = fig,
       width = 9.2, 
       height = 3.6,
       dpi = 600)
