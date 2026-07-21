# This file: - Completes final data prep steps
#            - Runs the Stan model

# Load event times and model ---------------------------------------------------

# Load model 
model <- cmdstanr::cmdstan_model("./code/stan-survival-model.stan")

# Load the event time database
dat <- read.csv("./data/df-event-times.csv")

# Final data tweaks for fitting ------------------------------------------------

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


# Prep data for stan ----------------------------------------------------------

dat.stan <- prepare_stan_data_joint(dat)

# Rescale location-specific doses so they range from 0 to 1
dat.stan$dose_nose <- dat.stan$dose_nose / max(dat.stan$dose_nose)
dat.stan$dose_throat <- dat.stan$dose_throat / max(dat.stan$dose_throat)
dat.stan$dose_trachea <- dat.stan$dose_trachea / max(dat.stan$dose_trachea)
dat.stan$dose_lung <- dat.stan$dose_lung / max(dat.stan$dose_lung)
dat.stan$dose_gi <- dat.stan$dose_gi / max(dat.stan$dose_gi)

# Set unavailable peak times to 9999 so it does not affect fitting (Stan doesn't accept NAs)
dat.stan$peak_observed_time[is.na(dat.stan$peak_observed_time)] <- 9999 

# Save Stan-prepped data for later analysis & plotting
saveRDS(dat.stan, "./data/df-passed-to-Stan.RDS")


# Run the model  ---------------------------------------------------------------

# Set number of iterations
n_iter <- 4000 
n_chains <- 4 

# Model fitting
fit <- model$sample(
  data = dat.stan,
  chains = n_chains,
  iter_warmup = n_iter / 2,
  iter_sampling = n_iter / 2,
  parallel_chains = 6,
  max_treedepth = 12,
  refresh = 20); fit

# Save model
fit$save_object("./outputs/fits/fit-main.RDS")

# Save the model summary
fit.summary <- fit$summary()
write.csv(fit.summary, "./outputs/fits/fit-summary.csv", row.names = FALSE)

