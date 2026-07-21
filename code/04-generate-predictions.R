# This file: - Generates predictions for most subsequent analyses in figure files


# Prep -------------------------------------------------------------------------

# Load the model fit
fit <- readRDS("./outputs/fits/fit-main.RDS")

# Set maximum location-specific doses for rescaling
#  (note these were computed and taken from the 03-model-fitting file)
max_dose_nose <- 6.821186 
max_dose_throat <- 6.90768 
max_dose_trachea <- 7.350608
max_dose_lung <- 6.562293
max_dose_gi <- 7

# Get ready for parallel processing
library(doParallel)
registerDoParallel(cores = 6)


# Generate predictions for all cofactors ---------------------------------------

# Set the # of posterior draws per run
n_samps <- 50
n_runs <- 1000 / n_samps

# Set the seed for sampling the posterior
seed <- 13589

# Start by generating predictions for 10^1
df.pred <- get_metrics_across_cofactors_parallel(fit, n_samps, 
                                                 seed = seed,
                                                 tissue_options = 1:6,
                                                 route_options = 1:5,
                                                 dose_options = c(1.2, 7.4, 1:7),
                                                 assay_options = 1:4,
                                                 age_options = 1:3,
                                                 species_options = 1:3,
                                                 sex_options = c(0, 1),
                                                 rescale_doses = TRUE)

# Loop to generate predictions for other doses
for (run_num in 1:(n_runs - 1)) {
  df.pred.ii <- get_metrics_across_cofactors_parallel(fit, n_samps, 
                                                      seed = seed + run_num,
                                                      tissue_options = 1:6,
                                                      route_options = 1:5,
                                                      dose_options =  c(1.2, 7.4, 1:7),
                                                      assay_options = 1:4,
                                                      age_options = 1:3,
                                                      species_options = 1:3,
                                                      sex_options = c(0, 1),
                                                      rescale_doses = TRUE)
  
  df.pred <- rbind(df.pred, df.pred.ii)
  
  # Write it temporarily to save in case something fails
  write.csv(df.pred, 
            "./outputs/predictions/pred-across-cofactors-1000.csv",
            row.names = FALSE)
}

# Update calendar times and AUC
df.pred$duration_median <- df.pred$peak_median + df.pred$last_median
df.pred$last_median <- df.pred$first_pos_median + df.pred$peak_median + df.pred$last_median
df.pred$peak_median <- df.pred$first_pos_median + df.pred$peak_median
df.pred$auc <- abs(df.pred$auc) # take the abs value since some studies reported negative titer vals

# Save these predictions
write.csv(df.pred, 
          "./outputs/predictions/pred-across-cofactors-1000.csv",
          row.names = FALSE)


# Generate predictions for dose effects ----------------------------------------

# Set the # of posterior draws per run
n_samps <- 50
n_runs <- 500 / n_samps

# Set the seed for sampling the posterior
seed <- 13589

# Set doses to run
#doses <- c(seq(1.3, 1.9, 0.1), seq(2.1, 2.9, 0.1), seq(3.1, 3.9, 0.1), 
#           seq(4.1, 4.9, 0.1), seq(5.1, 5.9, 0.1),  seq(6.1, 6.9, 0.1), 
#           7.1, 7.2, 7.3)

# Start by generating predictions the first 50 samples
df.pred.dose <- get_metrics_across_cofactors_parallel(fit, n_samps, 
                                                      seed = seed,
                                                      tissue_options = c(1, 4, 6),
                                                      route_options = 1:5,
                                                      dose_options = seq(1.2, 7.4, 0.1),
                                                      assay_options = c(1, 4),
                                                      age_options = 1:3,
                                                      species_options = 1:3,
                                                      sex_options = c(0, 1),
                                                      rescale_doses = TRUE)

# Loop to generate predictions for other doses
for (run_num in 1:(n_runs - 1)) {
  df.pred.dose.ii <- get_metrics_across_cofactors_parallel(fit, n_samps, 
                                                           seed = seed + run_num,
                                                           tissue_options = c(1, 4, 6),
                                                           route_options = 1:5,
                                                           dose_options = seq(1.2, 7.4, 0.1),
                                                           assay_options = c(1, 4),
                                                           age_options = 1:3,
                                                           species_options = 1:3,
                                                           sex_options = c(0, 1),
                                                           rescale_doses = TRUE)
  
  df.pred.dose <- rbind(df.pred.dose, df.pred.dose.ii)
  
  # Write it temporarily to save in case something fails
  write.csv(df.pred.dose, 
            "./outputs/predictions/pred-across-doses-500.csv",
            row.names = FALSE)
}

# Update calendar times and AUC
df.pred.dose$duration_median <- df.pred.dose$peak_median + df.pred.dose$last_median
df.pred.dose$last_median <- df.pred.dose$first_pos_median + df.pred.dose$peak_median + df.pred.dose$last_median
df.pred.dose$peak_median <- df.pred.dose$first_pos_median + df.pred.dose$peak_median
df.pred.dose$auc <- abs(df.pred.dose$auc) # take the abs value since some studies reported negative titer vals

# Save these predictions
write.csv(df.pred.dose, 
          "./outputs/predictions/pred-across-doses-450.csv",
          row.names = FALSE)

saveRDS(df.pred.dose, 
          "./outputs/predictions/pred-across-doses-450.RDS")


# Pairwise differences, all cofactors, higher doses ---------------------------

df.pred <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Mean difference, 10^4 and 10^7 only 
df.pred.main <- subset(df.pred, dose_total %in% c(4, 7))
colnames(df.pred.main)[colnames(df.pred.main) == "duration_median"] <- "duration"

# Calculate all possible pairwise differences
route.df <- get_pairwise_differences(df.pred.main, cofactor = "route")
dose.df <- get_pairwise_differences(df.pred.main, cofactor = "dose")
sex.df <- get_pairwise_differences(df.pred.main, cofactor = "sex")
age.df <- get_pairwise_differences(df.pred.main, cofactor = "age")
species.df <- get_pairwise_differences(df.pred.main, cofactor = "species")

# Save these pairwise differences. Very large so saving as RDS
saveRDS(route.df, "./outputs/predictions/pred-route-differences-1000.RDS")
saveRDS(dose.df, "./outputs/predictions/pred-dose-differences-1000.RDS")
saveRDS(sex.df, "./outputs/predictions/pred-sex-differences-1000.RDS")
saveRDS(age.df, "./outputs/predictions/pred-age-differences-1000.RDS")
saveRDS(species.df, "./outputs/predictions/pred-species-differences-1000.RDS")


# Pairwise differences, all cofactors, lower doses ----------------------------

df.pred <- fread("./outputs/predictions/pred-across-cofactors-1000.csv")

# Mean difference, 10^1 and 10^4 only 
df.pred.main <- subset(df.pred, dose_total %in% c(1, 4, 7))
colnames(df.pred.main)[colnames(df.pred.main) == "duration_median"] <- "duration"

# Subset to fewer samples, to lower computational cost for now
set.seed(155)
samps <- sample(df.pred.main$sample_num, 200, replace = TRUE)

# Subset to these samples
df.pred.main <- subset(df.pred.main, sample_num %in% samps)

# Calculate & save all possible pairwise differences
route.low.df <- get_pairwise_differences(df.pred.main, cofactor = "route")
saveRDS(route.low.df, "./outputs/predictions/pred-route-differences-lowdoses-200.RDS")

dose.low.df <- get_pairwise_differences(df.pred.main, cofactor = "dose")
saveRDS(dose.low.df, "./outputs/predictions/pred-dose-differences-lowdoses-200.RDS")

sex.low.df <- get_pairwise_differences(df.pred.main, cofactor = "sex")
saveRDS(sex.low.df, "./outputs/predictions/pred-sex-differences-lowdoses-200.RDS")

age.low.df <- get_pairwise_differences(df.pred.main, cofactor = "age")
saveRDS(age.low.df, "./outputs/predictions/pred-age-differences-lowdoses-200.RDS")

species.low.df <- get_pairwise_differences(df.pred.main, cofactor = "species")
saveRDS(species.low.df, "./outputs/predictions/pred-species-differences-lowdoses-200.RDS")



# Predictions for survival curves ----------------------------------------------

# Draws per panel
n_draws <- 500

# Generate survival curves for non-AE exposures (larger doses only)
df.estimates.not.ae <- get_params_for_survival_curves(fit, n_draws,
                                                      route_options = c(1, 2, 3, 5),
                                                      tissue_options = c(1:4, 6),
                                                      dose_options = c(4, 5, 6, 7),
                                                      assay_options = c(1, 4))

# Generate survival curves for AE exposures (lower doses only)
df.estimates.ae <- get_params_for_survival_curves(fit, n_draws,
                                                  route_options = 4,
                                                  tissue_options = c(1:4, 6),
                                                  dose_options = c(1, 2, 3, 4),
                                                  assay_options = c(1, 4))

# Combine them
df.estimates <- rbind(df.estimates.not.ae, df.estimates.ae)

# Save them
write.csv(df.estimates, 
          "./outputs/predictions/pred-for-survival-curves-500.csv",
          row.names = FALSE)


