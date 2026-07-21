# This file: - Extracts event times from the full database

# Load the data  ---------------------------------------------------------------

# Full database
dat.all <- read.csv("./data/database-clean.csv")

# Remove individuals without dose information, intracranial inoculations, & other tissues
dat.all <- subset(dat.all, !is.na(dose_nose) & 
                           route_idx != 6 &
                           location_grp %in% c("URT", "LRT", "GI"))


# Adjust CT values -------------------------------------------------------------
# CT values need to be rescaled for consistency with viral loads 
#      (larger numbers must be larger viral loads in the later processing)

# Remove the SD from the CT values
dat.all$value[str_detect(dat.all$value, " [+]")] <- substr(dat.all$value[str_detect(dat.all$value, " [+]")], 1, 5)

# Subtract the CT values from 40 (so larger values indicate larger titers)
dat.all$value[dat.all$unit_rep %in% c( "Ct value" ,"Ct value (+/- SD)") &
                !is.na(as.numeric(dat.all$value))] <- 
  40 - as.numeric(dat.all$value[dat.all$unit_rep %in% c( "Ct value" ,"Ct value (+/- SD)") &
                                  !is.na(as.numeric(dat.all$value))])

# Bixler has days that are "9 or 10", we will use 9.5 for these
dat.all$day_post_infection[dat.all$day_post_infection == "9 or 10"] <- 9.5


# Individuals with ID names ----------------------------------------------------

# Remove the individuals without ID names
dat.id <- subset(dat.all, indiv_rep %notin% c("No ID", "MM-N-14(-1/-2)", "MM-O-14(-1/-2)"))

# Extract the event times, note that default warnings on NAs introduced by coercion & 
#    returning -Inf are expected and accounted for
dat.id <- extract_times_with_ids(dat.id)

# Extract each individual dataframe
dat.id.times <- dat.id$dat.times # df with only the event times
dat.id.all <- dat.id$dat.all # df with all sampled times
dat.id.times$group_number <- NA # NA b/c only no ID df has to keep track of group #s

# Fix flagged individuals from one study manually (by manually checking data)
dat.id.times$first_lower_bound[dat.id.times$first_lower_bound == "FLAG"] <- 0
dat.id.times$first_upper_bound[dat.id.times$first_upper_bound == "FLAG"] <- 10
dat.id.times$peak_lower_bound[!is.na(dat.id.times$peak_lower_bound) & 
                                dat.id.times$peak_lower_bound == "FLAG"] <- 0
dat.id.times$peak_upper_bound[!is.na(dat.id.times$peak_upper_bound) & 
                                dat.id.times$peak_upper_bound == "FLAG"] <- 10
dat.id.times$last_lower_bound[!is.na(dat.id.times$last_lower_bound) & 
                                dat.id.times$last_lower_bound == "FLAG"] <- NA
dat.id.times$last_upper_bound[!is.na(dat.id.times$last_upper_bound) & 
                                dat.id.times$last_upper_bound == "FLAG"] <- NA


# Individuals without ID names -------------------------------------------------

dat.noid <- subset(dat.all, indiv_rep %in% c("No ID"))
dat.noid <- subset(dat.noid, !str_detect(unit_rep, "Peak")) # Remove samples only reported as a peak because we don't have a censored interval to bound them

# Find article and location_grp/location_subgrps with ID and No ID data, because
#   these individuals will be double counted on accident unless removed
for (article.ii in unique(dat.all$article)) {
  dat.art <- subset(dat.all, article == article.ii)
  
  if ("No ID" %in% unique(dat.art$indiv_rep) & 
      length(unique(dat.art$indiv_rep)) > 1) {
    
    cat("\n", article.ii, ":", unique(dat.art$indiv_rep),   "\n")
    
    for (loc.grp.ii in unique(dat.art$location_grp)) {
      dat.loc.grp <- subset(dat.art, location_grp == loc.grp.ii)
      
      for (loc.subgrp.ii in unique(dat.loc.grp$location_subgrp)) {
        dat.loc.subgrp <- subset(dat.loc.grp, location_subgrp == loc.subgrp.ii)
        
        for (assay.ii in unique(dat.loc.subgrp$rna_type)) {
          dat.assay <- subset(dat.loc.subgrp, rna_type == assay.ii)
          
          if ("No ID" %in% unique(dat.assay$indiv_rep) & 
              length(unique(dat.assay$indiv_rep)) > 1) {
            cat(article.ii, loc.grp.ii, loc.subgrp.ii, "has ID and No ID data.\n")
          }
        }
      }
    }
  }
}

# Remove these individuals as necessary
dat.noid <- subset(dat.noid, !(article == "Salguero et al. 2021" & indiv_rep == "No ID"))
dat.noid <- subset(dat.noid, !(article == "Woolsey et al. 2020" & indiv_rep == "No ID" &
                                 location_grp %in% c("URT", "LRT", "GI", "Oral")))

# Remove invasive nose samples from Berry et al. 2022 and Liu et al. 2022
dat.noid <- subset(dat.noid, !(article %in% c("Berry et al. 2022", "Liu et al. 2022") &
                               sample_type == "Invasive" & location_subgrp == "Nose"))

# Run the function to extract times
dat.noid.times <- extract_times_without_ids(dat.noid)

# Add indiv_sample column to match with dataframe for ID individuals
dat.noid.times$indiv_sample <- paste(dat.noid.times$indiv, 
                                     dat.noid.times$location_subgrp,
                                     dat.noid.times$assay_type,
                                     dat.noid.times$pcr_target_gene, 
                                     sep = "_")


# Combine ----------------------------------------------------------------------

dat.times <- rbind(dat.id.times, dat.noid.times)

# Convert bounds to numerics for ease of later use
dat.times$first_lower_bound <- as.numeric(dat.times$first_lower_bound)
dat.times$first_upper_bound <- as.numeric(dat.times$first_upper_bound)
dat.times$peak_lower_bound <- as.numeric(dat.times$peak_lower_bound)
dat.times$peak_observed_time <- as.numeric(dat.times$peak_observed_time)
dat.times$peak_upper_bound <- as.numeric(dat.times$peak_upper_bound)
dat.times$last_lower_bound <- as.numeric(dat.times$last_lower_bound)
dat.times$last_upper_bound <- as.numeric(dat.times$last_upper_bound)

# Fill in any event times (e.g., know that peak needs to happen between first and last)
dat.times <- fill_in_event_times(dat.times)
dat.id.all <- fill_in_event_times(dat.id.all)

# Add indicator for peak titer
dat.times$has_titer <- 1
dat.times$has_titer[is.na(dat.times$peak_observed_titer)] <- 0
dat.times$peak_observed_titer[is.na(dat.times$peak_observed_titer)] <- 9999


# Save -------------------------------------------------------------------------

write.csv(dat.times, file = "./data/df-event-times.csv", row.names = FALSE)
write.csv(dat.id.all, file = "./data/df-all-times.csv", row.names = FALSE)


