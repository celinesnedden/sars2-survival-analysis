# All custom functions for our analyses are below. 

`%notin%` <- Negate(`%in%`)

# Function for plotting
plot_data_histograms <- function(data, cofactor) {
  
  # Creates a histogram with ggplot that stratifies by the indicated cofactor
  # data: the set of data to plot
  # cofactor: the cofactor used to stratify the plot
  # The function will also print the sample sizes to the console for visual confirmation
  
  
  if (cofactor == "route"){
    data$new_route_name <- data$route_name
    data$cofactor_stratify <- data$new_route_name
    data$title <- "Exposure Route"
  }
  else if (cofactor == "route_full"){
    data$new_route_name <- data$inoc_route_grp
    data$new_route_name[data$new_route_name %in% c("IN, IB", "IN, IT", "IT, IN")] <- "IN + IT/IB"
    data$new_route_name[data$new_route_name %in% c("IN, OR")] <- "IN + OR"
    data$new_route_name[data$new_route_name %in% c("Multi-route (3+)")] <- "Multi-route" 
    data$new_route_name <- factor(data$new_route_name,
                                  levels = c("IC", "OC", "IN", "OR",
                                             "IN + OR", "IT", "IB", 
                                             "IN + IT/IB", "AE", "Multi-route", 
                                             "IG"))
    data$cofactor_stratify <- data$new_route_name
    data$title <- "Exposure Route"
  }
  else if (cofactor == "dose") {
    data$cofactor_stratify <- floor(log10(as.numeric(data$inoc_dose_total_pfu)))
    data$cofactor_stratify[is.na(data$cofactor_stratify)] <- "Unknown"
    data$title <- "Exposure Dose"
  }
  else if (cofactor == "sex") {
    data$cofactor_stratify[data$sex_name == "Female"] <- "Female"
    data$cofactor_stratify[data$sex_name == "Male"] <- "Male"
    data$cofactor_stratify[data$sex_name == "Unknown"] <- "Unknown"
    data$cofactor_stratify <- factor(data$cofactor_stratify, 
                                     levels = c("Male", "Female", "Unknown"))
    data$title <- "Sex"
  }
  else if (cofactor == "age") {
    data$cofactor_stratify[data$age_name == "Juvenile"] <- "Juvenile"
    data$cofactor_stratify[data$age_name == "Adult"] <- "Adult"
    data$cofactor_stratify[data$age_name == "Geriatric"] <- "Geriatric"
    data$cofactor_stratify[data$age_name == "Unknown"] <- "Unknown"
    data$cofactor_stratify <- factor(data$cofactor_stratify, 
                                     levels = c("Juvenile", "Adult", "Geriatric", "Unknown"))
    data$title <- "Age Class"
  }
  else if (cofactor == "species") {
    data$cofactor_stratify[data$species_name == "Rhesus macaque"] <- "Rhesus"
    data$cofactor_stratify[data$species_name == "Cynomolgus macaque"] <- "Cynomolgus"
    data$cofactor_stratify[data$species_name == "African green monkey"] <- "African green"
    data$cofactor_stratify <- factor(data$cofactor_stratify, 
                                     levels = c("Rhesus", "Cynomolgus", "African green"))
    data$title <- "Species"
  }
  else if (cofactor == "assay_full") {
    data$cofactor_stratify <- data$rna_type
    data$cofactor_stratify[data$cofactor_stratify == "totRNA"] <- "Total RNA"
    data$cofactor_stratify[data$cofactor_stratify == "sgRNA"] <- "sgRNA"
    data$cofactor_stratify[data$cofactor_stratify == "gRNA"] <- "gRNA"
    data$cofactor_stratify[data$cofactor_stratify == "culture"] <- "Culture"
    data$cofactor_stratify[data$cofactor_stratify == "Unknown"] <- "Unknown"
    
    data$cofactor_stratify <- factor(data$cofactor_stratify, 
                                     levels = c("Total RNA", "gRNA", 
                                                "sgRNA", "Culture", 
                                                "Unknown"))
    data$title <- "Assay"
  }
  else if (cofactor == "assay") {
    data$cofactor_stratify[str_detect(data$assay_type, "RNA|Unknown")] <- "PCR"
    data$cofactor_stratify[!str_detect(data$assay_type, "RNA|Unknown")] <- "Culture"
    data$cofactor_stratify <- factor(data$cofactor_stratify, 
                                     levels = c("PCR", "Culture"))
    data$title <- "Assay"
  }
  else if (cofactor == "location") {
    data$cofactor_stratify <- data$location_grp
    data$cofactor_stratify[data$cofactor_stratify %notin% c("URT", "LRT", "GI")] <- "Other"
    #data$cofactor_stratify[data$cofactor_stratify %in% c("URT", "LRT")] <- "Respiratory"
    #data$cofactor_stratify[data$cofactor_stratify %in% c("GI")] <- "Gastrointestinal"
    
    data$title <- "Location"
    data$cofactor_stratify <- factor(data$cofactor_stratify,
                                     levels = c("URT", "LRT", "GI", "Other"))
    #levels = c("Respiratory", "Gastrointestinal", "Other"))
  }
  
  data.counts <- count(data, cofactor_stratify)
  print(data.counts)
  
  fig <- ggplot(data, aes(x = cofactor_stratify, 
                          y = (after_stat(count))/sum(after_stat(count)))) + 
    geom_bar(color = "black", linewidth = 0.2, width = 0.8) + 
    geom_text(stat='count', aes(label=after_stat(count)),
              size = 2, vjust = -0.3) +
    #geom_text(data = data.counts,
    #          aes(x = cofactor_stratify, y = 0.01, label = n), angle = 90,
    #          size = 2.5, hjust = 0, color = "grey55") +
    #scale_x_discrete(breaks = breaks, labels = names) +
    scale_y_continuous(limits = c(0, 0.95), expand = c(0, 0),
                       labels = scales::percent) +
    facet_wrap(.~ title) +
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.9, "line"),
          axis.title = element_blank(),
          axis.text = element_text(size = 8),
          axis.text.x = element_text(size = 7, hjust = 1, angle =45),
          axis.ticks.length = unit(0, "cm"),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          strip.text = element_text(face = "bold", size = 9, vjust = -1),
          strip.background = element_blank(),
          plot.background = element_rect(fill = "transparent", colour = NA_character_),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.2, linetype = "solid"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major = element_line(size = 0.08, linetype = 'solid',
                                          colour = "grey88"),
          plot.margin = margin(0,3,0,0,"pt"))
  
  if (cofactor == "dose") {
    fig <- fig + 
      scale_x_discrete(labels = c(expression(paste(10^1)),
                                  expression(paste(10^2)), expression(paste(10^3)), 
                                  expression(paste(10^4)), expression(paste(10^5)), 
                                  expression(paste(10^6)), expression(paste(10^7)), "Unknown")) +
      labs(y = "Observations") +
      theme(axis.title.y = element_text())
  }
  else {
    fig <- fig +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank())
  }
  
  return(fig)
  
} 


# Write a function to generate single-tissue inoculation predictions
generate_single_tissue_inoculations <- function(sample_pars,
                                                sample_nums,
                                                tissue_inoculated,
                                                dose_total,
                                                assay_idx,
                                                age_idx,
                                                sex_idx,
                                                species_idx,
                                                lab_idx,
                                                lab_effect) {
  
  # Set up dataframe to store predictions
  df.pred <- data.frame(percent = NA, median_first = NA, median_peak = NA,
                        mean_titer = NA, median_last = NA, sample_first = NA,
                        sample_peak = NA, sample_titer = NA, sample_last = NA,
                        shape_first = NA, scale_first = NA, shape_peak = NA,
                        scale_peak = NA, shape_last = NA, scale_last = NA,
                        tissue_idx = NA, sample_num = NA,
                        age_idx = NA, sp_idx = NA, sex_idx = NA, assay_idx = NA)
  
  # Set up doses based on the inoculated tissue
  if (tissue_inoculated == "Nose"){dose_nose <- dose_total} else {dose_nose <- 0}
  if (tissue_inoculated == "Throat"){dose_throat <- dose_total} else {dose_throat <- 0}
  if (tissue_inoculated == "Trachea"){dose_trachea <- dose_total} else {dose_trachea <- 0}
  if (tissue_inoculated == "Lung"){dose_lung <- dose_total} else {dose_lung <- 0}
  if (tissue_inoculated == "Upper GI"){dose_gi <- dose_total} else {dose_gi <- 0}
  
  # Set up inoculation flags
  if (tissue_inoculated == "Nose"){inoc_nose_idx <- 1} else {inoc_nose_idx <- 0}
  if (tissue_inoculated == "Throat"){inoc_throat_idx <- 1} else {inoc_throat_idx <- 0}
  if (tissue_inoculated == "Trachea"){inoc_trachea_idx <- 1} else {inoc_trachea_idx <- 0}
  if (tissue_inoculated == "Lung"){inoc_lung_idx <- 1} else {inoc_lung_idx <- 0}
  if (tissue_inoculated == "Upper GI"){inoc_gi_idx <- 1} else {inoc_gi_idx <- 0}
  
  # Set up route numbers
  if (tissue_inoculated == "Nose"){route_idx <- 1}
  if (tissue_inoculated == "Throat"){route_idx <- 1} 
  if (tissue_inoculated == "Trachea"){route_idx <- 2}
  if (tissue_inoculated == "Lung"){route_idx <- 2}
  if (tissue_inoculated == "Upper GI"){route_idx <- 5}
  
  # Inoculate straight into the nose
  for (ii in 1:length(sample_nums)) {
    
    for (age in age_idx) {
      
      for (species in species_idx) {
        
        for (sex in sex_idx) {
          
          for (assay in assay_idx) {
            
            # Nose
            nose_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 1, location_idx = 0, tissue_idx = 1, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = inoc_nose_idx, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            nose_row$tissue_idx <- 1
            
            # Throat
            throat_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 1, location_idx = 1, tissue_idx = 2, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = inoc_throat_idx, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            throat_row$tissue_idx <- 2
            
            # Trachea
            trachea_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 2, location_idx = 0, tissue_idx = 3, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = inoc_trachea_idx, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            trachea_row$tissue_idx <- 3
            
            # Lung
            lung_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 2, location_idx = 1, tissue_idx = 4, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = inoc_lung_idx, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            lung_row$tissue_idx <- 4
            
            # Upper GI
            upgi_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 3, location_idx = 0, tissue_idx = 5, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = inoc_gi_idx, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            upgi_row$tissue_idx <- 5
            
            # Lower GI
            logi_row <- 
              get_individual_parameter_sample(sample.pars, ii = sample_nums[ii], 
                                              organ_idx = 3, location_idx = 1, tissue_idx = 6, 
                                              route_idx = route_idx, 
                                              dose_nose = dose_nose, dose_throat = dose_throat, 
                                              dose_trachea = dose_trachea, dose_lung = dose_lung, 
                                              dose_gi = dose_gi, 
                                              inoc_idx = 0, 
                                              assay_idx = assay, age_idx = age, 
                                              sex_idx = sex, sp_idx = species, 
                                              lab_idx = lab_idx, lab_effect = lab_effect)
            logi_row$tissue_idx <- 6
            
            
            # Combine them all into the existing dataframe, including the sample number
            next_rows <- rbind(nose_row, throat_row, trachea_row,
                               lung_row, upgi_row, logi_row)
            next_rows$sample_num <- sample_nums[ii]
            next_rows$age_idx <- age
            next_rows$sex_idx <- sex
            next_rows$sp_idx <- species
            next_rows$assay_idx <- assay
            
            df.pred <- rbind(df.pred, next_rows)
            
            
          }
        }
        
      }
      
    }
    
  }
  
  # Set column with the tissue inoculated
  df.pred$tissue_inoculated <- tissue_inoculated
  
  # Add covariates
  df.pred$dose_total <- dose_total
  #df.pred$assay_idx <- assay_idx
  #df.pred$age_idx <- age_idx
  #df.pred$sex_idx <- sex_idx
  #df.pred$sp_idx <- species_idx
  
  # Return it
  return(df.pred[-1, ])
  
}

# Function to create the adjacency matrices
generate_adjacency_matrix <- function(df, tissue_inoc) {
  
  # Empty adjacency matrix to fill 
  df.adj <- data.frame(Nose = rep(0, 6), 
                       Throat = rep(0, 6), 
                       Trachea = rep(0, 6), 
                       Lung = rep(0, 6), 
                       `Upper GI` = rep(0, 6), 
                       `Lower GI` = rep(0, 6))
  rownames(df.adj) <- c("Nose", "Throat", "Trachea", 
                        "Lung", "Upper GI", "Lower GI")
  
  # Subset to the correct set of samples from target inoculated tissues
  df <- subset(df, tissue_inoculated == tissue_inoc)
  
  # Subset to only the tissues that are predicted to test positive
  df <- subset(df, ever_positive == 1)
  
  # Set numeric for the inoculated tissue
  df$tissue_inoc_idx <- as.numeric(df$tissue_inoculated)
  
  # Set the sample & inoculated tissue to loop over
  df$inoc_sample <- paste0(df$tissue_inoculated, "-", df$sample_num, "-",
                           df$assay_idx, "-", df$age_idx, "-",
                           df$sex_idx, "-", df$sp_idx)
  
  # Loop over all the sampled tissues for each sample to fill in the 
  #    adjacency matrix
  for (samp in unique(df$inoc_sample)) {
    
    # Subset to the target set of tissues
    df.indiv <- subset(df, inoc_sample == samp)
    
    # Subset out tissues testing positive within the first day that aren't 
    #   the inoculated tissue
    df.indiv.0 <- subset(df.indiv, time_floor == 0 & tissue_idx != tissue_inoc_idx)
    tissue_inoc <- unique(df.indiv$tissue_inoc_idx)
    
    # Connect (i.e., increment) inoculated tissue to all other day 0 tissues
    if (nrow(df.indiv.0) >= 1) {
      for (ii in 1:nrow(df.indiv.0)){
        tissue_next <- df.indiv.0$tissue_idx[ii]
        df.adj[tissue_inoc, tissue_next] <- df.adj[tissue_inoc, tissue_next] + 1
      }
    }
    
    # Remove inoculated tissues for the rest of the checks, if there were other
    #   day 0 zero tissues 
    if (nrow(df.indiv.0) >= 1) {
      df.indiv <- subset(df.indiv, tissue_idx != tissue_inoc_idx)
    }
    
    # Get all subsequent days with positive tissues, ordered 
    dpis <- sort(unique(df.indiv$time_floor))
    
    # Connect all day X tissues with all tissues from the next positive day,
    #   by looping over the available DPIs
    if (length(dpis) > 1) {
      
      # Loop over all but the last DPI that becomes detectable
      for (jj in 1:(length(dpis)-1)) {
        
        df.indiv.this <- subset(df.indiv, time_floor == dpis[jj]) # this time point
        df.indiv.next <- subset(df.indiv, time_floor == dpis[jj + 1]) # next time point
        
        # Connect (i.e., increment) all tissues at this time point with all
        #   tissues that test positive at the next time point
        for (kk in 1:nrow(df.indiv.this)) {
          tissue_this <- df.indiv.this$tissue_idx[kk]
          
          for (ll in 1:nrow(df.indiv.next)) {
            tissue_next <- df.indiv.next$tissue_idx[ll]
            
            # If time point is flagged as greater than 20, then it's a residual
            #   never positive and will be excluded
            if (dpis[jj + 1] < 20) {
              df.adj[tissue_this, tissue_next] <- df.adj[tissue_this, tissue_next] + 1
            }
          }
        }
      }
    }
  }
  
  return(df.adj)
  
}


get_metrics_for_specific_labs <- function(fit, n_draws,
                                          seed = NA,
                                          dose_nose = NA,
                                          dose_throat = NA,
                                          dose_trachea = NA,
                                          dose_lung = NA,
                                          dose_gi = NA,
                                          dose_total = NA,
                                          tissue = 1,
                                          assay = 1,
                                          sex = 0,
                                          age = 2,
                                          species = 1,
                                          route = 3,
                                          lab = 1) {
  
  # Track time needed to generate predictions
  start_time <- Sys.time()

  # Get the parameter samples
  sample.pars <- get_all_parameter_samples(fit)
  
  # Set up data frame to store predictions
  pred.df  <- data.frame(sample_num = numeric(),
                         organ_group = character(),
                         organ_idx = numeric(),
                         tissue_idx = numeric(),
                         dose_nose = numeric(),
                         dose_throat = numeric(),
                         dose_trachea = numeric(),
                         dose_lung = numeric(),
                         dose_gi = numeric(),
                         dose_total = numeric(),
                         route_idx = numeric(),
                         sp_idx = numeric(),
                         age_idx = numeric(),
                         sex_idx = numeric(),
                         assay_idx = numeric(),
                         location_idx = numeric(),
                         lab_idx = numeric(),
                         percent_positive = numeric(),
                         first_pos_median = numeric(),
                         peak_median = numeric(),
                         titer_mean = numeric(),
                         last_median = numeric())
  
  # Get random sample numbers to loop over
  if (!is.na(seed)){
    set.seed(seed)
    sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  }
  else {
    sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  }
  
  # Set organ groups and location_idx
  tissue.ii <- tissue
  if (tissue.ii == 1) {
    organ_group <- "URT"
    organ.ii <- 1
    loc.ii <- 0
  }
  else if (tissue.ii == 2) {
    organ_group <- "URT"
    organ.ii <- 1
    loc.ii <- 1
  }
  else if (tissue.ii == 3) {
    organ_group <- "LRT"
    organ.ii <- 2
    loc.ii <- 0
  }
  else if (tissue.ii == 4) {
    organ_group <- "LRT"
    organ.ii <- 2
    loc.ii <- 1
  }
  else if (tissue.ii == 5) {
    organ_group <- "GI"
    organ.ii <- 3
    loc.ii <- 0
  }
  else if (tissue.ii == 6) {
    organ_group <- "GI"
    organ.ii <- 3
    loc.ii <- 1
  }
  
  # Set inoculated category
  inoc <- 0 # set zeros first, then override with 1s if necessary
  if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0){inoc <- 1}
  if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1){inoc <- 1}
  if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0){inoc <- 1}
  if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1){inoc <- 1}
  if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0){inoc <- 1}
  
  
  for (ii in sample_num) {

    sampled_metrics <- get_individual_parameter_sample(sample.pars,
                                                       ii,
                                                       organ.ii,
                                                       loc.ii,
                                                       tissue.ii,
                                                       route,
                                                       dose_nose,
                                                       dose_throat,
                                                       dose_trachea,
                                                       dose_lung,
                                                       dose_gi,
                                                       inoc,
                                                       assay,
                                                       age,
                                                       sex,
                                                       species,
                                                       lab,
                                                       lab_effect = "Yes")
    
    # Add estimates 
    new.obs <- data.frame(sample_num = ii,
                          organ_group = organ_group,
                          organ_idx = organ.ii,
                          tissue_idx = tissue.ii,
                          dose_nose = dose_nose,
                          dose_throat = dose_throat,
                          dose_trachea = dose_trachea,
                          dose_lung = dose_lung,
                          dose_gi = dose_gi,
                          dose_total = dose_total,
                          route_idx = route,
                          sp_idx = species,
                          age_idx = age,
                          sex_idx = sex,
                          assay_idx = assay,
                          location_idx = loc.ii,
                          lab_idx = lab,
                          percent_positive = sampled_metrics$percent,
                          first_pos_median = sampled_metrics$median_first,
                          peak_median = sampled_metrics$median_peak,
                          titer_mean = sampled_metrics$mean_titer,
                          last_median = sampled_metrics$median_last
    )
    
    pred.df <- rbind(pred.df, new.obs)
    
  }
  
  # Calculate AUC values
  pred.df <- calculate_auc(pred.df)

  # Print out computation time
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  cat("\n")
  print(paste("Total computation time:", round(elapsed_time, 2), "minutes"))
  
  # Return collection of predictions
  return(pred.df)
}


get_all_parameter_samples <- function(fit){
  
  pars <- c("percent_intercept", "percent_location", #"percent_inoc",
            "percent_route", "percent_dose",
            "percent_sex", "percent_age", "percent_species",
            "percent_assay", "percent_lab",
            
            "shape_intercept_first",
            
            "median_intercept_first", "median_location_first",
            "median_route_first", "median_dose_first", 
            "median_sex_first", "median_age_first", "median_species_first",
            "median_assay_first", "median_lab_first", 
            
            "shape_intercept_peak",
            
            "median_intercept_peak", "median_location_peak",
            "median_route_peak", "median_dose_peak", 
            "median_sex_peak", "median_age_peak", "median_species_peak",
            "median_assay_peak", "median_lab_peak", 
            "median_time_peak",
            
            "titer_intercept", "titer_location", #"titer_inoc",
            "titer_route", "titer_dose",
            "titer_sex", "titer_age", "titer_species",
            "titer_assay",  
            "titer_time",
            "titer_true_sd", 
            
            "titer_lab",
            
            "shape_intercept_last",
            
            "median_intercept_last", "median_location_last",
            "median_route_last", "median_dose_last", 
            "median_sex_last", "median_age_last", "median_species_last",
            "median_assay_last", "median_lab_last", 
            "median_time_last"
  )
  
  # Extract in wide format, to get correlated samples
  sample.pars <- fit$draws(pars, format = "df")
  
  return(sample.pars)
  
}

get_individual_parameter_sample <- function(sample.pars,
                                            ii, # row to get value from
                                            organ_idx,
                                            location_idx,
                                            tissue_idx,
                                            route_idx,
                                            dose_nose,
                                            dose_throat,
                                            dose_trachea,
                                            dose_lung,
                                            dose_gi,
                                            inoc_idx,
                                            assay_idx,
                                            age_idx,
                                            sex_idx,
                                            sp_idx,
                                            lab_idx,
                                            lab_effect) {
  
  # sample.pars: a wide dataframe of parameter samples, obtained from
  #              the function get_all_parameter_samples
  
  # PERCENT POSITIVE SAMPLES
  
  percent_intercept_column <- paste0("percent_intercept[", organ_idx, "]") 
  percent_intercept <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_intercept_column)])
  
  percent_location_column <- paste0("percent_location[", organ_idx, "]") 
  percent_location <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_location_column)])
  
  percent_route_column <- paste0("percent_route[", route_idx, ",", tissue_idx, "]") 
  percent_route <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_route_column)])
  
  percent_dose_nose_column <- paste0("percent_dose[1,", tissue_idx, "]") 
  percent_dose_nose <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_dose_nose_column)])
  
  percent_dose_throat_column <- paste0("percent_dose[2,", tissue_idx, "]") 
  percent_dose_throat <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_dose_throat_column)])
  
  percent_dose_trachea_column <- paste0("percent_dose[3,", tissue_idx, "]") 
  percent_dose_trachea <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_dose_trachea_column)])
  
  percent_dose_lung_column <- paste0("percent_dose[4,", tissue_idx, "]") 
  percent_dose_lung <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_dose_lung_column)])
  
  percent_dose_gi_column <- paste0("percent_dose[5,", tissue_idx, "]") 
  percent_dose_gi <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_dose_gi_column)])
  
  percent_sex_column <- paste0("percent_sex[", organ_idx, "]") 
  percent_sex <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_sex_column)])
  
  percent_age_column <- paste0("percent_age[", age_idx, ",", organ_idx, "]")
  percent_age <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_age_column)])
  
  percent_species_column <- paste0("percent_species[", sp_idx, ",", organ_idx, "]")
  percent_species <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_species_column)])
  
  percent_assay_column <- paste0("percent_assay[", assay_idx, ",", organ_idx, "]")
  percent_assay <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_assay_column)])
  
  percent_lab_column <- paste0("percent_lab[", lab_idx, ",", organ_idx,  "]")
  percent_lab <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == percent_lab_column)])
  
  if (lab_effect == "No") {
    percent_lab <- 0
  }
  
  percent_trans <- percent_intercept + 
    percent_location * location_idx +
    percent_dose_nose * dose_nose +
    percent_dose_throat * dose_throat +
    percent_dose_trachea * dose_trachea +
    percent_dose_lung * dose_lung +
    percent_dose_gi * dose_gi +
    percent_route + 
    percent_age +
    percent_species +
    percent_sex * sex_idx +
    percent_assay + 
    percent_lab
  
  percent <- exp(percent_trans) / (1 + exp(percent_trans))
  
  
  # FIRST POSITIVIE
  
  ## Get shape parameter 
  
  shape_first_column <- paste0("shape_intercept_first[", organ_idx, "]")
  shape_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == shape_first_column)])
  
  
  ## Get scale parameter 
  
  median_intercept_column <- paste0("median_intercept_first[", organ_idx, "]")
  median_intercept_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
  
  median_location_column <-paste0("median_location_first[", organ_idx, "]") 
  median_location_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
  
  median_route_column <- paste0("median_route_first[", route_idx, ",", tissue_idx, "]") 
  median_route_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
  
  median_dose_nose_column <- paste0("median_dose_first[1,", tissue_idx, "]") 
  median_dose_nose_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
  
  median_dose_throat_column <- paste0("median_dose_first[2,", tissue_idx, "]") 
  median_dose_throat_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
  
  median_dose_trachea_column <- paste0("median_dose_first[3,", tissue_idx, "]") 
  median_dose_trachea_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
  
  median_dose_lung_column <- paste0("median_dose_first[4,", tissue_idx, "]") 
  median_dose_lung_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
  
  median_dose_gi_column <- paste0("median_dose_first[5,", tissue_idx, "]") 
  median_dose_gi_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
  
  median_sex_column <- paste0("median_sex_first[", organ_idx, "]") 
  median_sex_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
  
  median_age_column <- paste0("median_age_first[", age_idx, ",", organ_idx, "]")
  median_age_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
  
  median_species_column <- paste0("median_species_first[", sp_idx, ",", organ_idx,  "]")
  median_species_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
  
  median_assay_column <- paste0("median_assay_first[", assay_idx, ",", organ_idx, "]")
  median_assay_first <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
  
  median_lab_column <- paste0("median_lab_first[", lab_idx, ",", organ_idx,  "]")
  median_lab_first <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
  
  if (lab_effect == "No") {
    median_lab_first <- 0
  }
  
  # Calculate sample median
  
  median_first <- exp(-1/10 * (median_intercept_first + 
                                 median_location_first * location_idx +  
                                 median_route_first +  
                                 median_dose_nose_first * dose_nose +  
                                 median_dose_throat_first * dose_throat + 
                                 median_dose_trachea_first * dose_trachea + 
                                 median_dose_lung_first * dose_lung + 
                                 median_dose_gi_first * dose_gi + 
                                 median_species_first + 
                                 median_age_first +  
                                 median_sex_first * sex_idx +  
                                 median_assay_first + 
                                 median_lab_first
  ))
  
  scale_first <- median_first / (log(2)^(1 / shape_first))
  first_pos_sample <- rweibull(1, shape = shape_first, scale = scale_first)
  
  
  # PEAK TIME
  
  ## Get shape parameter 
  
  shape_peak_column <- paste0("shape_intercept_peak[", organ_idx, "]")
  shape_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == shape_peak_column)])
  
  ## Get scale parameter 
  
  median_intercept_column <- paste0("median_intercept_peak[", organ_idx, "]")
  median_intercept_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
  
  median_location_column <-paste0("median_location_peak[", organ_idx, "]") 
  median_location_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
  
  median_route_column <- paste0("median_route_peak[", route_idx, ",", tissue_idx, "]") 
  median_route_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
  
  median_dose_nose_column <- paste0("median_dose_peak[1,", tissue_idx, "]") 
  median_dose_nose_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
  
  median_dose_throat_column <- paste0("median_dose_peak[2,", tissue_idx, "]") 
  median_dose_throat_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
  
  median_dose_trachea_column <- paste0("median_dose_peak[3,", tissue_idx, "]") 
  median_dose_trachea_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
  
  median_dose_lung_column <- paste0("median_dose_peak[4,", tissue_idx, "]") 
  median_dose_lung_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
  
  median_dose_gi_column <- paste0("median_dose_peak[5,", tissue_idx, "]") 
  median_dose_gi_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
  
  median_sex_column <- paste0("median_sex_peak[", organ_idx, "]") 
  median_sex_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
  
  median_age_column <- paste0("median_age_peak[", age_idx, ",", organ_idx, "]")
  median_age_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
  
  median_species_column <- paste0("median_species_peak[", sp_idx, ",", organ_idx,  "]")
  median_species_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
  
  median_assay_column <- paste0("median_assay_peak[", assay_idx, ",", organ_idx, "]")
  median_assay_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
  
  median_lab_column <- paste0("median_lab_peak[", lab_idx, ",", organ_idx,  "]")
  median_lab_peak <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
  
  median_time_column <- paste0("median_time_peak[", inoc_idx + 1, ",", organ_idx, "]")
  median_time_peak <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
  
  if (lab_effect == "No") {
    median_lab_peak <- 0
  }
  
  # Calculate sample median
  median_peak <- exp(-1/10 * (median_intercept_peak + 
                                median_location_peak * location_idx +  
                                median_route_peak +  
                                median_dose_nose_peak * dose_nose +  
                                median_dose_throat_peak * dose_throat + 
                                median_dose_trachea_peak * dose_trachea + 
                                median_dose_lung_peak * dose_lung + 
                                median_dose_gi_peak * dose_gi + 
                                median_species_peak + 
                                median_age_peak +  
                                median_sex_peak * sex_idx +  
                                median_assay_peak + 
                                median_lab_peak +
                                median_time_peak * median_first
  ))
  
  scale_peak <- median_peak / (log(2)^(1 / shape_peak))
  
  # Draw an individual sample 
  median_peak_samp <- exp(-1/10 * (median_intercept_peak + 
                                     median_location_peak * location_idx +  
                                     median_route_peak +  
                                     median_dose_nose_peak * dose_nose +  
                                     median_dose_throat_peak * dose_throat + 
                                     median_dose_trachea_peak * dose_trachea + 
                                     median_dose_lung_peak * dose_lung + 
                                     median_dose_gi_peak * dose_gi + 
                                     median_species_peak + 
                                     median_age_peak +  
                                     median_sex_peak * sex_idx +  
                                     median_assay_peak + 
                                     median_lab_peak +
                                     median_time_peak * first_pos_sample))
  
  scale_peak_samp <- median_peak_samp / (log(2)^(1 / shape_peak))
  peak_sample <- rweibull(1, shape = shape_peak, scale = scale_peak_samp)
  
  
  # PEAK TITER
  
  ## Get parameters 
  
  titer_intercept_column <-paste0("titer_intercept[", organ_idx, "]") 
  titer_intercept <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_intercept_column)])
  
  titer_location_column <-paste0("titer_location[", organ_idx, "]") 
  titer_location <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_location_column)])
  
  #titer_inoc_column <-paste0("titer_inoc[", organ_idx, "]") 
  #titer_inoc <- as.numeric(
  #  sample.pars[ii, which(colnames(sample.pars) == titer_inoc_column)])
  
  titer_route_column <- paste0("titer_route[", route_idx, ",", tissue_idx, "]") 
  titer_route <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_route_column)])
  
  titer_dose_nose_column <- paste0("titer_dose[1,", tissue_idx, "]") 
  titer_dose_nose <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_dose_nose_column)])
  
  titer_dose_throat_column <- paste0("titer_dose[2,", tissue_idx, "]") 
  titer_dose_throat <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_dose_throat_column)])
  
  titer_dose_trachea_column <- paste0("titer_dose[3,", tissue_idx, "]") 
  titer_dose_trachea <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_dose_trachea_column)])
  
  titer_dose_lung_column <- paste0("titer_dose[4,", tissue_idx, "]") 
  titer_dose_lung <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_dose_lung_column)])
  
  titer_dose_gi_column <- paste0("titer_dose[5,", tissue_idx, "]") 
  titer_dose_gi <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_dose_gi_column)])
  
  titer_sex_column <- paste0("titer_sex[", organ_idx, "]")
  titer_sex <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_sex_column)])
  
  titer_age_column <- paste0("titer_age[", age_idx, ",", organ_idx, "]")
  titer_age <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_age_column)])
  
  titer_species_column <- paste0("titer_species[", sp_idx, ",", organ_idx, "]")
  titer_species <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_species_column)])
  
  titer_assay_column <- paste0("titer_assay[", assay_idx, ",", organ_idx, "]")
  titer_assay <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_assay_column)])
  
  if (lab_idx <= 140) {
    titer_lab_column <- paste0("titer_lab[", lab_idx, ",", organ_idx,  "]")
    titer_lab <- as.numeric(
      sample.pars[ii, which(colnames(sample.pars) == titer_lab_column)])
  }
  else {titer_lab <- 0}

  titer_time_column <- paste0("titer_time[", inoc_idx + 1, ",", organ_idx, "]")
  titer_time <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_time_column)])
  
  titer_sd_column <- paste0("titer_true_sd[", organ_idx, "]")
  titer_sd <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == titer_sd_column)])
  
  if (lab_effect == "No") {
    titer_lab <- 0
  }
  
  #print(assay_idx)
  #print(organ_idx)
  # Mean titer
  
  titer_mean <- titer_intercept +
    titer_location * location_idx + 
    #titer_inoc * inoc_idx + 
    titer_route +
    titer_dose_nose * dose_nose +
    titer_dose_throat * dose_throat +
    titer_dose_trachea * dose_trachea +
    titer_dose_lung * dose_lung +
    titer_dose_gi * dose_gi +
    titer_sex * sex_idx +
    titer_age +
    titer_species +
    titer_assay +
    titer_lab +
    titer_time * median_peak # median peak
  
  titer_sample <- rnorm(1, titer_mean, titer_sd) 
  
  
  if (lab_idx > 140) {
    titer_mean <- NA
    titer_sample <- NA}
  
  
  # LAST POSITIVE
  
  ## Get shape parameter
  
  shape_last_column <- paste0("shape_intercept_last[", organ_idx, "]")
  shape_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == shape_last_column)])
  
  ## Get scale parameter 
  
  median_intercept_column <- paste0("median_intercept_last[", organ_idx, "]")
  median_intercept_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
  
  median_location_column <-paste0("median_location_last[", organ_idx, "]") 
  median_location_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
  
  median_route_column <- paste0("median_route_last[", route_idx, ",", tissue_idx, "]") 
  median_route_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
  
  median_dose_nose_column <- paste0("median_dose_last[1,", tissue_idx, "]") 
  median_dose_nose_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
  
  median_dose_throat_column <- paste0("median_dose_last[2,", tissue_idx, "]") 
  median_dose_throat_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
  
  median_dose_trachea_column <- paste0("median_dose_last[3,", tissue_idx, "]") 
  median_dose_trachea_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
  
  median_dose_lung_column <- paste0("median_dose_last[4,", tissue_idx, "]") 
  median_dose_lung_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
  
  median_dose_gi_column <- paste0("median_dose_last[5,", tissue_idx, "]") 
  median_dose_gi_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
  
  median_sex_column <- paste0("median_sex_last[", organ_idx, "]") 
  median_sex_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
  
  median_age_column <- paste0("median_age_last[", age_idx, ",", organ_idx, "]")
  median_age_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
  
  median_species_column <- paste0("median_species_last[", sp_idx, ",", organ_idx,  "]")
  median_species_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
  
  median_assay_column <- paste0("median_assay_last[", assay_idx, ",", organ_idx, "]")
  median_assay_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
  
  median_lab_column <- paste0("median_lab_last[", lab_idx, ",", organ_idx,  "]")
  median_lab_last <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
  
  median_time_column <- paste0("median_time_last[", inoc_idx + 1, ",", organ_idx, "]")
  median_time_last <- as.numeric(
    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
  
  if (lab_effect == "No") {
    median_lab_last <- 0
  }
  
  median_last <- exp(-1/10 * (median_intercept_last + 
                                median_location_last * location_idx +  
                                median_route_last +  
                                median_dose_nose_last * dose_nose +  
                                median_dose_throat_last * dose_throat + 
                                median_dose_trachea_last * dose_trachea + 
                                median_dose_lung_last * dose_lung + 
                                median_dose_gi_last * dose_gi + 
                                median_species_last + 
                                median_age_last +  
                                median_sex_last * sex_idx +  
                                median_assay_last + 
                                median_lab_last +
                                median_time_last * titer_mean #median_peak
  ))
  
  scale_last <- median_last / (log(2)^(1 / shape_last))
  
  median_last_sample <- exp(-1/10 * (median_intercept_last + 
                                     median_location_last * location_idx +  
                                     median_route_last +  
                                     median_dose_nose_last * dose_nose +  
                                     median_dose_throat_last * dose_throat + 
                                     median_dose_trachea_last * dose_trachea + 
                                     median_dose_lung_last * dose_lung + 
                                     median_dose_gi_last * dose_gi + 
                                     median_species_last + 
                                     median_age_last +  
                                     median_sex_last * sex_idx +  
                                     median_assay_last + 
                                     median_lab_last +
                                     median_time_last * titer_sample #median_peak
  ))
  
  scale_last_sample <- median_last_sample / (log(2)^(1 / shape_last))
  last_sample <- rweibull(1, shape = shape_last, scale = scale_last_sample)
  

  #print(percent)
  #print(median_first)
  #print(median_peak)
  #print(mean_titer)
  #print(median_last)
  #print(sample_first)
  #print(sample_peak)
  #print(sample_titer)
  #print(sample_last)
  #print(shape_first) 
  #print(scale_first) 
  #print(shape_peak) 
  #print(scale_peak) 
  #print(shape_last) 
  #print(scale_last) 
  
  # Combine into one output
  output <- data.frame(percent = percent,
                       median_first = median_first,
                       median_peak = median_peak,
                       mean_titer = titer_mean,
                       median_last = median_last,
                       sample_first = first_pos_sample,
                       sample_peak = peak_sample,
                       sample_titer = titer_sample,
                       sample_last = last_sample,
                       shape_first = shape_first,
                       scale_first = scale_first,
                       shape_peak = shape_peak,
                       scale_peak = scale_peak,
                       shape_last = shape_last,
                       scale_last = scale_last)
  #print(output)
  
  return(output)
  
}

calculate_auc <- function(df){
  
  #df: data frame with model predictions
  #     must have columns for median first time, mean peak titer, and median last time
  
  #df$auc <- 1/2 * (df$last_median - df$first_pos_median) * df$titer_mean # This is wrong... since they're delay times 
  
  df$auc <- 1/2 * (df$peak_median + df$last_median) * df$titer_mean
  return(df)
  
}



get_metrics_across_cofactors <- function(fit, n_draws,
                                         dose_options = c(4, 7),
                                         route_options = 1:5,
                                         age_options = 1:3,
                                         sex_options = 0:1,
                                         species_options = 1:3,
                                         tissue_options = c(1, 4, 6),
                                         assay_options = c(1, 4),
                                         lab_effect = "No",
                                         lab_options = 1,
                                         rescale_doses = FALSE) {
  
  # fit: a stan model fit to take parameter samples from
  # n_draws: the number of samples to take for each cofactor combination
  
  # Get the parameter samples
  sample.pars <- get_all_parameter_samples(fit)
  
  # Set up data frame to store predictions
  predictions <- data.frame(sample_num = numeric(),
                            organ_group = character(),
                            organ_idx = numeric(),
                            tissue_idx = numeric(),
                            dose_nose = numeric(),
                            dose_throat = numeric(),
                            dose_trachea = numeric(),
                            dose_lung = numeric(),
                            dose_gi = numeric(),
                            dose_total = numeric(),
                            route_idx = numeric(),
                            sp_idx = numeric(),
                            age_idx = numeric(),
                            sex_idx = numeric(),
                            assay_idx = numeric(),
                            location_idx = numeric(),
                            lab_idx = numeric(),
                            percent_positive = numeric(),
                            first_pos_median = numeric(),
                            peak_median = numeric(),
                            titer_mean = numeric(),
                            last_median = numeric())
  
  # Get random sample numbers to loop over
  sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  
  for (total.dose.ii in dose_options) {
    
    total_dose <- total.dose.ii
    assign_dose_distribution(total_dose)
    
    cat("\nGenerating estimates for total dose", total.dose.ii, "\n")
    
    for (route.ii in route_options) {
      
      cat("Generating estimates for route", route.ii, "\n")
      
      dose_nose <- get(paste0("route", route.ii, ".dose.nose"))
      dose_throat <- get(paste0("route", route.ii, ".dose.throat"))
      dose_trachea <- get(paste0("route", route.ii, ".dose.trachea"))
      dose_lung <- get(paste0("route", route.ii, ".dose.lung"))
      dose_gi <- get(paste0("route", route.ii, ".dose.gi"))
      
      if (rescale_doses == TRUE) {
        dose_nose <- dose_nose / get("max_dose_nose")
        dose_throat <- dose_throat / get("max_dose_throat")
        dose_trachea <- dose_trachea / get("max_dose_trachea")
        dose_lung <- dose_lung / get("max_dose_lung")
        dose_gi <- dose_gi / get("max_dose_gi")
        
        print(dose_nose)
        
      }
      
      for (tissue.ii in tissue_options) {
        
        cat("Generating estimates for tissue", tissue.ii, "\n")
        
        # Set organ groups and location_idx
        if (tissue.ii == 1) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 0
        }
        else if (tissue.ii == 2) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 1
        }
        else if (tissue.ii == 3) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 0
        }
        else if (tissue.ii == 4) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 1
        }
        else if (tissue.ii == 5) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 0
        }
        else if (tissue.ii == 6) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 1
        }
        
        # Set inoculated category
        inoc <- 0 # set zeros first, then override with 1s if necessary
        if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0){inoc <- 1}
        if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1){inoc <- 1}
        if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0){inoc <- 1}
        if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1){inoc <- 1}
        if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0){inoc <- 1}
        
        for (assay.ii in assay_options) {
          
          for (sex.ii in sex_options) {
            
            for (age.ii in age_options) {
              
              for (sp.ii in species_options) {
                
                for (lab.ii in lab_options) {
                  
                  for (ii in sample_num) {
                    
                    sampled_metrics <- get_individual_parameter_sample(sample.pars,
                                                                       ii,
                                                                       organ.ii,
                                                                       loc.ii,
                                                                       tissue.ii,
                                                                       route.ii,
                                                                       dose_nose,
                                                                       dose_throat,
                                                                       dose_trachea,
                                                                       dose_lung,
                                                                       dose_gi,
                                                                       inoc,
                                                                       assay.ii,
                                                                       age.ii,
                                                                       sex.ii,
                                                                       sp.ii,
                                                                       lab.ii,
                                                                       lab_effect)
                    
                    #print(sampled_metrics)
                    
                    # Add estimates
                    
                    new.obs <- data.frame(sample_num = ii,
                                          organ_group = organ_group,
                                          organ_idx = organ.ii,
                                          tissue_idx = tissue.ii,
                                          dose_nose = dose_nose,
                                          dose_throat = dose_throat,
                                          dose_trachea = dose_trachea,
                                          dose_lung = dose_lung,
                                          dose_gi = dose_gi,
                                          dose_total = total_dose,
                                          route_idx = route.ii,
                                          sp_idx = sp.ii,
                                          age_idx = age.ii,
                                          sex_idx = sex.ii,
                                          assay_idx = assay.ii,
                                          location_idx = loc.ii,
                                          lab_idx = lab.ii,
                                          percent_positive = sampled_metrics$percent,
                                          first_pos_median = sampled_metrics$median_first,
                                          peak_median = sampled_metrics$median_peak,
                                          titer_mean = sampled_metrics$mean_titer,
                                          last_median = sampled_metrics$median_last
                    )
                    
                    predictions <- rbind(predictions, new.obs)
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Calculate AUC values
  predictions <- calculate_auc(predictions)
  
  if (lab_effect == "No"){
    predictions$lab_idx <- NA
  }
  
  # Return collection of predictions
  return(predictions)
}



get_params_for_survival_curves <- function(fit, n_draws,
                                           dose_options = c(4, 7),
                                           route_options = 1:5,
                                           age_options = 1:3,
                                           sex_options = 0:1,
                                           species_options = 1:3,
                                           tissue_options = c(1, 4, 6),
                                           assay_options = c(1, 4),
                                           lab_effect = "No",
                                           lab_options = 1) {
  
  # Get parameter samples
  sample.pars <- get_all_parameter_samples(fit)
  sample_num <- sample(seq_len(nrow(sample.pars)), n_draws, replace = TRUE)
  
  # Create full design grid (replaces nested loops)
  design_grid <- expand.grid(
    total_dose = dose_options,
    route.ii = route_options,
    tissue.ii = tissue_options,
    assay.ii = assay_options,
    sex.ii = sex_options,
    age.ii = age_options,
    sp.ii = species_options,
    lab.ii = lab_options,
    KEEP.OUT.ATTRS = FALSE
  )
  
  n_combos <- nrow(design_grid)
  total_rows <- n_combos * n_draws
  
  out_list <- vector("list", total_rows)
  counter <- 1
  
  for (row in seq_len(n_combos)) {
    
    combo <- design_grid[row, ]
    
    total_dose <- combo$total_dose
    route.ii  <- combo$route.ii
    tissue.ii <- combo$tissue.ii
    assay.ii  <- combo$assay.ii
    sex.ii    <- combo$sex.ii
    age.ii    <- combo$age.ii
    sp.ii     <- combo$sp.ii
    lab.ii    <- combo$lab.ii
    
    # Assign dose distribution once per combo
    assign_dose_distribution(total_dose)
    
    dose_nose    <- get(paste0("route", route.ii, ".dose.nose"))
    dose_throat  <- get(paste0("route", route.ii, ".dose.throat"))
    dose_trachea <- get(paste0("route", route.ii, ".dose.trachea"))
    dose_lung    <- get(paste0("route", route.ii, ".dose.lung"))
    dose_gi      <- get(paste0("route", route.ii, ".dose.gi"))
    
    # Tissue mapping (vectorized replacement for if/else ladder)
    organ_lookup <- list(
      `1` = list(group="URT", organ=1, loc=0),
      `2` = list(group="URT", organ=1, loc=1),
      `3` = list(group="LRT", organ=2, loc=0),
      `4` = list(group="LRT", organ=2, loc=1),
      `5` = list(group="GI",  organ=3, loc=0),
      `6` = list(group="GI",  organ=3, loc=1)
    )
    
    organ_group <- organ_lookup[[as.character(tissue.ii)]]$group
    organ.ii    <- organ_lookup[[as.character(tissue.ii)]]$organ
    loc.ii      <- organ_lookup[[as.character(tissue.ii)]]$loc
    
    # Inoculation flag
    inoc <- 0
    if (dose_nose > 0    & organ_group == "URT" & loc.ii == 0) inoc <- 1
    if (dose_throat > 0  & organ_group == "URT" & loc.ii == 1) inoc <- 1
    if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0) inoc <- 1
    if (dose_lung > 0    & organ_group == "LRT" & loc.ii == 1) inoc <- 1
    if (dose_gi > 0      & organ_group == "GI"  & loc.ii == 0) inoc <- 1
    
    # Precompute scaled doses
    dose_nose_s    <- dose_nose    / 6.821186
    dose_throat_s  <- dose_throat  / 6.90768
    dose_trachea_s <- dose_trachea / 7.350608
    dose_lung_s    <- dose_lung    / 6.562293
    dose_gi_s      <- dose_gi      / 7
    
    # Loop only over sampled posterior draws
    for (ii in sample_num) {
      
      sampled_metrics <- get_individual_parameter_sample(
        sample.pars,
        ii,
        organ.ii,
        loc.ii,
        tissue.ii,
        route.ii,
        dose_nose_s,
        dose_throat_s,
        dose_trachea_s,
        dose_lung_s,
        dose_gi_s,
        inoc,
        assay.ii,
        age.ii,
        sex.ii,
        sp.ii,
        lab.ii,
        lab_effect
      )
      
      out_list[[counter]] <- data.frame(
        sample_num = ii,
        organ_group = organ_group,
        organ_idx = organ.ii,
        tissue_idx = tissue.ii,
        dose_nose = dose_nose,
        dose_throat = dose_throat,
        dose_trachea = dose_trachea,
        dose_lung = dose_lung,
        dose_gi = dose_gi,
        dose_total = total_dose,
        route_idx = route.ii,
        sp_idx = sp.ii,
        age_idx = age.ii,
        sex_idx = sex.ii,
        assay_idx = assay.ii,
        location_idx = loc.ii,
        lab_idx = lab.ii,
        percent_positive = sampled_metrics$percent,
        first_shape = sampled_metrics$shape_first,
        first_scale = sampled_metrics$scale_first,
        peak_shape = sampled_metrics$shape_peak,
        peak_scale = sampled_metrics$scale_peak,
        last_shape = sampled_metrics$shape_last,
        last_scale = sampled_metrics$scale_last
      )
      
      counter <- counter + 1
    }
  }
  
  predictions <- data.table::rbindlist(out_list)
  
  if (lab_effect == "No") {
    predictions$lab_idx <- NA
  }
  
  return(predictions)
}


get_individual_times <- function(fit, n_draws,
                                 dose_options = c(4, 7),
                                 route_options = 1:5,
                                 age_options = 1:3,
                                 sex_options = 0:1,
                                 species_options = 1:3,
                                 tissue_options = c(1, 4, 6),
                                 assay_options = c(1, 4),
                                 lab_effect = "No",
                                 lab_options = 1,
                                 rescale_doses = TRUE) {
  
  # fit: a stan model fit to take parameter samples from
  # n_draws: the number of samples to take for each cofactor combination
  
  # Get the parameter samples
  sample.pars <- get_all_parameter_samples(fit)
  
  # Set up data frame to store predictions
  predictions <- data.frame(sample_num = numeric(),
                            organ_group = character(),
                            organ_idx = numeric(),
                            tissue_idx = numeric(),
                            dose_nose = numeric(),
                            dose_throat = numeric(),
                            dose_trachea = numeric(),
                            dose_lung = numeric(),
                            dose_gi = numeric(),
                            dose_total = numeric(),
                            route_idx = numeric(),
                            sp_idx = numeric(),
                            age_idx = numeric(),
                            sex_idx = numeric(),
                            assay_idx = numeric(),
                            location_idx = numeric(),
                            lab_idx = numeric(),
                            percent_positive = numeric(),
                            first_sample = numeric(),
                            titer_sample = numeric(),
                            peak_sample = numeric(),
                            last_sample = numeric())
  
  # Get random sample numbers to loop over
  sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  
  for (total.dose.ii in dose_options) {
    
    total_dose <- total.dose.ii
    assign_dose_distribution(total_dose)
    
    cat("\nGenerating estimates for total dose", total.dose.ii, "\n")
    
    for (route.ii in route_options) {
      
      cat("Generating estimates for route", route.ii, "\n")
      
      dose_nose <- get(paste0("route", route.ii, ".dose.nose"))
      dose_throat <- get(paste0("route", route.ii, ".dose.throat"))
      dose_trachea <- get(paste0("route", route.ii, ".dose.trachea"))
      dose_lung <- get(paste0("route", route.ii, ".dose.lung"))
      dose_gi <- get(paste0("route", route.ii, ".dose.gi"))
      
      if (rescale_doses == TRUE) {
        dose_nose <- dose_nose / get("max_dose_nose")
        dose_throat <- dose_throat / get("max_dose_throat")
        dose_trachea <- dose_trachea / get("max_dose_trachea")
        dose_lung <- dose_lung / get("max_dose_lung")
        dose_gi <- dose_gi / get("max_dose_gi")
        
        print(dose_nose)
        
      }
      
      for (tissue.ii in tissue_options) {
        
        cat("Generating estimates for tissue", tissue.ii, "\n")
        
        # Set organ groups and location_idx
        if (tissue.ii == 1) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 0
        }
        else if (tissue.ii == 2) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 1
        }
        else if (tissue.ii == 3) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 0
        }
        else if (tissue.ii == 4) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 1
        }
        else if (tissue.ii == 5) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 0
        }
        else if (tissue.ii == 6) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 1
        }
        
        # Set inoculated category
        inoc <- 0 # set zeros first, then override with 1s if necessary
        if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0){inoc <- 1}
        if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1){inoc <- 1}
        if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0){inoc <- 1}
        if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1){inoc <- 1}
        if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0){inoc <- 1}
        
        for (assay.ii in assay_options) {
          
          for (sex.ii in sex_options) {
            
            for (age.ii in age_options) {
              
              for (sp.ii in species_options) {
                
                for (lab.ii in lab_options) {
                  
                  for (ii in sample_num) {
                    
                    sampled_metrics <- get_individual_parameter_sample(sample.pars,
                                                                       ii,
                                                                       organ.ii,
                                                                       loc.ii,
                                                                       tissue.ii,
                                                                       route.ii,
                                                                       dose_nose,
                                                                       dose_throat,
                                                                       dose_trachea,
                                                                       dose_lung,
                                                                       dose_gi,
                                                                       inoc,
                                                                       assay.ii,
                                                                       age.ii,
                                                                       sex.ii,
                                                                       sp.ii,
                                                                       lab.ii,
                                                                       lab_effect)
                    
                    #print(sampled_metrics)
                    
                    # Add estimates
                    
                    new.obs <- data.frame(sample_num = ii,
                                          organ_group = organ_group,
                                          organ_idx = organ.ii,
                                          tissue_idx = tissue.ii,
                                          dose_nose = dose_nose,
                                          dose_throat = dose_throat,
                                          dose_trachea = dose_trachea,
                                          dose_lung = dose_lung,
                                          dose_gi = dose_gi,
                                          dose_total = total_dose,
                                          route_idx = route.ii,
                                          sp_idx = sp.ii,
                                          age_idx = age.ii,
                                          sex_idx = sex.ii,
                                          assay_idx = assay.ii,
                                          location_idx = loc.ii,
                                          lab_idx = lab.ii,
                                          percent_positive = sampled_metrics$percent,
                                          first_sample = sampled_metrics$sample_first,
                                          peak_sample = sampled_metrics$sample_peak,
                                          titer_sample = sampled_metrics$sample_titer,
                                          last_sample = sampled_metrics$sample_last
                    )
                    
                    predictions <- rbind(predictions, new.obs)
                    
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  if (lab_effect == "No"){
    predictions$lab_idx <- NA
  }
  
  # Return collection of predictions
  return(predictions)
}


get_metrics_across_cofactors_parallel <- function(fit, n_draws,
                                                  seed = NA,
                                                  dose_options = c(4, 7),
                                                  tissue_options = c(1, 4, 6),
                                                  assay_options = c(1, 4),
                                                  sex_options = c(0, 1),
                                                  age_options = c(1, 2, 3),
                                                  species_options = c(1, 2, 3),
                                                  route_options = c(1, 2, 3, 4, 5),
                                                  lab_effect = "No",
                                                  lab_options = 1,
                                                  rescale_doses = TRUE) {
  
  start_time <- Sys.time()
  
  sample.pars <- get_all_parameter_samples(fit)
  
  if (!is.na(seed)) set.seed(seed)
  sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)

  all_results <- list()
  dose_counter <- 1
  
  for (total_dose in dose_options) {
    
    assign_dose_distribution(total_dose)
    
    dose_nose_vec     <- sapply(1:5, function(i) get(paste0("route", i, ".dose.nose")))
    dose_throat_vec   <- sapply(1:5, function(i) get(paste0("route", i, ".dose.throat")))
    dose_trachea_vec  <- sapply(1:5, function(i) get(paste0("route", i, ".dose.trachea")))
    dose_lung_vec     <- sapply(1:5, function(i) get(paste0("route", i, ".dose.lung")))
    dose_gi_vec       <- sapply(1:5, function(i) get(paste0("route", i, ".dose.gi")))
    
    if (rescale_doses) {
      dose_nose_vec    <- dose_nose_vec / get("max_dose_nose")
      dose_throat_vec  <- dose_throat_vec / get("max_dose_throat")
      dose_trachea_vec <- dose_trachea_vec / get("max_dose_trachea")
      dose_lung_vec    <- dose_lung_vec / get("max_dose_lung")
      dose_gi_vec      <- dose_gi_vec / get("max_dose_gi")
    }
    
    # ---- STEP 1: Create full grid ----
    grid <- expand.grid(
      route = route_options,
      tissue = tissue_options,
      assay = assay_options,
      sex = sex_options,
      age = age_options,
      species = species_options,
      lab = lab_options,
      sample = sample_num,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    # ---- STEP 2: Parallelize across full grid ----
    results_list <- foreach(i = 1:nrow(grid),
                            .export = c("get_individual_parameter_sample")) %dopar% {
                              
                              row <- grid[i,]
                              
                              route.ii  <- row$route
                              tissue.ii <- row$tissue
                              assay.ii  <- row$assay
                              sex.ii    <- row$sex
                              age.ii    <- row$age
                              sp.ii     <- row$species
                              lab.ii    <- row$lab
                              ii        <- row$sample
                              
                              # Organ logic
                              if (tissue.ii %in% c(1,2)) { organ_group <- "URT"; organ.ii <- 1 }
                              if (tissue.ii %in% c(3,4)) { organ_group <- "LRT"; organ.ii <- 2 }
                              if (tissue.ii %in% c(5,6)) { organ_group <- "GI";  organ.ii <- 3 }
                              
                              loc.ii <- ifelse(tissue.ii %% 2 == 1, 0, 1)
                              
                              dose_nose     <- dose_nose_vec[route.ii]
                              dose_throat   <- dose_throat_vec[route.ii]
                              dose_trachea  <- dose_trachea_vec[route.ii]
                              dose_lung     <- dose_lung_vec[route.ii]
                              dose_gi       <- dose_gi_vec[route.ii]
                              
                              inoc <- 0
                              if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0) inoc <- 1
                              if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1) inoc <- 1
                              if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0) inoc <- 1
                              if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1) inoc <- 1
                              if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0) inoc <- 1
                              
                              sampled_metrics <- get_individual_parameter_sample(
                                sample.pars, ii,
                                organ.ii, loc.ii, tissue.ii,
                                route.ii,
                                dose_nose, dose_throat,
                                dose_trachea, dose_lung,
                                dose_gi, inoc,
                                assay.ii, age.ii, sex.ii,
                                sp.ii, lab.ii, lab_effect
                              )
                              
                              data.frame(
                                sample_num = ii,
                                organ_group = organ_group,
                                organ_idx = organ.ii,
                                tissue_idx = tissue.ii,
                                dose_nose = dose_nose,
                                dose_throat = dose_throat,
                                dose_trachea = dose_trachea,
                                dose_lung = dose_lung,
                                dose_gi = dose_gi,
                                dose_total = total_dose,
                                route_idx = route.ii,
                                sp_idx = sp.ii,
                                age_idx = age.ii,
                                sex_idx = sex.ii,
                                assay_idx = assay.ii,
                                location_idx = loc.ii,
                                lab_idx = lab.ii,
                                percent_positive = sampled_metrics$percent,
                                first_pos_median = sampled_metrics$median_first,
                                peak_median = sampled_metrics$median_peak,
                                titer_mean = sampled_metrics$mean_titer,
                                last_median = sampled_metrics$median_last
                              )
                            }
    
    all_results[[dose_counter]] <- bind_rows(results_list)
    dose_counter <- dose_counter + 1
    
    current_elapsed <- difftime(Sys.time(), start_time, units = "mins")
    cat("Finished generating estimates for a total dose of", total_dose, "log10 pfu;\n",
        "elapsed computation time so far is", round(current_elapsed, 2), "minutes\n")
  }

  predictions <- bind_rows(all_results)
  predictions <- calculate_auc(predictions)
  
  if (lab_effect == "No") predictions$lab_idx <- NA
  
  end_time <- Sys.time()
  cat("\nTotal computation time:",
      round(difftime(end_time, start_time, units = "mins"), 2),
      "minutes\n")
  
  return(predictions)
}




get_pairwise_differences <- function(df, cofactor = "route") {
  
  
  if (cofactor != "lab"){
    df <- subset(df, select = -c(lab_idx))
  }
  
  # Remove columns that aren't relevant for cofactor comparison
  df <- subset(df, select = -c(organ_group,
                               organ_idx,
                               dose_nose,
                               dose_throat,
                               dose_trachea,
                               dose_lung,
                               dose_gi,
                               location_idx 
  ))
  
  # Get pairwise differences and keep track of which cofactors were compared
  if (cofactor == "route") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, dose_total, sp_idx, age_idx, sex_idx, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(route_idx), 2, function(x) {
        diff_auc <- auc[route_idx == x[1]] - auc[route_idx == x[2]]
        diff_duration <- duration[route_idx == x[1]] - duration[route_idx == x[2]]
        diff_percent <- percent_positive[route_idx == x[1]] - percent_positive[route_idx == x[2]]
        diff_first <- first_pos_median[route_idx == x[1]] - first_pos_median[route_idx == x[2]]
        diff_peak <- peak_median[route_idx == x[1]] - peak_median[route_idx == x[2]]
        diff_titer <- titer_mean[route_idx == x[1]] - titer_mean[route_idx == x[2]]
        diff_last <- last_median[route_idx == x[1]] - last_median[route_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "dose") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, sp_idx, age_idx, sex_idx, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(dose_total), 2, function(x) {
        diff_auc <- auc[dose_total == x[1]] - auc[dose_total == x[2]]
        diff_duration <- duration[dose_total == x[1]] - duration[dose_total == x[2]]
        diff_percent <- percent_positive[dose_total == x[1]] - percent_positive[dose_total == x[2]]
        diff_first <- first_pos_median[dose_total == x[1]] - first_pos_median[dose_total == x[2]]
        diff_peak <- peak_median[dose_total == x[1]] - peak_median[dose_total == x[2]]
        diff_titer <- titer_mean[dose_total == x[1]] - titer_mean[dose_total == x[2]]
        diff_last <- last_median[dose_total == x[1]] - last_median[dose_total == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "sex") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, sp_idx, age_idx, dose_total, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(sex_idx), 2, function(x) {
        diff_auc <- auc[sex_idx == x[1]] - auc[sex_idx == x[2]]
        diff_duration <- duration[sex_idx == x[1]] - duration[sex_idx == x[2]]
        diff_percent <- percent_positive[sex_idx == x[1]] - percent_positive[sex_idx == x[2]]
        diff_first <- first_pos_median[sex_idx == x[1]] - first_pos_median[sex_idx == x[2]]
        diff_peak <- peak_median[sex_idx == x[1]] - peak_median[sex_idx == x[2]]
        diff_titer <- titer_mean[sex_idx == x[1]] - titer_mean[sex_idx == x[2]]
        diff_last <- last_median[sex_idx == x[1]] - last_median[sex_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "age") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, sp_idx, sex_idx, dose_total, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(age_idx), 2, function(x) {
        diff_auc <- auc[age_idx == x[1]] - auc[age_idx == x[2]]
        diff_duration <- duration[age_idx == x[1]] - duration[age_idx == x[2]]
        diff_percent <- percent_positive[age_idx == x[1]] - percent_positive[age_idx == x[2]]
        diff_first <- first_pos_median[age_idx == x[1]] - first_pos_median[age_idx == x[2]]
        diff_peak <- peak_median[age_idx == x[1]] - peak_median[age_idx == x[2]]
        diff_titer <- titer_mean[age_idx == x[1]] - titer_mean[age_idx == x[2]]
        diff_last <- last_median[age_idx == x[1]] - last_median[age_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "species") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, age_idx, sex_idx, dose_total, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(sp_idx), 2, function(x) {
        diff_auc <- auc[sp_idx == x[1]] - auc[sp_idx == x[2]]
        diff_duration <- duration[sp_idx == x[1]] - duration[sp_idx == x[2]]
        diff_percent <- percent_positive[sp_idx == x[1]] - percent_positive[sp_idx == x[2]]
        diff_first <- first_pos_median[sp_idx == x[1]] - first_pos_median[sp_idx == x[2]]
        diff_peak <- peak_median[sp_idx == x[1]] - peak_median[sp_idx == x[2]]
        diff_titer <- titer_mean[sp_idx == x[1]] - titer_mean[sp_idx == x[2]]
        diff_last <- last_median[sp_idx == x[1]] - last_median[sp_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "assay") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, age_idx, sex_idx, dose_total, sp_idx) %>%
      summarise(pairwise_diff = combn(unique(assay_idx), 2, function(x) {
        diff_auc <- auc[assay_idx == x[1]] - auc[assay_idx == x[2]]
        diff_duration <- duration[assay_idx == x[1]] - duration[assay_idx == x[2]]
        diff_percent <- percent_positive[assay_idx == x[1]] - percent_positive[assay_idx == x[2]]
        diff_first <- first_pos_median[assay_idx == x[1]] - first_pos_median[assay_idx == x[2]]
        diff_peak <- peak_median[assay_idx == x[1]] - peak_median[assay_idx == x[2]]
        diff_titer <- titer_mean[assay_idx == x[1]] - titer_mean[assay_idx == x[2]]
        diff_last <- last_median[assay_idx == x[1]] - last_median[assay_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  else if (cofactor == "tissue") {
    df.pair <- df %>%
      group_by(sample_num, assay_idx, route_idx, age_idx, sex_idx, dose_total, sp_idx) %>%
      summarise(pairwise_diff = combn(unique(tissue_idx), 2, function(x) {
        diff_auc <- auc[tissue_idx == x[1]] - auc[tissue_idx == x[2]]
        diff_duration <- duration[tissue_idx == x[1]] - duration[tissue_idx == x[2]]
        diff_percent <- percent_positive[tissue_idx == x[1]] - percent_positive[tissue_idx == x[2]]
        diff_first <- first_pos_median[tissue_idx == x[1]] - first_pos_median[tissue_idx == x[2]]
        diff_peak <- peak_median[tissue_idx == x[1]] - peak_median[tissue_idx == x[2]]
        diff_titer <- titer_mean[tissue_idx == x[1]] - titer_mean[tissue_idx == x[2]]
        diff_last <- last_median[tissue_idx == x[1]] - last_median[tissue_idx == x[2]]
        data.frame(cofactor1 = x[1], 
                   cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration,
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  
  else if (cofactor == "lab") {
    df.pair <- df %>%
      group_by(sample_num, tissue_idx, route_idx, age_idx, sex_idx, dose_total, sp_idx, assay_idx) %>%
      summarise(pairwise_diff = combn(unique(lab_idx), 2, function(x) {
        diff_auc <- auc[lab_idx == x[1]] - auc[lab_idx == x[2]]
        diff_duration <- duration[lab_idx == x[1]] - duration[lab_idx == x[2]]
        diff_percent <- percent_positive[lab_idx == x[1]] - percent_positive[lab_idx == x[2]]
        diff_first <- first_pos_median[lab_idx == x[1]] - first_pos_median[lab_idx == x[2]]
        diff_peak <- peak_median[lab_idx == x[1]] - peak_median[lab_idx == x[2]]
        diff_titer <- titer_mean[lab_idx == x[1]] - titer_mean[lab_idx == x[2]]
        diff_last <- last_median[lab_idx == x[1]] - last_median[lab_idx == x[2]]
        data.frame(cofactor1 = x[1], cofactor2 = x[2], 
                   diff_percent = diff_percent,
                   diff_first = diff_first,
                   diff_peak = diff_peak,
                   diff_titer = diff_titer,
                   diff_last = diff_last,
                   diff_duration = diff_duration, 
                   diff_auc = diff_auc)
      }, simplify = FALSE)) %>%
      unnest(pairwise_diff)
  }
  
  df.pair$cofactor <- cofactor
  return(df.pair)
}


# Extracts all event times for individuals with ID names
extract_times_with_ids <- function(dat) {
  
  dat.times <- dat[1, ]
  dat.times$first_upper_bound <- NA
  dat.times$first_lower_bound <- NA
  dat.times$peak_upper_bound <- NA
  dat.times$peak_observed_time <- NA
  dat.times$peak_observed_titer <- NA
  dat.times$peak_lower_bound <- NA
  dat.times$last_upper_bound <- NA
  dat.times$last_lower_bound <- NA
  
  dat.all <- dat[1,]
  dat.all$first_upper_bound <- NA
  dat.all$first_lower_bound <- NA
  dat.all$peak_upper_bound <- NA
  dat.all$peak_observed_time <- NA
  dat.all$peak_observed_titer <- NA
  dat.all$peak_lower_bound <- NA
  dat.all$last_upper_bound <- NA
  dat.all$last_lower_bound <- NA
  
  dat.times <- dat.times[0,]
  dat.all <- dat.all[0,]
  
  # Loop over individual ID names, locations, and assay types 
  
  for (indiv.ii in unique(dat$indiv)) {
    dat.indiv <- subset(dat, indiv == indiv.ii)
    
    for (loc.grp.ii in unique(dat.indiv$location_grp)) {
      dat.loc.grp <- subset(dat.indiv, location_grp == loc.grp.ii) 
      
      for (loc.subgrp.ii in unique(dat.loc.grp$location_subgrp)) {
        dat.loc.subgrp <- subset(dat.loc.grp, location_subgrp == loc.subgrp.ii)
        
        for (assay.ii in unique(dat.loc.subgrp$rna_type)) {
          dat.assay.high <- subset(dat.loc.subgrp, rna_type == assay.ii)
          
          for (tg.ii in unique(dat.assay.high$pcr_target_gene)) {
            
            indiv_sample <- paste(indiv.ii, loc.subgrp.ii,
                                  assay.ii, tg.ii, sep = "_")
            
            first_upper_bound <- NA
            first_lower_bound <- NA
            peak_upper_bound <- NA
            peak_observed_time <- NA
            peak_observed_titer <- NA
            peak_lower_bound <- NA
            last_upper_bound <- NA
            last_lower_bound <- NA
            
            
            # Have to deal with annoying NAs
            if (is.na(tg.ii)) {
              dat.assay <- subset(dat.assay.high, is.na(pcr_target_gene))
            }
            else {
              dat.assay <- subset(dat.assay.high, pcr_target_gene == tg.ii)
            }
            
            # Check which sample types are available 
            
            # Invasive only data -----------------------------------------------
            
            if ("Non-invasive" %notin% dat.assay$sample_type & 
                "Invasive" %in% dat.assay$sample_type) {
              
              #cat("Running on", indiv_sample, " for Inv only \n")
              
              # Invasive only data can't contribute a peak time
              
              # Throw an error if there's more than 1 DPI (can't happen for invasive)
              if (length(unique(dat.assay$day_post_infection)) > 1){
                cat("ERROR: More than 1 DPI for invasive only sample of individual ", 
                    indiv.ii, "\n")
              }
              
              # Check if there is at least one positive among the invasives
              if (1 %in% unique(dat.assay$pos_value)) {
                # Lower bound is zero since there are no other sampling days
                first_lower_bound <- 0
                # Upper bound is the necropsy day
                first_upper_bound <- unique(as.numeric(dat.assay$day_post_infection))
                last_lower_bound <- first_upper_bound
                last_upper_bound <- Inf
              }
              else {
                # Lower bound is the day of necropsy
                first_lower_bound <- unique(as.numeric(dat.assay$day_post_infection))
                # Upper bound is infinity
                first_upper_bound <- Inf
                
                # But if the necropsy day is >15 days after inoculation, don't include
                #if (first_lower_bound >= 15) {
                #  first_lower_bound <- NA
                #  first_upper_bound < NA
                #}
                
              }
              
              
            }
            
            # Non-invasive only data -------------------------------------------
            
            else if ("Non-invasive" %in% dat.assay$sample_type & 
                     "Invasive" %notin% dat.assay$sample_type) {
              
              #cat("Running on", indiv_sample, " for Inv only \n")
              
              # Check if there are any positives
              if (1 %in% unique(dat.assay$pos_value)) {
                
                # If the DPI wasn't reported explicitly, flag it
                if (NA %in% as.numeric(dat.assay$day_post_infection)) {
                  
                  if (length(unique(as.numeric(dat.assay$day_post_infection))) > 1) {
                    cat("ERROR: Only some DPI are coded as ranges, while others are coded explicitely\n")
                  }
                  
                  first_upper_bound <- "FLAG"
                  first_lower_bound <- "FLAG"
                  peak_upper_bound <- "FLAG"
                  peak_lower_bound <- "FLAG"
                  last_upper_bound <- "FLAG"
                  last_lower_bound <- "FLAG"
                  
                }
                else {
                  
                  ## FIRST POSITIVE TIMES
                  
                  # Find the upper bound (smallest DPI where the individual is positive)
                  first_upper_bound <- min(as.numeric(dat.assay$day_post_infection[dat.assay$pos_value == 1]))
                  
                  # Find the lower bound (previous sampling time)
                  dpis_tested <- unique(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) < first_upper_bound])
                  if (length(dpis_tested) == 0) {
                    # Lower bound is inoculation time if no other sampling days 
                    first_lower_bound <- 0 
                  }
                  else {
                    # Lower bound is maximum sampling time that's less than
                    #     the first observed positive
                    first_lower_bound <- max(as.numeric(dpis_tested))
                  }
                  
                  
                  ## PEAK TIMES: must have quantitative information
                  
                  # Must have some quantitative information
                  if (FALSE %in% is.na(unique(as.numeric(dat.assay$value)))) {
                    
                    # Must have at least two sample times or 
                    #   have explicit peak range available
                    if ((length(unique(dat.assay$day_post_infection)) > 1) | 
                        (TRUE %in% str_detect(dat.assay$unit_subgrp, "Peak|peak"))) {
                      
                      # Throw error if somehow only one sample available
                      if (nrow(dat.assay) == 1) {
                        cat(blue("NOTE: only one sample time available:",
                                 dat.assay$day_post_infection, "\n"))
                      }
                      
                      # Find the day with the peak observed viral load
                      # May need to restrict to > day 1 if inoculated tissue
                      peak_load <- max(as.numeric(dat.assay$value), na.rm = TRUE)
                      peak_dpi <- dat.assay$day_post_infection[!is.na(as.numeric(dat.assay$value)) &
                                                                 as.numeric(dat.assay$value) == peak_load]
                      
                      # Check whether multiple days have the same viral load 
                      if (length(peak_dpi) > 1) {
                        cat(blue("Note: multiple days have the same peak load, selecting the first one.\n"))
                        peak_dpi <- min(as.numeric(peak_dpi))
                      }
                      
                      # Check whether it's the same day as the first positive
                      #if (peak_dpi == first_upper_bound){
                      #  cat(red(indiv_sample, ": Peak time is also the first positive day. Restrict this?!\n"))
                      #}
                      
                      # Find the lower & upper bounds as the previous & next sample times
                      peak_lower_bound <- max(as.numeric(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) < as.numeric(peak_dpi)]),
                                              na.rm = TRUE)
                      peak_upper_bound <- min(as.numeric(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) > as.numeric(peak_dpi)]),
                                              na.rm = TRUE) # will be infinity if there is no such time
                      peak_observed_time <- peak_dpi
                      peak_observed_titer <- peak_load
                      
                      # If there is no previous observation, then 
                      #    it's the inoculation time
                      if (peak_lower_bound == -Inf) {
                        peak_lower_bound <- 0
                      }
                    }
                  }
                  
                  
                  ## LAST POSITIVE TIMES
                  
                  # Lower bound is the last day they're observed positive 
                  last_lower_bound <- max(as.numeric(dat.assay$day_post_infection[dat.assay$pos_value == 1])) 
                  
                  # Check whether the last positive is also the first positive
                  if (last_lower_bound == first_upper_bound) {
                    # If the LB is the same day as the 1st positive, then there's
                    #    not enough information, so we exclude this time
                    #cat(blue("NOTE: the last positive (", last_lower_bound, 
                    #         ") is the first positive (", first_upper_bound,
                    #         ") for NI data!\n\n"))
                  } 
                  
                  if (last_lower_bound < first_upper_bound) {
                    # Throw an error if the last positive happens before the 
                    #     the first positive 
                    cat(red("ERROR: Last positive time happens before the first positive!!!!\n\n"))
                  } 
                  else { 
                    # Find all days tested after the last positive time
                    later_dpis_tested <- unique(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) > last_lower_bound])
                    later_dpis_tested <- later_dpis_tested[!is.na(later_dpis_tested)]
                    
                    # Check if there are any later negative samples
                    if (length(later_dpis_tested ) == 0) {
                      # If not, then it's right censored
                      last_upper_bound <- Inf
                    }
                    else {
                      last_upper_bound <- min(as.numeric(later_dpis_tested))
                    }
                  }
                }
              }
              
              # If there are no positives
              else {
                
                if (NA %in% as.numeric(dat.assay$day_post_infection)) {
                  
                  if (length(unique(as.numeric(dat.assay$day_post_infection))) > 1) {
                    cat("ERROR: Only some DPI are coded as ranges, while others are coded explicitely\n")
                  }
                  
                  first_upper_bound <- Inf
                  first_lower_bound <- "FLAG"
                  peak_upper_bound <- NA
                  peak_lower_bound <- NA
                  last_upper_bound <- NA
                  last_lower_bound <- NA
                }
                else {
                  first_lower_bound <- max(as.numeric(dat.assay$day_post_infection))
                  first_upper_bound <- Inf
                  peak_upper_bound <- NA
                  peak_lower_bound <- NA
                  last_upper_bound <- NA
                  last_lower_bound <- NA
                }
                
              }
            }
            
            # Both  ------------------------------------------------------------
            
            else if ("Non-invasive" %in% dat.assay$sample_type & 
                     "Invasive" %in% dat.assay$sample_type) {
              
              # Check if there is at least one positive
              if (1 %in% unique(dat.assay$pos_value)) {
                
                # Non-invasives take priority, because of higher resolution on timing
                if ("Non-invasive" %in% subset(dat.assay, pos_value == 1)$sample_type) {
                  dat.ni <- subset(dat.assay, sample_type == "Non-invasive")
                  
                  # If DPI not reported explicitly, flag it
                  if (NA %in% as.numeric(dat.ni$day_post_infection)) {
                    
                    if (length(unique(dat.ni$day_post_infection)) > 1) {
                      cat("ERROR: Only some DPI are coded as ranges, while others are coded explicitely\n")
                    }
                    first_upper_bound <- "FLAG"
                    first_lower_bound <- "FLAG"
                  }
                  else {
                    
                    ## FIRST POSITIVE
                    first_upper_bound <- min(as.numeric(dat.ni$day_post_infection[dat.ni$pos_value == 1]))
                    
                    # Find the number of days tested before the first positive
                    dpis_tested <- unique(dat.ni$day_post_infection[as.numeric(dat.ni$day_post_infection) < first_upper_bound])
                    
                    # Set the previous sampling time
                    if (length(dpis_tested) == 0) {
                      first_lower_bound <- 0 
                    }
                    else {
                      first_lower_bound <- max(as.numeric(dpis_tested))
                    }
                    
                    ## PEAK TIMES: must have quantitative information
                    
                    # Must have some quantitative information
                    if (FALSE %in% is.na(unique(as.numeric(dat.ni$value)))) {
                      
                      # Must have at least two sample times or 
                      #   have explicit peak range available
                      if ((length(unique(dat.ni$day_post_infection)) > 1) | 
                          (TRUE %in% str_detect(dat.ni$unit_subgrp, "Peak|peak"))) {
                        
                        # Throw error if somehow only one sample available
                        if (nrow(dat.ni) == 1) {
                          cat(blue("NOTE: only one sample time available:",
                                   dat.ni$day_post_infection, "\n"))
                        }
                        
                        # Find the day with the peak observed viral load
                        # May need to restrict to > day 1 if inoculated tissue
                        peak_load <- max(as.numeric(dat.ni$value), na.rm = TRUE)
                        peak_dpi <- dat.ni$day_post_infection[!is.na(as.numeric(dat.ni$value)) &
                                                                as.numeric(dat.ni$value) == peak_load]
                        
                        # Check whether multiple days have the same viral load 
                        if (length(peak_dpi) > 1) {
                          cat(blue("Note: multiple days have the same peak load, selecting the first one.\n"))
                          peak_dpi <- min(as.numeric(peak_dpi))
                        }
                        
                        # Check whether it's the same day as the first positive
                        #if (peak_dpi == first_upper_bound){
                        #  cat(red(indiv_sample, ": Peak time is also the first positive day. Restrict this?!\n"))
                        #}
                        
                        # Find the lower & upper bounds as the previous & next sample times
                        peak_lower_bound <- max(as.numeric(dat.ni$day_post_infection[as.numeric(dat.ni$day_post_infection) < as.numeric(peak_dpi)]),
                                                na.rm = TRUE)
                        peak_upper_bound <- min(as.numeric(dat.ni$day_post_infection[as.numeric(dat.ni$day_post_infection) > as.numeric(peak_dpi)]),
                                                na.rm = TRUE) # will be infinity if there is no such time
                        peak_observed_time <- peak_dpi
                        peak_observed_titer <- peak_load
                        
                        # If there is no previous observation, then 
                        #    it's the inoculation time
                        if (peak_lower_bound == -Inf) {
                          peak_lower_bound <- 0
                        }
                      }
                    }
                    
                    ## LAST POSITIVE TIMES
                    
                    # Lower bound is the last day they're observed positive 
                    #    for this we include their invasive samples!
                    last_lower_bound <- max(as.numeric(dat.assay$day_post_infection[dat.assay$pos_value == 1])) 
                    
                    # check if the last positive is the same as the first positive
                    if (last_lower_bound == first_upper_bound) {
                      # If the LB is the same day as the 1st positive, then there's
                      #    not enough information, so we exclude this time
                      #cat(blue("NOTE: the last positive (", last_lower_bound, 
                      #         ") is the first positive (", first_upper_bound,
                      #         ") for NI + I data!\n\n"))
                      ##last_lower_bound <- NA
                      #last_upper_bound <- NA
                    } 
                    
                    if (last_lower_bound < first_upper_bound) {
                      # Throw an error if the last positive happens before the 
                      #     the first positive 
                      cat(red("ERROR: Last positive time happens before the first positive!!!!\n\n"))
                    } 
                    else { 
                      # Find all days tested after the last positive time
                      later_dpis_tested <- unique(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) > last_lower_bound])
                      later_dpis_tested <- later_dpis_tested[!is.na(later_dpis_tested)]
                      
                      # Check if there are any later negative samples
                      if (length(later_dpis_tested ) == 0) {
                        # If not, then it's right censored
                        last_upper_bound <- Inf
                      }
                      else {
                        last_upper_bound <- min(as.numeric(later_dpis_tested))
                      }
                    }
                  }
                }
                
                # If no positive non-invasive, must have positive invasive
                else {
                  
                  # Can only contribute first positives (not peak / last positives)
                  first_upper_bound <- max(as.numeric(dat.assay$day_post_infection))
                  
                  # Find the previous sampling time (if it exists),
                  #    otherwise use the inoculation time
                  if (length(unique(dat.assay$day_post_infection)) > 1) {
                    first_lower_bound <- max(as.numeric(dat.assay$day_post_infection[as.numeric(dat.assay$day_post_infection) < first_upper_bound]))
                  }
                  else {
                    first_lower_bound <- 0
                  }
                }
              }
              
              # If there's no positives, use the last day post infection for R censoring
              else {
                first_lower_bound <- max(as.numeric(dat.assay$day_post_infection))
                first_upper_bound <- Inf
                
                if (NA %in% unique(as.numeric(dat.assay$day_post_infection))) {
                  cat("The unique sampling days to consider are:", 
                      unique(dat.assay$day_post_infection),
                      "\n")
                }
              }
            }
            
            # Handle oddity
            if (is.na(first_lower_bound) & first_upper_bound == Inf){
              first_upper_bound <- NA
            }
            
            dat.assay$first_upper_bound <- first_upper_bound
            dat.assay$first_lower_bound <- first_lower_bound
            dat.assay$peak_upper_bound <- peak_upper_bound
            dat.assay$peak_observed_time <- peak_observed_time
            dat.assay$peak_observed_titer <- peak_observed_titer
            dat.assay$peak_lower_bound <- peak_lower_bound
            dat.assay$last_upper_bound <- last_upper_bound
            dat.assay$last_lower_bound <- last_lower_bound
            dat.assay$indiv_sample <- indiv_sample
            
            # Add data to dataframe
            dat.times <- rbind(dat.times, dat.assay[1, ])
            dat.all <- rbind(dat.all, dat.assay)
            
          }
        }
      }
    }
  }
  
  colnames(dat.times)[which(colnames(dat.times) == "rna_type")] <- "assay_type"
  
  both_dfs <- list(dat.times = dat.times,
                   dat.all = dat.all)
  return(both_dfs)
}

# Extracts all event times for individuals without ID names (relies on helper fns)
extract_times_without_ids <- function(dat.noid) {
  
  dat.noid$first_lower_bound <- NA
  dat.noid$first_upper_bound <- NA
  dat.noid$peak_lower_bound <- NA
  dat.noid$peak_observed_time <- NA
  dat.noid$peak_observed_titer <- NA
  dat.noid$peak_upper_bound <- NA
  dat.noid$last_lower_bound <- NA
  dat.noid$last_upper_bound <- NA
  dat.noid$group_number <- NA
  
  dat.times <- dat.noid[0, ]
  
  group_number <- 0
  
  # Loop over every article, species, treatment group, sex, location, and assay ------
  for (article.ii in unique(dat.noid$article)) {
    dat.article <- subset(dat.noid, article == article.ii)
    
    cat("\nNow extracting for: ", article.ii, "\n")
    
    for (species.ii in unique(dat.article$animal_species)) {
      dat.species <- subset(dat.article, animal_species == species.ii)
      
      for (treat.ii in unique(dat.species$treatment)) {
        dat.treat <- subset(dat.species, treatment == treat.ii)
        
        for (sex.ii in unique(dat.treat$sex)) {
          dat.sex <- subset(dat.treat, sex == sex.ii)
          
          for (age.ii in unique(dat.sex$age_class_grp)) {
            dat.age <- subset(dat.sex, age_class_grp == age.ii)
            
            for (route.ii in unique(dat.age$inoc_route_rep)){
              dat.route <- subset(dat.age, inoc_route_rep == route.ii)
              
              for (dose.ii in unique(dat.route$inoc_dose_total_pfu)) {
                dat.dose <- subset(dat.route, inoc_dose_total_pfu == dose.ii)
                
                for (loc.grp.ii in unique(dat.dose$location_grp)) {
                  dat.loc.grp <- subset(dat.dose, location_grp == loc.grp.ii) 
                  
                  for (loc.subgrp.ii in unique(dat.loc.grp$location_subgrp)) {
                    dat.loc.subgrp <- subset(dat.loc.grp, location_subgrp == loc.subgrp.ii)
                    
                    for (assay.ii in unique(dat.loc.subgrp$rna_type)) {
                      dat.assay.high <- subset(dat.loc.subgrp, rna_type == assay.ii)
                      
                      for (tg.ii in unique(dat.assay.high$pcr_target_gene)) {
                        
                        
                        if (is.na(tg.ii)) {
                          dat.assay <- subset(dat.assay.high, is.na(pcr_target_gene))
                        }
                        else {
                          dat.assay <- subset(dat.assay.high, pcr_target_gene == tg.ii)
                        }
                        
                        # Add the group number
                        group_number <- group_number + 1
                        dat.assay$group_number <- group_number
                        
                        #cat("... and for group #", group_number, "\n")
                        
                        # Check which sample types are available ---------------
                        
                        ## Only Invasive data ----------------------------------
                        if ("Non-invasive" %notin% dat.assay$sample_type) {
                          
                          # Can only contribute first positive times
                          dat.first <- get_first_positives_only_invasives_noid(dat.assay)
                          dat.times <- rbind(dat.times, dat.first)
                          
                        }
                        
                        ## Only Non-invasive data ------------------------------
                        else if ("Invasive" %notin% dat.assay$sample_type) {
                          
                          dat.first <- get_first_positives_only_noninvasives_noid(dat.assay)
                          dat.peak <- get_peak_times_noid(subset(dat.assay, !is.na(as.numeric(day_post_infection))))
                          
                          if (nrow(dat.peak) > 0) {
                            dat.consolidated <- consolidate_times_noid(dat.first, dat.peak)
                          }
                          else {
                            dat.consolidated <- dat.first
                          }
                          
                          dat.times <- rbind(dat.times, dat.consolidated)
                          
                        }
                        
                        ## Both Invasive & Non-Invasive data -------------------
                        else {
                          
                          dat.first <- get_first_positives_both_noid(dat.assay)
                          dat.peak <- get_peak_times_noid(subset(dat.assay, !is.na(as.numeric(day_post_infection))))
                          
                          if (nrow(dat.peak) > 0) {
                            dat.consolidated <- consolidate_times_noid(dat.first, dat.peak)
                          }
                          else {
                            dat.consolidated <- dat.first
                          }
                          
                          dat.times <- rbind(dat.times, dat.consolidated)
                          
                        }
                      } # the parenthesis above this one closes off the Inv & NonIinv samples 
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  colnames(dat.times)[which(colnames(dat.times) == "rna_type")] <- "assay_type"
  
  # Handle the few individuals with ranges, being maximally conservative
  dat.times$first_lower_bound[dat.times$first_lower_bound == "Unknown (7/8)"] <- 7
  dat.times$first_lower_bound[dat.times$first_lower_bound == "Unknown (6-8)"] <- 6
  
  dat.times$first_upper_bound[dat.times$first_upper_bound == "Unknown (7/8)"] <- 8
  dat.times$first_upper_bound[dat.times$first_upper_bound == "Unknown (6-8)"] <- 8
  
  # Make the bounds numerics
  dat.times$first_lower_bound <- as.numeric(dat.times$first_lower_bound)
  dat.times$first_upper_bound <- as.numeric(dat.times$first_upper_bound)
  
  
  return(dat.times)
}

# Extracts first positives for no ID individuals, when they have both
#     invasive and non-invasive data
get_first_positives_both_noid <- function(dat) {
  
  dat.inv <- subset(dat, sample_type == "Invasive")
  dat.noninv <- subset(dat, sample_type == "Non-invasive")
  
  dat.inv.new <- get_first_positives_only_invasives_noid(dat.inv)
  dat.noninv.new <- get_first_positives_only_noninvasives_noid(dat.noninv)
  
  cat(red("Invasive dataframe gives", nrow(dat.inv.new), "observations",
      "while the non-invasive dataframe gives", nrow(dat.noninv.new),
      "observations.\n\n"))
  
  # If there's no noninvasive data with known sampling times, choose the
  #    invasive data
  if (nrow(dat.noninv.new) == 0) {
    dat.new <- dat.inv.new
    cat(blue("Proceeding with the invasive data frame...\n"))
  }
  ## If somehow the invasive data has more samples... could try to get more info...
  ##    but for simplicity right now, not going to do this
  #else if (nrow(dat.inv.new) > nrow(dat.noninv.new)) {
  #  number_inv_positives <- nrow(subset(dat.inv.new, first_upper_bound != Inf))
  #  number_noninv_positives <- nrow(subset(dat.noninv.new, first_upper_bound != Inf))
  #  cat(blue("The invasive data has", number_inv_positives, "positives",
  #           "while the noninvasive data has", number_noninv_positives, "\n"))
  #}
  # Otherwise, we prefer the noninvasive data
  else {
    dat.new <- dat.noninv.new
  }
  
  return(dat.new)
  
}

# Extracts first positives for no ID individuals, when they only have
#     non-invasive data
get_first_positives_only_noninvasives_noid <- function(dat) {
  
  dat.fill <- dat[0, ]
  
  # Remove unknowns, because they complicate extracting timing
  dat <- subset(dat, !str_detect(day_post_infection, "Unknown"))
  
  # Check for whether there are ever any positives...
  #     If there's none, find the day with the most late negatives
  if (nrow(dat) > 0 & 1 %notin% dat$pos_value) {
    
    # Fill in a dataframe for the number of individuals tested 
    #   on each sample day
    dpi.df <- data.frame(sample_rep = character(),
                           dpi = numeric(),
                           n_neg = numeric(),
                           n_total = numeric())
    unique_dpis <- unique(dat$day_post_infection) # unique sampling times
    
    for (dpi.ii in unique_dpis) {
      dpi.sub <- subset(dat, day_post_infection == dpi.ii)
      
      # If multiple sample types contribute, need to pick the one with the
      #   most samples for any given day
      if (length(unique(dpi.sub$sample_rep)) > 1) {
        dpi.sub.s1 <- subset(dpi.sub, sample_rep == unique(dpi.sub$sample_rep)[1])
        dpi.sub.s2 <- subset(dpi.sub, sample_rep == unique(dpi.sub$sample_rep)[2])
        if (sum(dpi.sub.s1$n_neg) >= sum(dpi.sub.s2$n_neg)) {
          dpi.sub <- dpi.sub.s1
        }
        else {dpi.sub <- dpi.sub.s2}
      }
      
      dpi.df <- rbind(dpi.df, 
                        data.frame(sample_rep = unique(dpi.sub$sample_rep),
                                   dpi = dpi.ii,
                                   n_neg = nrow(subset(dpi.sub, pos_value == 0)),
                                   n_total = nrow(dpi.sub)))
    }
    
    dpi.df <- dpi.df[order(dpi.df$n_neg, dpi.df$dpi, decreasing = TRUE), ]
    
    number_of_negatives <- dpi.df$n_neg[1]
    
    # Add to dataframe
    dat.new <- dat[rep(1, number_of_negatives), ]
    dat.new$first_lower_bound <- dpi.df$dpi[1] # last observation day with the most negatives
    dat.new$first_upper_bound <- Inf
    
    print(dpi.df)
    cat("We are adding", number_of_negatives, "right censored individuals.\n")
  
  }
  
  else if (nrow(dat) > 0 & 1 %in% dat$pos_value) {
    
    # Create a data frame with sampling times & # positives
    dpi.df <- data.frame(sample_rep = character(),
                         sample_type = character(),
                         dpi = numeric(),
                         n_pos = numeric(),
                         n_total = numeric())
    
    # Fill it in!
    for (dpi.ii in unique(dat$day_post_infection)) {
      
      dat.dpi <- subset(dat, day_post_infection == dpi.ii)
      
      # Some studies have multiple sample names...
      #    We will take the one with the most positives
      if (length(unique(dat.dpi$sample_rep)) > 1) {
        dat.dpi.s1 <- subset(dat.dpi, sample_rep == unique(dat.dpi$sample_rep)[1])
        dat.dpi.s2 <- subset(dat.dpi, sample_rep == unique(dat.dpi$sample_rep)[2])
        
        if (sum(dat.dpi.s1$pos_value) >= sum(dat.dpi.s2$pos_value)) {
          dat.dpi <- dat.dpi.s1
        }
        else {dat.dpi <- dat.dpi.s2}
      }
      
      dpi.df <- rbind(dpi.df,
                        data.frame(sample_rep = unique(dat.dpi$sample_rep),
                                   sample_type = "Non-invasive",
                                   dpi = dpi.ii,
                                   n_pos = sum(dat.dpi$pos_value),
                                   n_total = nrow(dat.dpi)))
    }
    
    dpi.df <- dpi.df[order(as.numeric(dpi.df$dpi)), ]
    #print(dpi.df)
    
    # Now extract times based on that dataframe...
    # If there's only one row, that's what we must use. Or if the 1st day has
    #    the maximum number of positives
    if (nrow(dpi.df) == 1 | dpi.df$n_pos[1] == max(dpi.df$n_pos)) {
      
      number_of_positives <- dpi.df$n_pos[1]
      number_of_negatives <- dpi.df$n_total[1] - dpi.df$n_pos[1]
      
      positive_lower_bounds <- rep(0, number_of_positives)
      positive_upper_bounds <- rep(dpi.df$dpi[1], number_of_positives)
      
      negative_lower_bounds <- rep(dpi.df$dpi[1], number_of_negatives)
      negative_upper_bounds <- rep(Inf, number_of_negatives)
      
      #cat(blue("Adding", number_of_positives, "positives and",
      #         number_of_negatives, "negatives.\n"))
      #cat(blue("The positive upper bounds are", 
      #         rep(dpi.df$dpi[1], number_of_positives), "\n"))
    }
    # If there's more than one row, loop through them to get as many positives
    #    as possible
    else {
      
      # To keep track of positives & negatives & max possible observations
      number_of_positives <- 0
      number_of_negatives <- 0
      max_tested <- max(dpi.df$n_total)
      
      # A list to store the bounds
      positive_lower_bounds <- c()
      positive_upper_bounds <- c()
      negative_lower_bounds <- c()
      negative_upper_bounds <- c()
      row_num <- 1

      # Find the first row with positives to start looping over
      first_positive_row <- which(dpi.df$n_pos > 0)[1] 
      
      # Find the previous day sampled
      if (first_positive_row == 1){
        previous_sampling_day <- 0
      } 
      else {
        previous_sampling_day <- dpi.df$dpi[first_positive_row - 1]
      }
      
      
      # Loop over all sample days, testing whether there's ever more positive
      #   individuals than have been previously observed
      for (row_num in first_positive_row:nrow(dpi.df)) {
        
        row_positives <- dpi.df$n_pos[row_num]
        
        if (row_positives > number_of_positives) {
          
          # Find the number of new positives and the current DPI
          number_new_positives <- row_positives - number_of_positives
          current_dpi <- dpi.df$dpi[row_num]
          
          # Add lower and upper bounds to list
          positive_lower_bounds <- c(positive_lower_bounds,
                                     rep(previous_sampling_day, number_new_positives))
          positive_upper_bounds <- c(positive_upper_bounds,
                                     rep(current_dpi, number_new_positives))
          
          ## Print out what's happening to double check
          #cat(blue("We added", number_new_positives, 
          #         "new positives with the following bounds: [",
          #         previous_sampling_day, ",", current_dpi, "]\n"))
          
          
          # Check whether fewer individuals were tested on this day than the next day
          #    because if so, then need to keep the same lower bound for future positives
          if (row_num != nrow(dpi.df) & dpi.df$n_total[row_num] < dpi.df$n_total[row_num + 1]) {
            # don't change the previous sampling day
          }
          else {
            # Set previous sampling day to the current day (to bound future positives)
            previous_sampling_day <- current_dpi
          }
        
          # Reset the total number of observed positives
          number_of_positives <- number_of_positives + number_new_positives
          

        }
      }
      
      # Now determine whether any right censored individuals can be added
      #     This only makes sense if we haven't already observed all possible positives
      if (number_of_positives < max_tested) {
        
        # Calculate how many total positives are ever observed
        total_ever_positive <- sum(dpi.df$n_pos)
        
        # If this is fewer than the maximum number of individuals tested,
        #    then we can use the last positive day as the lower bound for those
        #    individuals
        if (total_ever_positive < max_tested) {
          
          number_of_negatives <- max_tested - total_ever_positive
          last_day_with_positives <- max(as.numeric(dpi.df$dpi[dpi.df$n_pos > 0]))
          
          negative_lower_bounds <- c(negative_lower_bounds, 
                                     rep(last_day_with_positives, number_of_negatives))
          negative_upper_bounds <- c(negative_upper_bounds,
                                     rep(Inf, number_of_negatives))
          
          #cat(blue("We are adding", number_of_negatives, "negatives with bounds: [", 
          #         last_day_with_positives, ",", Inf, "]\n"))
          
        }
        
        # If we still have individuals unaccounted for (either positive or negative)
        #    we can right censor them based on the first day there are ever positives
        if (sum(number_of_negatives, number_of_positives) < max_tested) {
          
          number_unaccounted_for <- max_tested - sum(number_of_negatives, number_of_positives)
          first_ever_positive_day <- min(as.numeric(dpi.df$dpi[dpi.df$n_pos > 0]))
          
          number_of_negatives <- number_of_negatives + number_unaccounted_for
          negative_lower_bounds <- c(negative_lower_bounds, 
                                     rep(first_ever_positive_day, number_unaccounted_for))
          negative_upper_bounds <- c(negative_upper_bounds,
                                     rep(Inf, number_unaccounted_for))
          
          #cat("We are adding", number_unaccounted_for, "with bounds: [",
          #    first_ever_positive_day, ",", Inf, "]\n")
        }
        
      }
      
    }

    # Add to dataframe
    dat.new.pos <- dat[rep(1, number_of_positives), ]
    dat.new.pos$first_lower_bound <- positive_lower_bounds
    dat.new.pos$first_upper_bound <- positive_upper_bounds
    
    if (number_of_negatives > 0) {
      dat.new.neg <- dat[rep(1, number_of_negatives), ]
      dat.new.neg$first_lower_bound <- negative_lower_bounds
      dat.new.neg$first_upper_bound <- negative_upper_bounds
    }
    else {
      dat.new.neg <- dat[0, ]
    }
    
    dat.new <- rbind(dat.new.pos, dat.new.neg)
    
    cat("We are adding", number_of_positives, "interval censored individuals",
        "and", number_of_negatives, "right censored individuals.\n")
    
  }
  
  else {
    cat(blue("We have no days with known sampling time! Returning an empty dataframe."))
    dat.new <- dat[0, ]
  }
  
  return(dat.new)
 
}

# Extracts first positives for no ID individuals, when they only have
#     non-invasive data
get_first_positives_only_invasives_noid <- function(dat) {
  
  dat.fill <- dat[0, ]
  
  # Check each DPI, which have different individuals by design 
  for (dpi.ii in unique(dat$day_post_infection)) {
    
    dat.dpi <- subset(dat, day_post_infection == dpi.ii)
    
    # Fill in a data frame to find the tissue in this category
    #   with the highest number of positive individuals 
    #   This will be what we use to count the pos\neg distribution
    
    tissue.df <- data.frame(tissue = character(),
                            n_pos = numeric(),
                            n_total = numeric())
    
    # Loop over all tissues
    for (tissue.ii in unique(dat.dpi$sample_rep)) {
      
      dat.tissue <- subset(dat.dpi, sample_rep == tissue.ii)
      
      n_pos <- sum(dat.tissue$pos_value) # num positive
      n_total <- nrow(dat.tissue) # num total individuals
      
      tissue.df <- rbind(tissue.df,
                         data.frame(tissue = tissue.ii,
                                    n_pos = n_pos,
                                    n_total = n_total))
    }
    
    tissue.df$n_pos <- as.numeric(tissue.df$n_pos)
    tissue.df$n_total <- as.numeric(tissue.df$n_total)
    
    # Order the tissues based on # positive and # tested. Proceed with 1st row.
    tissue.df <- tissue.df[order(tissue.df$n_pos, tissue.df$n_total, decreasing = TRUE), ]
    #print(tissue.df)
    
    # Set variables for maximum positives, maximum tested, and chosen tissue
    max_n_positive <- tissue.df$n_pos[1]
    max_n_tested <- tissue.df$n_total[1]
    max_tissue <- tissue.df$tissue[1]
    
    # Check whether we can add negative individuals
    #   This is only possible if there's only one tested tissue or if the sum of
    #   positive individuals from all other locations is less than the total
    #   number of tested individuals in our chosen max tissue
    
    total_n_positive <- sum(as.numeric(tissue.df$n_pos))
    total_n_tested <- max(tissue.df$n_total)
    
    if (nrow(tissue.df) == 1 | sum(tissue.df$n_pos[2:nrow(tissue.df)]) == 0) {
      max_n_negative <- max_n_tested - max_n_positive
    }
    else if (total_n_positive < total_n_tested) {
      max_n_negative <- total_n_tested - total_n_positive
    }
    else {
      max_n_negative <- 0
    }
    
    #cat(blue("Adding", max_n_positive, "positives and",
    #         max_n_negative, "negatives\n\n"))
    
    # Now create new data frame for adding these individuals
    tissue_row_number <- which(dat.dpi$sample_rep == max_tissue)[1]
    df.new <- dat.dpi[rep(tissue_row_number, max_n_positive + max_n_negative), ]
    df.new$pos_value <- c(rep(1, max_n_positive), rep(0, max_n_negative))
    df.new$first_lower_bound <- c(rep(0, max_n_positive), rep(dpi.ii, max_n_negative))
    df.new$first_upper_bound <- c(rep(dpi.ii, max_n_positive), rep(Inf, max_n_negative))

    dat.fill <- rbind(dat.fill, df.new)
    
  }
  return(dat.fill)
}

# Extracts peak time for no ID individuals
get_peak_times_noid <- function(dat) {
  
  dat <- subset(dat, sample_type == "Non-invasive")
  
  # Requirements: (1) quantitative titers
  #               (2) at least one positive
  #               (3) at least two sample times (including one positive)
  
  # Check these requirements:
  
  # Any positives?
  if (1 %in% unique(dat$pos_value)){has_positive <- 1} 
  else {has_positive <- 0}
  
  # Any quantitative information?
  if (FALSE %in% is.na(unique(as.numeric(dat$value)))) {has_quantitative <- 1}
  else {has_quantitative <- 0}
  
  # Are there enough sample times?
  if ((length(unique(dat$day_post_infection)) > 1)) {has_enough <- 1}
  else {has_enough <- 0}
  # REMOVED THIS CONDITION: | (TRUE %in% str_detect(dat$unit_subgrp, "Peak|peak"))
  # SHOULD ADD SPECIFICALLY BACK
  
  dat.new <- dat[0, ]
  
  # Proceed if all three are true
  if (sum(has_positive, has_quantitative, has_enough) == 3) {
    
    dat$day_post_infection <- as.numeric(dat$day_post_infection)
    
    # Find the day on which, across all individuals, the largest titer was observed
    max_observed_titer <- max(as.numeric(dat$value), na.rm = TRUE)
    day_of_max_titer <- dat$day_post_infection[!is.na(as.numeric(dat$value)) &
                                               as.numeric(dat$value) == max_observed_titer]
    
    #print(max_observed_titer)
    #print(day_of_max_titer)
    
    # Check if multiple days have this max titer (this is rare!)
    #     If so, take the earliest one
    if (length(unique(day_of_max_titer)) >= 2) {
      cat(blue("Note: multiple days (", paste0(unique(day_of_max_titer), collapse = ", "), 
               ") have the same max titer, selecting the earliest one.\n"))
      day_of_max_titer <- min(as.numeric(day_of_max_titer), na.rm = TRUE)
    }
    if (length(day_of_max_titer) > 1 & length(unique(day_of_max_titer)) == 1) {
      cat(blue("Note: the same day has two of the same peak values.\n"))
      day_of_max_titer <- unique(day_of_max_titer)
    }

    
    # Check whether multiple individuals have this peak on that day
    num_with_max_titer <- sum(dat$day_post_infection == day_of_max_titer &
                                dat$value == max_observed_titer)
    
    if (num_with_max_titer > 1) {
      cat("There are",  num_with_max_titer, "individuals at the peak on this day.\n")
    }
    
    #print(dat)
    #print(day_of_max_titer)
    #print(min(dat$day_post_infection))
    
    # Get the previous and next sample days to bound this peak time
    #   If the peak is the first sample day, then the lower bound is zero
    if (day_of_max_titer == min(dat$day_post_infection, na.rm = TRUE)) {
      day_before_max_titer <- 0
    }
    else {
      day_before_max_titer <- max(dat$day_post_infection[dat$day_post_infection < day_of_max_titer], na.rm = TRUE)
    }
    
    # If the peak is the last observed day, then the upper bound is Inf
    if (day_of_max_titer == max(dat$day_post_infection, na.rm = TRUE)) {
      day_after_max_titer <- Inf
    }
    else {
      day_after_max_titer <- min(dat$day_post_infection[dat$day_post_infection > day_of_max_titer], na.rm = TRUE)
    }
    
    # Create a dataframe with one row for this maximum peak
    dat.max <- dat[1, ]
    dat.max$peak_lower_bound <- day_before_max_titer
    dat.max$peak_observed_time <- day_of_max_titer
    dat.max$peak_observed_titer <- max_observed_titer
    dat.max$peak_upper_bound <- day_after_max_titer
    
    cat("The day of the max titer is", day_of_max_titer, "which has the bounds: [",
        day_before_max_titer, ",", day_after_max_titer, "].\n")
    
    
    # Now, check for whether we can extract any other peak times
    
    # Set all below LOD values to -5
    dat$value[dat$pos_value == 0] <- -5
    dat$value <- as.numeric(dat$value)
    
    # Start by checking viral loads on the same max titer day
    all_titers_on_max_day <- dat$value[dat$day_post_infection == day_of_max_titer &
                                       dat$pos_value == 1 &
                                       dat$value < max_observed_titer]
    
    # Set data frame to store
    dat.other <- dat[0, ]
    
    # If there is at least one other observation, continue with searching
    if (length(all_titers_on_max_day) >= 1){
      
      # Now create a data frame with the maximum observed viral load
      #   on every test day
      dpi.df <- data.frame(dpi = numeric(),
                           dpi_max_titer = numeric(),
                           peak_day = numeric())
      
      dpi_options <- as.numeric(unique(dat$day_post_infection))
      
      for (dpi.ii in sort(dpi_options)) {
        
        # Subset to the day in question
        dat.dpi <- subset(dat, as.numeric(day_post_infection) == dpi.ii)
        
        # Find the maximum load at this dpi - if all < LOD will be -5
        max_dpi_titer <- max(as.numeric(dat.dpi$value),  na.rm = TRUE)
        
        # Indicator for whether its the peak day
        if (day_of_max_titer != dpi.ii) {is_peak_dpi <- 0} else {is_peak_dpi <- 1}
        
        # Add to the data frame
        dpi.df <- rbind(dpi.df, 
                        data.frame(dpi = dpi.ii,
                                   dpi_max_titer = max_dpi_titer,
                                   peak_day = is_peak_dpi,
                                   n_total = nrow(dat.dpi)))
        
        
      }
      #print(dpi.df)
      
      # For all titers on the max day, check if they're bigger 
      #    than any previous max
      for (observed_titer in all_titers_on_max_day[!is.na(all_titers_on_max_day)]) {
        
        # Find the max on all previous or later days, adjusting for whether its
        #   the first or last sample day
        if (dpi.df$peak_day[1] == 1) {max_on_all_previous_days <- -5}
        else {max_on_all_previous_days <- max(dpi.df$dpi_max_titer[dpi.df$dpi < day_of_max_titer], na.rm = TRUE)}
        
        if (dpi.df$peak_day[nrow(dpi.df)] == 1){max_on_all_later_days <- -5}
        else{max_on_all_later_days <- max(dpi.df$dpi_max_titer[dpi.df$dpi > day_of_max_titer], na.rm = TRUE)}
        
        # Compare to the observed titer that we're looping over
        is_bigger_than_previous <- observed_titer >= max_on_all_previous_days
        is_bigger_than_later <- observed_titer >= max_on_all_later_days
        
        #print(observed_titer)
        #print(max_on_all_previous_days)
        #print(max_on_all_later_days)
        
        #print(is_bigger_than_previous)
        #print(is_bigger_than_later)
        
        
        # Check if it's bigger than the previous and max times and there's enough samples
        if (is_bigger_than_previous & is_bigger_than_later & nrow(dpi.df) > 1) {
          
          # If it's not the last sample day, we can add an interval censored
          #     time (note: time 0 counts as a sample day)
          if (dpi.df$peak_day[nrow(dpi.df)] != 1) {
            
            dat.current <- dat[1, ]
            dat.current$peak_lower_bound <- day_before_max_titer
            dat.current$peak_observed_time <- day_of_max_titer
            dat.current$peak_observed_titer <- observed_titer
            dat.current$peak_upper_bound <- day_after_max_titer
            
            dat.other <- rbind(dat.other, dat.current)
            
            #cat(observed_titer, "is bigger than all the other max titers,",
            #    "so we add another observations with bounds: [",
            #    day_before_max_titer, ",", day_after_max_titer, "]\n")
            
            
          }
          
          # If it's the last sample day, we can add another right censored time
          else if (dpi.df$peak_day[nrow(dpi.df)] == 1) {
            
            dat.current <- dat[1, ]
            dat.current$peak_lower_bound <- day_before_max_titer
            dat.current$peak_observed_time <- day_of_max_titer
            dat.current$peak_observed_titer <- observed_titer
            dat.current$peak_upper_bound <- Inf
            
            dat.other <- rbind(dat.other, dat.current)
            
            #cat(observed_titer, "is bigger than all the other max titers,",
            #    "so we add another observations with bounds: [",
            #    day_before_max_titer, ",", Inf, "]\n")
          }
          
        }
        
        # If its larger than the previous but not the later max values, can add
        #    some right censored times, unless at some point all individuals are negative
        #    afterwards
        else if (is_bigger_than_previous & !is_bigger_than_later & nrow(dpi.df) > 1) {
          
          # Check whether any later times are ever all negative
          ever_all_negative <- dpi.df$dpi_max_titer[dpi.df$dpi > day_of_max_titer] == -5
          
          # If so, add interval censored time
          if (TRUE %in% ever_all_negative) {
            
            first_all_neg_dpi <- min(dpi.df$dpi[dpi.df$dpi > day_of_max_titer &
                                                dpi.df$dpi_max_titer == -5], na.rm = TRUE)
            
            dat.current <- dat[1, ]
            dat.current$peak_lower_bound <- day_before_max_titer
            dat.current$peak_observed_time <- day_of_max_titer
            dat.current$peak_observed_titer <- observed_titer
            dat.current$peak_upper_bound <- first_all_neg_dpi
            
            dat.other <- rbind(dat.other, dat.current)
            
            #print(dpi.df)
            #cat("Adding an interval censored observation, with bounds: [",
            #    day_before_max_titer, ",", first_all_neg_dpi, "]\n")
          }
          
          # If there's never an all-negative day, but the peak day is not 
          #   the first sample day, then add a right censored time
          else if (dpi.df$peak_day[1] != 1) {
            
            dat.current <- dat[1, ]
            dat.current$peak_lower_bound <- day_before_max_titer
            dat.current$peak_observed_time <- day_of_max_titer
            dat.current$peak_observed_titer <- observed_titer
            dat.current$peak_upper_bound <- Inf
            
            dat.other <- rbind(dat.other, dat.current)
            
            #print(dpi.df)
            #cat("Adding a right censored observation, with bounds: [",
            #    day_before_max_titer, ",", first_all_neg_dpi, "]\n")
            
          }
          # We can't infer bounds if the peak time is the very first sampling time
          #    and there are only ever positives!
          
        }
        
        else if (!is_bigger_than_previous & is_bigger_than_later) {
          
          # If it's not the last possible sampling time, then we can get an interval
          #   censored time
          if (dpi.df$peak_day[nrow(dpi.df)] != 1) {
            
            # Find the row number of the peak, which I'll use to reference
            peak_row <- which(dpi.df$peak_day == 1)
            
            if (peak_row == 2) {
              previous_day <- 0
            }
            else {
              # This actually isn't happening in the data
              previous_day <- dpi.df$dpi[peak_row - 2]
            }
            
            dat.current <- dat[1, ]
            dat.current$peak_lower_bound <- previous_day
            dat.current$peak_observed_time <- day_of_max_titer
            dat.current$peak_observed_titer <- observed_titer
            dat.current$peak_upper_bound <- day_after_max_titer
            
            dat.other <- rbind(dat.other, dat.current)
            
            cat(blue("Adding an individual with bounds: [", 
                     previous_day, ",", day_after_max_titer, "]\n"))
            
          }
          # If the first sampling time is all below the LOD,
          #    then know it's right censored with that time
          else if (dpi.df$dpi_max_titer[1] == -5){
            
            # Hmmmm.... I think I'm going to skip this
            #    Basically just "any time after the first positive...
            #    which I add later anyway!
            
          }
          
        }
        
        # Can't infer any bounds if its not larger than any other times
        else if (!is_bigger_than_previous & !is_bigger_than_later) {
        
          }
      
      }
    }
    
    # Combine the dataframes together
    dat.new <- rbind(dat.max, dat.other)
  }
  
  return(dat.new)
}

# Combines the first positive and peak times for no ID individuals into one
#     consistent data frame
consolidate_times_noid <- function(dat.first, dat.peak) {
  
  # Excluding the right censored first positives (since they are by definition)
  #    not observed to be positive
  dat.first.right.cens <- subset(dat.first, !is.na(first_upper_bound) & first_upper_bound == Inf)
  dat.first.int.cens <- subset(dat.first, !is.na(first_upper_bound) & first_upper_bound != Inf)
  dat.peak <- subset(dat.peak, !is.na(peak_lower_bound))
  
  # Set bounds on interval censored individuals to the max possible range
  dat.first.int.cens$first_lower_bound <- min(dat.first.int.cens$first_lower_bound)
  dat.first.int.cens$first_upper_bound <- max(dat.first.int.cens$first_upper_bound)
  
  # For each row, when possible, couple the first positive range with an observed
  #     peak time range. 
  for (row_num in 1:nrow(dat.first.int.cens)) {
    # If the row exists in the peak dataframe, use those values
    if (!is.na(dat.peak$peak_lower_bound[row_num])) {
      dat.first.int.cens$peak_lower_bound[row_num] <- dat.peak$peak_lower_bound[row_num]
      dat.first.int.cens$peak_observed_time[row_num] <- dat.peak$peak_observed_time[row_num]
      dat.first.int.cens$peak_observed_titer[row_num] <- dat.peak$peak_observed_titer[row_num]
      dat.first.int.cens$peak_upper_bound[row_num] <- dat.peak$peak_upper_bound[row_num]
    }
    else {
      dat.first.int.cens$peak_lower_bound[row_num] <- dat.first.int.cens$first_lower_bound[row_num]
      dat.first.int.cens$peak_upper_bound[row_num] <- Inf
    }
  }
  
  dat.new <- rbind(dat.first.int.cens, dat.first.right.cens)
  
  colnames(dat.new)[colnames(dat.new) == "assay_type"] <- "rna_type"
  
  return(dat.new)

}

# Extracts bounds for peak and last positive times, when other observations allow
fill_in_event_times <- function(dat) {
  
  dat$ever_positive <- NA
  dat$has_peak <- NA
  dat$has_last_positive <- NA
  
  # Set FLAGs to NA
  dat$first_lower_bound[dat$first_lower_bound == "FLAG"] <- NA
  dat$first_upper_bound[dat$first_upper_bound == "FLAG"] <- NA
  dat$peak_lower_bound[dat$peak_lower_bound == "FLAG"] <- NA
  dat$peak_upper_bound[dat$peak_upper_bound == "FLAG"] <- NA
  dat$last_lower_bound[dat$last_lower_bound == "FLAG"] <- NA
  dat$last_upper_bound[dat$last_upper_bound == "FLAG"] <- NA
  
  for (row_num in 1:nrow(dat)) {
    first_lower <- dat$first_lower_bound[row_num]
    first_upper <- dat$first_upper_bound[row_num]
    peak_lower <- dat$peak_lower_bound[row_num]
    peak_observed <- dat$peak_observed_time[row_num]
    peak_upper <- dat$peak_upper_bound[row_num]
    last_lower <- dat$last_lower_bound[row_num]
    last_upper <- dat$last_upper_bound[row_num]
    
    # Check if ever positive
    if (!is.na(first_upper) & first_upper != Inf) {ever_positive <- 1}
    else {ever_positive <- 0}
    
    # Now check if peak time can be bounded based on first and last positive times
    if (ever_positive == 1 & is.na(peak_lower) & !is.na(last_lower)) {
      peak_lower <- first_lower
      peak_upper <- last_upper
    }
    # If there's no last positive time, can at least bound peak time
    else if (ever_positive == 1 & is.na(peak_lower) & is.na(last_lower)) {
      peak_lower <- first_lower
      peak_upper <- Inf
    }
    
    # Now check if there's no last positive time for individuals with
    #   peak times (last positive has to be later than observed time for peak)
    if (ever_positive == 1 & !is.na(peak_lower) & 
        is.na(last_lower) & !is.na(peak_observed)) {
      last_lower <- peak_observed
      last_upper <- Inf
    }
    else if (ever_positive == 1 & !is.na(peak_lower) & 
        is.na(last_lower) & is.na(peak_observed)) {
      last_lower <- peak_lower
      last_upper <- Inf
    }
    
    # Now add indicator for whether has peak info
    if (ever_positive == 1 & !is.na(peak_lower)) {has_peak <- 1}
    else {has_peak <- 0}
    
    # And add indicator for whether has last positive info
    if (ever_positive == 1 & !is.na(last_lower)) {has_last_positive <- 1}
    else {has_last_positive <- 0}
    
    # Write true entries
    dat$first_lower_bound[row_num] <- first_lower
    dat$first_upper_bound[row_num] <- first_upper
    dat$peak_lower_bound[row_num] <- peak_lower 
    dat$peak_upper_bound[row_num] <- peak_upper 
    dat$last_lower_bound[row_num] <- last_lower 
    dat$last_upper_bound[row_num] <- last_upper
    
    dat$ever_positive[row_num] <- ever_positive
    dat$has_peak[row_num] <- has_peak
    dat$has_last_positive[row_num] <- has_last_positive
    
  }
  
  # Check to make sure all individuals with peak time have a last positive time
  if (0 %in% unique(dat$has_last_positive[dat$has_peak == 1])) {
    cat("WARNING: Some individuals have bounded peak times",
        "but not last positive times.\n")
  }
  
  # Remove rows with no data
  dat <- subset(dat, !(is.na(first_lower_bound) & 
                         is.na(peak_lower_bound) & 
                         is.na(last_lower_bound)))
  
  return(dat)
}


assign_location_names <- function(df, organ_group = FALSE){
  
  if ("organ_group" %notin% colnames(df)) {
    df$organ_group <- df$location_grp
  }
  
  df$location_name <- NA
  
  df$location_name[df$location_idx == 0 & df$organ_group == "URT"] <- "Nose"
  df$location_name[df$location_idx == 1 & df$organ_group == "URT"] <- "Throat"
  df$location_name[df$location_idx == -9999 & df$organ_group == "URT"] <- "Nose"
  
  df$location_name[df$location_idx == 0 & df$organ_group == "LRT"] <- "Trachea"
  df$location_name[df$location_idx == 1 & df$organ_group == "LRT"] <- "Lung"
  df$location_name[df$location_idx == -9999 & df$organ_group == "LRT"] <- "Lung"
  
  df$location_name[df$location_idx == 0 & df$organ_group == "GI"] <- "Upper GI"
  df$location_name[df$location_idx == 1 & df$organ_group == "GI"] <- "Lower GI"
  df$location_name[df$location_idx == -9999 & df$organ_group == "GI"] <- "Lower GI"
  
  df$location_name <- factor(df$location_name,
                             levels = c("Nose", 
                                        "Throat",
                                        "Trachea",
                                        "Lung",
                                        "Upper GI",
                                        "Lower GI"))
  
  return(df)
}

assign_route_names <- function(df) {
  
  df$route_name <- NA
  
  df$route_name[df$route_idx == 1] <- "IN"
  df$route_name[df$route_idx == 2] <- "IT"
  df$route_name[df$route_idx == 3] <- "IN + IT"
  df$route_name[df$route_idx == 4] <- "AE"
  df$route_name[df$route_idx == 5] <- "IG"
  
  df$route_name <- factor(df$route_name, levels = c("IN",
                                                    "IT",
                                                    "IN + IT",
                                                    "AE",
                                                    "IG"))
  
  return(df)
  
}

assign_tissue_names <- function(df) {
  
  df$tissue_name <- NA
  
  df$tissue_name[df$tissue_idx == 1] <- "Nose"
  df$tissue_name[df$tissue_idx == 2] <- "Throat"
  df$tissue_name[df$tissue_idx == 3] <- "Trachea"
  df$tissue_name[df$tissue_idx == 4] <- "Lung"
  df$tissue_name[df$tissue_idx == 5] <- "Upper GI"
  df$tissue_name[df$tissue_idx == 6] <- "Lower GI"
  
  df$tissue_name <- factor(df$tissue_name, levels = c("Nose",
                                                      "Throat",
                                                      "Trachea",
                                                      "Lung",
                                                      "Upper GI",
                                                      "Lower GI"))
  
  return(df)
  
}

assign_assay_names <- function(df, long = TRUE) {
  
  if (long == FALSE) {
    df$assay_name <- NA
    df$assay_name[df$assay_idx == 1] <- "Total RNA"
    df$assay_name[df$assay_idx == 2] <- "gRNA"
    df$assay_name[df$assay_idx == 3] <- "sgRNA"
    df$assay_name[df$assay_idx == 4] <- "Culture"
    df$assay_name[df$assay_idx == -9999] <- "Unknown"
    df$assay_name <- factor(df$assay_name,
                            levels = c("Total RNA", "gRNA", "sgRNA", "Culture", "Unknown")) 
  }
  else {
    #df$assay_name <- NA
    #df$assay_name[df$assay_idx == 1] <- "T↑" 
    #df$assay_name[df$assay_idx == 2] <- "T↓" 
    #df$assay_name[df$assay_idx == 3] <- "SG↑"
    #df$assay_name[df$assay_idx == 4] <- "SG↓"
    #df$assay_name[df$assay_idx == 5] <- "C↑"
    #df$assay_name[df$assay_idx == 6] <- "C↓"
    #df$assay_name <- factor(df$assay_name,
    #                        levels = c("T↑" ,"T↓",
    #                                   "SG↑", "SG↓",
    #                                   "C↑", "C↓")) 
    df$assay_name <- NA
    df$assay_name[df$assay_idx == 1] <- "Total RNA" 
    df$assay_name[df$assay_idx == 2] <- "gRNA" 
    df$assay_name[df$assay_idx == 3] <- "sgRNA" 
    df$assay_name[df$assay_idx == 4] <- "Culture"
    df$assay_name <- factor(df$assay_name,
                            levels = c("Total RNA",
                                       "gRNA",
                                       "sgRNA",
                                       "Culture")) 
  }
  
  return(df)
}

assign_species_names <- function(df) {
  
  df$species_name <- NA
  df$species_name[df$sp_idx == 1] <- "Rhesus macaque"
  df$species_name[df$sp_idx == 2] <- "Cynomolgus macaque"
  df$species_name[df$sp_idx == 3] <- "African green monkey"
  df$species_name <- factor(df$species_name,
                            levels = c("Rhesus macaque", 
                                       "Cynomolgus macaque",
                                       "African green monkey"))
  
  return(df)
}

assign_age_names <- function(df) {
  
  df$age_name <- NA
  df$age_name[df$age_idx == 1] <- "Juvenile"
  df$age_name[df$age_idx == 2] <- "Adult"
  df$age_name[df$age_idx == 3] <- "Geriatric"
  df$age_name[df$age_idx == -9999] <- "Unknown"
  df$age_name <- factor(df$age_name,
                        levels = c("Juvenile", 
                                   "Adult",
                                   "Geriatric",
                                   "Unknown"))
  
  return(df)
}

assign_sex_names <- function(df) {
  
  df$sex_name <- NA
  df$sex_name[df$sex_idx == 0] <- "Female"
  df$sex_name[df$sex_idx == 1] <- "Male"
  df$sex_name[df$sex_idx == -9999] <- "Unknown"
  df$sex_name <- factor(df$sex_name,
                        levels = c("Female", 
                                   "Male",
                                   "Unknown"))
  
  return(df)
}

assign_all_names <- function(df) {
  df <- assign_location_names(df)
  df <- assign_route_names(df)
  df <- assign_assay_names(df)
  df <- assign_species_names(df)
  df <- assign_age_names(df)
  df <- assign_sex_names(df)
  
  if ("tissue_idx" %notin% colnames(df)) {
    df$tissue_idx <- as.numeric(df$location_name)
  }
  
  df <- assign_tissue_names(df)
  return(df)
}

assign_metric_names_units <- function(df) {
  
  df$metric_name <- NA
  
  df$metric_name[df$metric == "diff_percent"] <- "Probability of positivity (%)"
  df$metric_name[df$metric == "diff_first"] <- "Time to detectability (days)"
  df$metric_name[df$metric == "diff_peak"] <- "Time to peak titer (days)"
  df$metric_name[df$metric == "diff_titer"] <- "Peak titer (log10 pfu)"
  df$metric_name[df$metric == "diff_last"] <- "Time to undetectability (days)"
  df$metric_name[df$metric == "diff_auc"] <- "AUC"
  df$metric_name[df$metric == "diff_duration"] <- "Duration of infection (days)"
  
  df$metric_name <- factor(df$metric_name, levels = rev(c("Probability of positivity (%)",
                                                          "Time to detectability (days)",
                                                          "Time to peak titer (days)",
                                                          "Peak titer (log10 pfu)",
                                                          "Time to undetectability (days)",
                                                          "Duration of infection (days)",
                                                          "AUC")))
  
  return(df)
}

prepare_stan_data_joint <- function(dat) {
  
  # Subset to target location & exclude individuals without dose info
  dat <- subset(dat, !is.na(dose_nose) & location_grp %in% c("URT", "LRT", "GI"))
  
  # Set upper GI to have no peak or last info
  dat$has_peak[dat$location_grp == "GI" & dat$location_idx == 0] <- 0
  dat$peak_lower_bound[dat$location_grp == "GI" & dat$location_idx == 0] <- NA
  dat$peak_upper_bound[dat$location_grp == "GI" & dat$location_idx == 0] <- NA
  dat$has_last_positive[dat$location_grp == "GI" & dat$location_idx == 0] <- 0
  dat$last_lower_bound[dat$location_grp == "GI" & dat$location_idx == 0] <- NA
  dat$last_upper_bound[dat$location_grp == "GI" & dat$location_idx == 0] <- NA
  
  # NEW: If lower and upper bounds are 0 and infinity, set both to NA
  #.         and set has_XX to 0.
  #dat$peak_lower_bound[dat$peak_lower_bound == 0 & dat$peak_upper_bound == Inf] <- NA
  #dat$peak_upper_bound[is.na(dat$peak_lower_bound) & dat$peak_upper_bound == Inf] <- NA
  #dat$last_lower_bound[dat$last_lower_bound == 0 & dat$last_upper_bound == Inf] <- NA
  #dat$last_upper_bound[is.na(dat$last_lower_bound) & dat$last_upper_bound == Inf] <- NA
  #dat$has_peak[is.na(dat$peak_lower_bound) & is.na(dat$peak_upper_bound)] <- 0
  #dat$has_last_positive[is.na(dat$last_lower_bound) & is.na(dat$last_upper_bound)] <- 0
  
  # Set NA lower and upper bounds to values for stan
  dat$peak_lower_bound[is.na(dat$peak_lower_bound)] <- 9999
  dat$peak_upper_bound[is.na(dat$peak_upper_bound)] <- 10000
  dat$last_lower_bound[is.na(dat$last_lower_bound)] <- 9999
  dat$last_upper_bound[is.na(dat$last_upper_bound)] <- 10000
  
  # Classify last positive censor type
  dat$last_positive_type[dat$last_upper_bound == Inf] <- 0
  dat$last_positive_type[dat$last_upper_bound != Inf] <- 1
  
  # Set right censored upper bounds to something besides infinity
  dat$first_upper_bound[dat$first_upper_bound == Inf] <- 45 #50
  dat$peak_upper_bound[dat$peak_upper_bound == Inf] <- 50 #60
  dat$last_upper_bound[dat$last_upper_bound == Inf] <- 55 #70
  
  # Change value for samples without titer
  dat$peak_observed_titer[dat$has_titer == 0] <- 0
  
  # Update lab groups based on how many are available per study
  labs_with_titer <- unique(dat$lab_group[dat$has_titer == 1])
  labs_without_titer <- unique(dat$lab_group[dat$lab_group %notin% labs_with_titer])
  
  dat$lab_factor <- factor(dat$lab_group, levels = c(labs_with_titer, labs_without_titer))
  dat$lab_idx <- as.numeric(dat$lab_factor)
  
  # Set organ location numbers
  dat$organ_location[dat$location_grp == "URT"] <- 1
  dat$organ_location[dat$location_grp == "LRT"] <- 2
  dat$organ_location[dat$location_grp == "GI"] <- 3
  
  # Set tissue location numbers
  dat$tissue_location[dat$location_grp == "URT" & dat$location_idx == 0] <- 1
  dat$tissue_location[dat$location_grp == "URT" & dat$location_idx == 1] <- 2
  dat$tissue_location[dat$location_grp == "LRT" & dat$location_idx == 0] <- 3
  dat$tissue_location[dat$location_grp == "LRT" & dat$location_idx == 1] <- 4
  dat$tissue_location[dat$location_grp == "GI" & dat$location_idx == 0] <- 5
  dat$tissue_location[dat$location_grp == "GI" & dat$location_idx == 1] <- 6
  dat$tissue_location[dat$location_idx == -9999] <- -9999
  
  # Fix AE inoculation values
  dat$inoc_dose_AE_pfu[is.na(dat$inoc_dose_AE_pfu)] <- 0

  # Prep data for Stan
  dat.stan <- list(N = nrow(dat),
                   L_assay = max(dat$assay_idx),
                   L_lab = max(dat$lab_idx),
                   L_lab_titer = max(dat$lab_idx[dat$has_titer == 1]),
                   ever_positive = dat$ever_pos,
                   has_peak = dat$has_peak,
                   has_titer = dat$has_titer,
                   has_last_positive = dat$has_last_positive,
                   last_positive_type = dat$last_positive_type,
                   organ_location = dat$organ_location,
                   tissue_location = dat$tissue_location,
                   location = dat$location_idx,
                   inoc = dat$inoc_idx,
                   route = dat$route_idx,
                   dose_nose = dat$dose_nose,
                   dose_throat = dat$dose_throat,
                   dose_trachea = dat$dose_trachea,
                   dose_lung = dat$dose_lung,
                   dose_gi = dat$dose_gi,
                   dose_total = dat$inoc_dose_total_pfu,
                   dose_IN = 0, #dat$dose_IN,
                   dose_OR = 0, #dat$dose_OR,
                   dose_IT = 0, #dat$dose_IT,
                   dose_IB = 0, #dat$dose_IB,
                   dose_IG = 0, #dat$dose_IG,
                   dose_AE = dat$inoc_dose_AE_pfu,
                   dose_OC = 0, #dat$dose_OC,
                   assay = dat$assay_idx,
                   assay_full = dat$assay_full,
                   lab = dat$lab_idx,
                   sex = dat$sex_idx,
                   age = dat$age_idx,
                   species = dat$sp_idx,
                   prob_male = dat$prob_male,
                   prob_female = dat$prob_female,
                   prob_geriatric = dat$prob_geriatric,
                   first_lower_bounds = as.numeric(dat$first_lower_bound),
                   first_upper_bounds = as.numeric(dat$first_upper_bound),
                   peak_lower_bounds = as.numeric(dat$peak_lower_bound),
                   peak_observed_time = as.numeric(dat$peak_observed_time), # NOT USED FOR FITTING
                   peak_upper_bounds = as.numeric(dat$peak_upper_bound),
                   last_lower_bounds = as.numeric(dat$last_lower_bound),
                   last_upper_bounds = as.numeric(dat$last_upper_bound),
                   peak_observed_titer = as.numeric(dat$peak_observed_titer),
                   indiv_idx = dat$indiv_idx)
  
  return(dat.stan)
  
}


# Assigns doses based on our standard assumptions about distribution for the
#     5 route categories we studied
assign_dose_distribution <- function(total_dose) {
  
  total_dose <- 10^(total_dose)
  
  # URT exposure; 50% nose, 50% throat
  route1.dose.nose <<- log10(total_dose * 1/2)
  route1.dose.throat <<- log10(total_dose * 1/2)
  route1.dose.trachea <<- 0
  route1.dose.lung <<- 0
  route1.dose.gi <<- 0
  
  # LRT exposure; 100% trachea
  route2.dose.nose <<- 0
  route2.dose.throat <<- 0
  route2.dose.trachea <<- log10(total_dose)
  route2.dose.lung <<- 0
  route2.dose.gi <<- 0
  
  # URT + LRT (Liquid) exposure; 25% nose, 25% throat, 50% trachea
  route3.dose.nose <<- log10(total_dose * 1/4)
  route3.dose.throat <<- log10(total_dose * 1/4)
  route3.dose.trachea <<- log10(total_dose * 1/2)
  route3.dose.lung <<- 0
  route3.dose.gi <<- 0
  
  # URT + LRT (Aerosol) exposure; 25% across nose, throat, trachea, lung
  route4.dose.nose <<- log10(total_dose * 1/4)
  route4.dose.throat <<- log10(total_dose * 1/4)
  route4.dose.trachea <<- log10(total_dose * 1/4)
  route4.dose.lung <<- log10(total_dose * 1/4)
  route4.dose.gi <<- 0
  
  # GI exposure; 100% GI
  route5.dose.nose <<- 0
  route5.dose.throat <<- 0
  route5.dose.trachea <<- 0
  route5.dose.lung <<- 0
  route5.dose.gi <<- log10(total_dose)
  
}


get_model_estimates_joint <- function(fit, 
                                      draws = "all",
                                      n_draws = NA,
                                      dose_options = c(3, 5, 7),
                                      sp_options = 1,
                                      sp_standard = 1,
                                      age_options = 2,
                                      age_standard = 2,
                                      sex_options = 0,
                                      sex_standard = 0,
                                      assay_options = 1, 
                                      lab_standard = 1,
                                      lab_effect = "No") {
  
  # Extract all parameters
  pars <- c("percent_intercept", "percent_location", #"percent_inoc",
            "percent_route", "percent_dose",
            "percent_sex", "percent_age", "percent_species",
            "percent_assay", "percent_lab",
            
            "shape_intercept_first",
            
            "median_intercept_first", "median_location_first",
            "median_route_first", "median_dose_first", 
            "median_sex_first", "median_age_first", "median_species_first",
            "median_assay_first", "median_lab_first", 
            
            "shape_intercept_peak",
            
            "median_intercept_peak", "median_location_peak",
            "median_route_peak", "median_dose_peak", 
            "median_sex_peak", "median_age_peak", "median_species_peak",
            "median_assay_peak", "median_lab_peak", 
            "median_time_peak",
            
            "titer_intercept", "titer_location", "titer_inoc",
            "titer_route", "titer_dose",
            "titer_sex", "titer_age", "titer_species",
            "titer_assay", "titer_lab", 
            "titer_time",
            "titer_intercept_sd", 
            
            "shape_intercept_last",
            
            "median_intercept_last", "median_location_last",
            "median_route_last", "median_dose_last", 
            "median_sex_last", "median_age_last", "median_species_last",
            "median_assay_last", "median_lab_last", 
            "median_time_last"
  )
  
  
  # Need wide format samples to get correlated values
  sample.pars <- fit$draws(pars, format = "df")
  
  p.fits <- data.frame(sample_num = NA,
                       organ_group = NA,
                       organ_idx = NA,
                       dose_nose = NA,
                       dose_throat = NA,
                       dose_trachea = NA,
                       dose_lung = NA,
                       dose_gi = NA,
                       dose_total = NA,
                       route_idx = NA,
                       sp_idx = NA,
                       age_idx = NA,
                       sex_idx = NA,
                       assay_idx = NA,
                       location_idx = NA,
                       lab_idx = NA,
                       percent_positive = NA,
                       first_pos_scale = NA,
                       first_pos_shape = NA,
                       first_pos_median = NA,
                       first_pos_sample = NA,
                       peak_scale = NA,
                       peak_shape = NA,
                       peak_median = NA,
                       peak_sample = NA,
                       titer_mean = NA,
                       titer_sample = NA,
                       last_scale = NA,
                       last_shape = NA,
                       last_median = NA,
                       last_sample = NA
  )
  
  # Line num to loop through vectors of randomly sampled cofactors
  line_num <- 1
  
  # Get row numbers to loop through
  if (draws == "random") {
    sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  }
  else if (draws == "all") {
    sample_num <- 1:nrow(sample.pars)
    cat("Using all parameter samples to generate estimates -- this may take some time!\n\n")
  }
  
  for (total.dose.ii in dose_options) {
    
    total_dose <- total.dose.ii
    assign_dose_distribution(total_dose)
    
    cat("\nGenerating estimates for total dose", total.dose.ii, "\n")
    
    for (route.ii in 1:5) {
      
      cat("Generating estimates for route", route.ii, "\n")
      
      dose_nose <- get(paste0("route", route.ii, ".dose.nose"))
      dose_throat <- get(paste0("route", route.ii, ".dose.throat"))
      dose_trachea <- get(paste0("route", route.ii, ".dose.trachea"))
      dose_lung <- get(paste0("route", route.ii, ".dose.lung"))
      dose_gi <- get(paste0("route", route.ii, ".dose.gi"))
      
      for (tissue.ii in 1:6) {
        
        # Set organ groups and location_idx
        if (tissue.ii == 1) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 0
        }
        else if (tissue.ii == 2) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 1
        }
        else if (tissue.ii == 3) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 0
        }
        else if (tissue.ii == 4) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 1
        }
        else if (tissue.ii == 5) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 0
        }
        else if (tissue.ii == 6) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 1
        }
        
        # Set inoculated category
        if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0){inoc <- 1}
        if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1){inoc <- 1}
        if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0){inoc <- 1}
        if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1){inoc <- 1}
        if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0){inoc <- 0}
        
        for (assay.ii in assay_options) {
          
          for (sex.ii in sex_options) {
            
            for (age.ii in age_options) {
              
              for (sp.ii in sp_options) {
                
                for (ii in sample_num) {
                  
                  # First positive  -------------------------
                  
                  ### Cure fraction -------
                  
                  percent_intercept_column <-paste0("percent_intercept[", organ.ii, "]") 
                  percent_intercept <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_intercept_column)])
                  
                  percent_location_column <-paste0("percent_location[", organ.ii, "]") 
                  percent_location <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_location_column)])
                  
                  #percent_inoc_column <-paste0("percent_inoc[", organ.ii, "]") 
                  #percent_inoc <- as.numeric(
                  #  sample.pars[ii, which(colnames(sample.pars) == percent_inoc_column)])
                  
                  percent_route_column <- paste0("percent_route[", route.ii, ",", organ.ii, "]") 
                  percent_route <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_route_column)])
                  
                  percent_dose_nose_column <- paste0("percent_dose[1,", tissue.ii, "]") 
                  percent_dose_nose <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_nose_column)])
                  
                  percent_dose_throat_column <- paste0("percent_dose[2,", tissue.ii, "]") 
                  percent_dose_throat <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_throat_column)])
                  
                  percent_dose_trachea_column <- paste0("percent_dose[3,", tissue.ii, "]") 
                  percent_dose_trachea <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_trachea_column)])
                  
                  percent_dose_lung_column <- paste0("percent_dose[4,", tissue.ii, "]") 
                  percent_dose_lung <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_lung_column)])
                  
                  percent_dose_gi_column <- paste0("percent_dose[5,", tissue.ii, "]") 
                  percent_dose_gi <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_gi_column)])
                  
                  percent_sex_column <- paste0("percent_sex[", organ.ii, "]") 
                  percent_sex <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_sex_column)])
                  
                  percent_age_column <- paste0("percent_age[", age.ii, ",", organ.ii, "]")
                  percent_age <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_age_column)])
                  
                  percent_species_column <- paste0("percent_species[", sp.ii, ",", organ.ii, "]")
                  percent_species <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_species_column)])
                  
                  percent_assay_column <- paste0("percent_assay[", assay.ii, ",", organ.ii, "]")
                  percent_assay <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_assay_column)])
                  
                  percent_lab_column <- paste0("percent_lab[", lab_standard, ",", organ.ii,  "]")
                  percent_lab <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_lab_column)])
                  
                  if (lab_effect == "No") {
                    percent_lab <- 0
                  }
                  
                  #print(percent_sex)
                  #print(percent_age)
                  #print(percent_species)
                  
                  percent_trans <- percent_intercept + 
                    #percent_inoc * inoc + 
                    percent_location * loc.ii +
                    percent_dose_nose * dose_nose +
                    percent_dose_throat * dose_throat +
                    percent_dose_trachea * dose_trachea +
                    percent_dose_lung * dose_lung +
                    percent_dose_gi * dose_gi +
                    percent_route + 
                    percent_age +
                    percent_species +
                    percent_sex * sex.ii +
                    percent_assay + 
                    percent_lab
                  
                  percent <- exp(percent_trans) / (1 + exp(percent_trans))
                  
                  
                  ## Get shape parameter ------
                  
                  shape_first_column <- paste0("shape_intercept_first[", organ.ii, "]")
                  #shape_first_column <- paste0("shape_intercept_first[", assay.ii, ",", organ.ii, "]")
                  shape_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_first_column)])
                  
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_first[", organ.ii, "]")
                  median_intercept_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                 
                  median_location_column <-paste0("median_location_first[", organ.ii, "]") 
                  median_location_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_first[", route.ii, ",", organ.ii, "]") 
                  median_route_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
              
                  median_dose_nose_column <- paste0("median_dose_first[1,", tissue.ii, "]") 
                  median_dose_nose_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_first[2,", tissue.ii, "]") 
                  median_dose_throat_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_first[3,", tissue.ii, "]") 
                  median_dose_trachea_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_first[4,", tissue.ii, "]") 
                  median_dose_lung_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_first[5,", tissue.ii, "]") 
                  median_dose_gi_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_first[", organ.ii, "]") 
                  median_sex_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_first[", age.ii, ",", organ.ii, "]")
                  median_age_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_first[", sp.ii, ",", organ.ii,  "]")
                  median_species_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_first[", assay.ii, ",", organ.ii, "]")
                  median_assay_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_first[", lab_standard, ",", organ.ii,  "]")
                  median_lab_first <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_first <- 0
                  }
                  
                  median_first <- exp(-1/10 * (median_intercept_first + 
                                                 median_location_first * loc.ii +  
                                                 median_route_first +  
                                                 median_dose_nose_first * dose_nose +  
                                                 median_dose_throat_first * dose_throat + 
                                                 median_dose_trachea_first * dose_trachea + 
                                                 median_dose_lung_first * dose_lung + 
                                                 median_dose_gi_first * dose_gi + 
                                                 median_species_first + 
                                                 median_age_first +  
                                                 median_sex_first * sex.ii +  
                                                 median_assay_first + 
                                                 median_lab_first
                  ))
                  
                  scale_first <- median_first / (log(2)^(1 / shape_first))
                  
                  ## Make calculation ------
                  
                  #first_pos_quantile <- scale_first * ((-log(1 - prob_level)) ^ (1/shape_first))
                  first_pos_sample <- rweibull(1, shape = shape_first, scale = scale_first)
                  
                  
                  # Peak time --------------------------------------------------
                  
                  ## Get shape parameter ------
                  
                  shape_peak_column <- paste0("shape_intercept_peak[", organ.ii, "]")
                  #shape_peak_column <- paste0("shape_intercept_peak[", assay.ii, ",", organ.ii, "]")
                  shape_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_peak_column)])
                  
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_peak[", organ.ii, "]")
                  median_intercept_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                  
                  median_location_column <-paste0("median_location_peak[", organ.ii, "]") 
                  median_location_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_peak[", route.ii, ",", organ.ii, "]") 
                  median_route_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
                  
                  median_dose_nose_column <- paste0("median_dose_peak[1,", tissue.ii, "]") 
                  median_dose_nose_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_peak[2,", tissue.ii, "]") 
                  median_dose_throat_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_peak[3,", tissue.ii, "]") 
                  median_dose_trachea_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_peak[4,", tissue.ii, "]") 
                  median_dose_lung_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_peak[5,", tissue.ii, "]") 
                  median_dose_gi_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_peak[", organ.ii, "]") 
                  median_sex_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_peak[", age.ii, ",", organ.ii, "]")
                  median_age_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_peak[", sp.ii, ",", organ.ii,  "]")
                  median_species_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_peak[", assay.ii, ",", organ.ii, "]")
                  median_assay_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_peak[", lab_standard, ",", organ.ii,  "]")
                  median_lab_peak <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  median_time_column <- paste0("median_time_peak[", inoc + 1, ",", organ.ii, "]")
                  median_time_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_peak <- 0
                  }
                  
                  median_peak <- exp(-1/10 * (median_intercept_peak + 
                                                median_location_peak * loc.ii +  
                                                median_route_peak +  
                                                median_dose_nose_peak * dose_nose +  
                                                median_dose_throat_peak * dose_throat + 
                                                median_dose_trachea_peak * dose_trachea + 
                                                median_dose_lung_peak * dose_lung + 
                                                median_dose_gi_peak * dose_gi + 
                                                median_species_peak + 
                                                median_age_peak +  
                                                median_sex_peak * sex.ii +  
                                                median_assay_peak + 
                                                median_lab_peak +
                                                median_time_peak * median_first
                  ))
                  
                  scale_peak <- median_peak / (log(2)^(1 / shape_peak))
                  
                  ## Make calculation ------
                  
                  #peak_quantile <- scale_peak * ((-log(1 - prob_level)) ^ (1/shape_peak))
                  
                  
                  ## Draw sample based on first sample time --------------------
                  
                  median_peak_samp <- exp(-1/10 * (median_intercept_peak + 
                                                   median_location_peak * loc.ii +  
                                                   median_route_peak +  
                                                   median_dose_nose_peak * dose_nose +  
                                                   median_dose_throat_peak * dose_throat + 
                                                   median_dose_trachea_peak * dose_trachea + 
                                                   median_dose_lung_peak * dose_lung + 
                                                   median_dose_gi_peak * dose_gi + 
                                                   median_species_peak + 
                                                   median_age_peak +  
                                                   median_sex_peak * sex.ii +  
                                                   median_assay_peak + 
                                                   median_lab_peak +
                                                   median_time_peak * first_pos_sample
                  ))
                  
                  scale_peak_samp <- median_peak_samp / (log(2)^(1 / shape_peak))
                  
                  peak_sample <- rweibull(1, shape = shape_peak, scale = scale_peak_samp)
                  
                  
                  # Peak titer -------------------------------------------------
                  
                  ## Get parameters --------------------------------------------
                  
                  titer_intercept_column <-paste0("titer_intercept[", organ.ii, "]") 
                  titer_intercept <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_intercept_column)])
                  
                  titer_location_column <-paste0("titer_location[", organ.ii, "]") 
                  titer_location <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_location_column)])
                  
                  titer_inoc_column <-paste0("titer_inoc[", organ.ii, "]") 
                  titer_inoc <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_inoc_column)])
                  
                  titer_route_column <- paste0("titer_route[", route.ii, ",", organ.ii, "]") 
                  titer_route <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_route_column)])
                  
                  titer_dose_nose_column <- paste0("titer_dose[1,", organ.ii, "]") 
                  titer_dose_nose <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_nose_column)])
                  
                  titer_dose_throat_column <- paste0("titer_dose[2,", organ.ii, "]") 
                  titer_dose_throat <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_throat_column)])
                  
                  titer_dose_trachea_column <- paste0("titer_dose[3,", organ.ii, "]") 
                  titer_dose_trachea <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_trachea_column)])
                  
                  titer_dose_lung_column <- paste0("titer_dose[4,", organ.ii, "]") 
                  titer_dose_lung <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_lung_column)])
                  
                  titer_dose_gi_column <- paste0("titer_dose[5,", organ.ii, "]") 
                  titer_dose_gi <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_gi_column)])
                  
                  titer_sex_column <- paste0("titer_sex[", organ.ii, "]")
                  titer_sex <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_sex_column)])
                  
                  titer_age_column <- paste0("titer_age[", age.ii, ",", organ.ii, "]")
                  titer_age <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_age_column)])
                  
                  titer_species_column <- paste0("titer_species[", sp.ii, ",", organ.ii, "]")
                  titer_species <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_species_column)])
                  
                  titer_assay_column <- paste0("titer_assay[", assay.ii, ",", organ.ii, "]")
                  titer_assay <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_assay_column)])
                  
                  titer_lab_column <- paste0("titer_lab[", lab_standard, ",", organ.ii,  "]")
                  titer_lab <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_lab_column)])
                  
                  titer_time_column <- paste0("titer_time[", inoc + 1, ",", organ.ii, "]")
                  titer_time <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_time_column)])
                  
                  titer_sd_column <- paste0("titer_intercept_sd[", organ.ii, "]")
                  titer_sd <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_sd_column)])
                  
                  if (lab_effect == "No") {
                    titer_lab <- 0
                  }
                  
                  
                  ## Make calculation ------------------------------------------
                  
                  titer_mean <- titer_intercept +
                    titer_location * loc.ii + 
                    titer_inoc * inoc + 
                    titer_route +
                    titer_dose_nose * dose_nose +
                    titer_dose_throat * dose_throat +
                    titer_dose_trachea * dose_trachea +
                    titer_dose_lung * dose_lung +
                    titer_dose_gi * dose_gi +
                    titer_sex * sex.ii +
                    titer_age +
                    titer_species +
                    titer_assay +
                    titer_lab +
                    titer_time * median_peak # median peak
                  
                  
                  titer_sample <- rnorm(1, titer_mean, titer_sd) 
                  
                  
                  # Last positive time -----------------------------------------
                  
                  ## Get shape parameter ------
                  
                  shape_last_column <- paste0("shape_intercept_last[", organ.ii, "]")
                  #shape_last_column <- paste0("shape_intercept_last[", assay.ii, ",", organ.ii, "]")
                  shape_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_last_column)])
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_last[", organ.ii, "]")
                  median_intercept_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                  
                  median_location_column <-paste0("median_location_last[", organ.ii, "]") 
                  median_location_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_last[", route.ii, ",", organ.ii, "]") 
                  median_route_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
                  
                  median_dose_nose_column <- paste0("median_dose_last[1,", tissue.ii, "]") 
                  median_dose_nose_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_last[2,", tissue.ii, "]") 
                  median_dose_throat_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_last[3,", tissue.ii, "]") 
                  median_dose_trachea_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_last[4,", tissue.ii, "]") 
                  median_dose_lung_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_last[5,", tissue.ii, "]") 
                  median_dose_gi_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_last[", organ.ii, "]") 
                  median_sex_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_last[", age.ii, ",", organ.ii, "]")
                  median_age_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_last[", sp.ii, ",", organ.ii,  "]")
                  median_species_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_last[", assay.ii, ",", organ.ii, "]")
                  median_assay_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_last[", lab_standard, ",", organ.ii,  "]")
                  median_lab_last <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  median_time_column <- paste0("median_time_last[", inoc + 1, ",", organ.ii, "]")
                  median_time_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_last <- 0
                  }
                  
                  median_last <- exp(-1/10 * (median_intercept_last + 
                                                median_location_last * loc.ii +  
                                                median_route_last +  
                                                median_dose_nose_last * dose_nose +  
                                                median_dose_throat_last * dose_throat + 
                                                median_dose_trachea_last * dose_trachea + 
                                                median_dose_lung_last * dose_lung + 
                                                median_dose_gi_last * dose_gi + 
                                                median_species_last + 
                                                median_age_last +  
                                                median_sex_last * sex.ii +  
                                                median_assay_last + 
                                                median_lab_last +
                                                median_time_last * median_peak
                  ))
                  
                  scale_last <- median_last / (log(2)^(1 / shape_last))
                  
                  ## Make calculation ------
                  
                  #last_quantile <- scale_last * ((-log(1 - prob_level)) ^ (1/shape_last))
                  
                  
                  ## Draw sample based on peak time ----------------------------
                  
                  median_last_samp <- exp(-1/10 * (median_intercept_last + 
                                                     median_location_last * loc.ii +  
                                                     median_route_last +  
                                                     median_dose_nose_last * dose_nose +  
                                                     median_dose_throat_last * dose_throat + 
                                                     median_dose_trachea_last * dose_trachea + 
                                                     median_dose_lung_last * dose_lung + 
                                                     median_dose_gi_last * dose_gi + 
                                                     median_species_last + 
                                                     median_age_last +  
                                                     median_sex_last * sex.ii +  
                                                     median_assay_last + 
                                                     median_lab_last +
                                                     median_time_last * peak_sample
                  ))
                  
                  scale_last_samp <- median_last_samp / (log(2)^(1 / shape_last))
                  
                  
                  last_sample <- rweibull(1, 
                                          shape = shape_last, 
                                          scale = scale_last_samp)
                  
                  #print(percent)
                  #print(median_first)
                  #print(median_peak)
                  #print(titer_mean)
                  #print(median_last)
                  
                  
                  # Add estimates -------------------------------------------
                  
                  new.obs <- data.frame(sample_num = ii,
                                        organ_group = organ_group,
                                        organ_idx = organ.ii,
                                        dose_nose = dose_nose,
                                        dose_throat = dose_throat,
                                        dose_trachea = dose_trachea,
                                        dose_lung = dose_lung,
                                        dose_gi = dose_gi,
                                        dose_total = total_dose,
                                        route_idx = route.ii,
                                        sp_idx = sp.ii,
                                        age_idx = age.ii,
                                        sex_idx = sex.ii,
                                        assay_idx = assay.ii,
                                        location_idx = loc.ii,
                                        lab_idx = lab_standard,
                                        percent_positive = percent,
                                        first_pos_scale = scale_first,
                                        first_pos_shape = shape_first,
                                        first_pos_median = median_first,
                                        first_pos_sample = first_pos_sample,
                                        peak_scale = scale_peak,
                                        peak_shape = shape_peak,
                                        peak_median = median_peak,
                                        peak_sample = peak_sample,
                                        titer_mean = titer_mean,
                                        titer_sample = titer_sample,
                                        last_scale = scale_last,
                                        last_shape = shape_last,
                                        last_median = median_last,
                                        last_sample = last_sample
                  )
                  
                  p.fits <- rbind(p.fits, new.obs)
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  p.fits <- p.fits[-1, ]
  return(p.fits)
}


get_dose_estimates_joint <- function(fit, 
                                     draws = "all",
                                     n_draws = NA,
                                     dose_options = seq(1, 7, 0.25),
                                     sp_options = 1,
                                     sp_standard = 1,
                                     age_options = 2,
                                     age_standard = 2,
                                     sex_options = 0,
                                     sex_standard = 0,
                                     assay_options = 1, 
                                     lab_standard = 1,
                                     lab_effect = "No") {
  
  # Extract all parameters
  pars <- c("percent_intercept", "percent_location", #"percent_inoc",
            "percent_route", "percent_dose",
            "percent_sex", "percent_age", "percent_species",
            "percent_assay", "percent_lab",
            
            "shape_intercept_first",
            
            "median_intercept_first", "median_location_first",
            "median_route_first", "median_dose_first", 
            "median_sex_first", "median_age_first", "median_species_first",
            "median_assay_first", "median_lab_first", 
            
            "shape_intercept_peak",
            
            "median_intercept_peak", "median_location_peak",
            "median_route_peak", "median_dose_peak", 
            "median_sex_peak", "median_age_peak", "median_species_peak",
            "median_assay_peak", "median_lab_peak", 
            "median_time_peak",
            
            "titer_intercept", "titer_location", "titer_inoc",
            "titer_route", "titer_dose",
            "titer_sex", "titer_age", "titer_species",
            "titer_assay", "titer_lab", 
            "titer_time",
            "titer_intercept_sd", 
            
            "shape_intercept_last",
            
            "median_intercept_last", "median_location_last",
            "median_route_last", "median_dose_last", 
            "median_sex_last", "median_age_last", "median_species_last",
            "median_assay_last", "median_lab_last", 
            "median_time_last"
  )
  
  
  # Need wide format samples to get correlated values
  sample.pars <- fit$draws(pars, format = "df")
  
  p.fits <- data.frame(sample_num = NA,
                       organ_group = NA,
                       organ_idx = NA,
                       dose_nose = NA,
                       dose_throat = NA,
                       dose_trachea = NA,
                       dose_lung = NA,
                       dose_gi = NA,
                       dose_total = NA,
                       route_idx = NA,
                       sp_idx = NA,
                       age_idx = NA,
                       sex_idx = NA,
                       assay_idx = NA,
                       location_idx = NA,
                       lab_idx = NA,
                       percent_positive = NA,
                       first_pos_median = NA,
                       peak_median = NA,
                       titer_mean = NA,
                       last_median = NA
  )
  
  # Line num to loop through vectors of randomly sampled cofactors
  line_num <- 1
  
  # Get row numbers to loop through
  if (draws == "random") {
    sample_num <- sample(1:nrow(sample.pars), n_draws, replace = TRUE)
  }
  else if (draws == "all") {
    sample_num <- 1:nrow(sample.pars)
    cat("Using all parameter samples to generate estimates -- this may take some time!\n\n")
  }
  
  #print(sample_num)
  
  for (total.dose.ii in dose_options) {
    
    total_dose <- total.dose.ii
    assign_dose_distribution(total_dose)
    
    cat("\nGenerating estimates for total dose", total.dose.ii, "\n")
    
    for (route.ii in 1:5) {
      
      cat("Generating estimates for route", route.ii, "\n")
      
      dose_nose <- get(paste0("route", route.ii, ".dose.nose"))
      dose_throat <- get(paste0("route", route.ii, ".dose.throat"))
      dose_trachea <- get(paste0("route", route.ii, ".dose.trachea"))
      dose_lung <- get(paste0("route", route.ii, ".dose.lung"))
      dose_gi <- get(paste0("route", route.ii, ".dose.gi"))
      
      for (tissue.ii in 1:6) {
        
        # Set organ groups and location_idx
        if (tissue.ii == 1) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 0
        }
        else if (tissue.ii == 2) {
          organ_group <- "URT"
          organ.ii <- 1
          loc.ii <- 1
        }
        else if (tissue.ii == 3) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 0
        }
        else if (tissue.ii == 4) {
          organ_group <- "LRT"
          organ.ii <- 2
          loc.ii <- 1
        }
        else if (tissue.ii == 5) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 0
        }
        else if (tissue.ii == 6) {
          organ_group <- "GI"
          organ.ii <- 3
          loc.ii <- 1
        }
        
        # Set inoculated category
        inoc <- 0
        if (dose_nose > 0 & organ_group == "URT" & loc.ii == 0){inoc <- 1}
        if (dose_throat > 0 & organ_group == "URT" & loc.ii == 1){inoc <- 1}
        if (dose_trachea > 0 & organ_group == "LRT" & loc.ii == 0){inoc <- 1}
        if (dose_lung > 0 & organ_group == "LRT" & loc.ii == 1){inoc <- 1}
        if (dose_gi > 0 & organ_group == "GI" & loc.ii == 0){inoc <- 1}
        
        for (assay.ii in assay_options) {
          
          for (sex.ii in sex_options) {
            
            for (age.ii in age_options) {
              
              for (sp.ii in sp_options) {
                
                # Get the median estimate -----------------
                
                # First positive  -------------------------
                
                ### Cure fraction -------
                
                percent_intercept_column <-paste0("percent_intercept[", organ.ii, "]")
                percent_intercept <- as.numeric(fit$summary(percent_intercept_column)[3])
                
                percent_location_column <-paste0("percent_location[", organ.ii, "]") 
                percent_location <- as.numeric(fit$summary(percent_location_column)[3])
               
                percent_route_column <- paste0("percent_route[", route.ii, ",", organ.ii, "]") 
                percent_route <- as.numeric(fit$summary(percent_route_column)[3])
                
                percent_dose_nose_column <- paste0("percent_dose[1,", tissue.ii, "]") 
                percent_dose_nose <- as.numeric(fit$summary(percent_dose_nose_column)[3])
                
                percent_dose_throat_column <- paste0("percent_dose[2,", tissue.ii, "]") 
                percent_dose_throat <- as.numeric(fit$summary(percent_dose_throat_column)[3])
                
                percent_dose_trachea_column <- paste0("percent_dose[3,", tissue.ii, "]") 
                percent_dose_trachea <- as.numeric(fit$summary(percent_dose_trachea_column)[3])
                
                percent_dose_lung_column <- paste0("percent_dose[4,", tissue.ii, "]") 
                percent_dose_lung <- as.numeric(fit$summary(percent_dose_lung_column)[3])
                
                percent_dose_gi_column <- paste0("percent_dose[5,", tissue.ii, "]") 
                percent_dose_gi <- as.numeric(fit$summary(percent_dose_gi_column)[3])
                
                percent_sex_column <- paste0("percent_sex[", organ.ii, "]") 
                percent_sex <- as.numeric(fit$summary(percent_sex_column)[3])
                
                percent_age_column <- paste0("percent_age[", age.ii, ",", organ.ii, "]")
                percent_age <- as.numeric(fit$summary(percent_age_column)[3])
                
                percent_species_column <- paste0("percent_species[", sp.ii, ",", organ.ii, "]")
                percent_species <- as.numeric(fit$summary(percent_species_column)[3])
                
                percent_assay_column <- paste0("percent_assay[", assay.ii, ",", organ.ii, "]")
                percent_assay <- as.numeric(fit$summary(percent_assay_column)[3])
                
                percent_lab_column <- paste0("percent_lab[", lab_standard, ",", organ.ii,  "]")
                percent_lab <- as.numeric(fit$summary(percent_lab_column)[3])
                
                if (lab_effect == "No") {
                  percent_lab <- 0
                }
                
                percent_trans <- percent_intercept + 
                  percent_location * loc.ii +
                  percent_dose_nose * dose_nose +
                  percent_dose_throat * dose_throat +
                  percent_dose_trachea * dose_trachea +
                  percent_dose_lung * dose_lung +
                  percent_dose_gi * dose_gi +
                  percent_route + 
                  percent_age +
                  percent_species +
                  percent_sex * sex.ii +
                  percent_assay + 
                  percent_lab
                
                percent <- exp(percent_trans) / (1 + exp(percent_trans))
                
                #print(percent)
                
                ## Get shape parameter ------
                
                shape_first_column <- paste0("shape_intercept_first[", organ.ii, "]")
                shape_first <- as.numeric(fit$summary(shape_first_column)[3])
                
                
                ## Get scale parameter ----
                
                median_intercept_column <- paste0("median_intercept_first[", organ.ii, "]")
                median_intercept_first <- as.numeric(fit$summary(median_intercept_column)[3])
                
                median_location_column <-paste0("median_location_first[", organ.ii, "]") 
                median_location_first <- as.numeric(fit$summary(median_location_column)[3])
                
                median_route_column <- paste0("median_route_first[", route.ii, ",", organ.ii, "]") 
                median_route_first <- as.numeric(fit$summary(median_route_column)[3])
                
                median_dose_nose_column <- paste0("median_dose_first[1,", tissue.ii, "]") 
                median_dose_nose_first <- as.numeric(fit$summary(median_dose_nose_column)[3])
                
                median_dose_throat_column <- paste0("median_dose_first[2,", tissue.ii, "]") 
                median_dose_throat_first <- as.numeric(fit$summary(median_dose_throat_column)[3])
                
                median_dose_trachea_column <- paste0("median_dose_first[3,", tissue.ii, "]") 
                median_dose_trachea_first <- as.numeric(fit$summary(median_dose_trachea_column)[3])
                
                median_dose_lung_column <- paste0("median_dose_first[4,", tissue.ii, "]") 
                median_dose_lung_first <- as.numeric(fit$summary(median_dose_lung_column)[3])
                
                median_dose_gi_column <- paste0("median_dose_first[5,", tissue.ii, "]") 
                median_dose_gi_first <- as.numeric(fit$summary(median_dose_gi_column)[3])
                
                median_sex_column <- paste0("median_sex_first[", organ.ii, "]") 
                median_sex_first <- as.numeric(fit$summary(median_sex_column)[3])
                
                median_age_column <- paste0("median_age_first[", age.ii, ",", organ.ii, "]")
                median_age_first <- as.numeric(fit$summary(median_age_column)[3])
                
                median_species_column <- paste0("median_species_first[", sp.ii, ",", organ.ii,  "]")
                median_species_first <- as.numeric(fit$summary(median_species_column)[3])
                
                median_assay_column <- paste0("median_assay_first[", assay.ii, ",", organ.ii, "]")
                median_assay_first <- as.numeric(fit$summary(median_assay_column)[3])
                
                median_lab_column <- paste0("median_lab_first[", lab_standard, ",", organ.ii,  "]")
                median_lab_first <- as.numeric(fit$summary(median_lab_column)[3])
                
                if (lab_effect == "No") {
                  median_lab_first <- 0
                }
                
                median_first <- exp(-1/10 * (median_intercept_first + 
                                               median_location_first * loc.ii +  
                                               median_route_first +  
                                               median_dose_nose_first * dose_nose +  
                                               median_dose_throat_first * dose_throat + 
                                               median_dose_trachea_first * dose_trachea + 
                                               median_dose_lung_first * dose_lung + 
                                               median_dose_gi_first * dose_gi + 
                                               median_species_first + 
                                               median_age_first +  
                                               median_sex_first * sex.ii +  
                                               median_assay_first + 
                                               median_lab_first
                ))
                
                #print(median_first)
                
                # Peak time --------------------------------------------------
                
                ## Get shape parameter ------
                
                shape_peak_column <- paste0("shape_intercept_peak[", organ.ii, "]")
                shape_peak <- as.numeric(fit$summary(shape_peak_column)[3])
                
                ## Get scale parameter ----
                
                median_intercept_column <- paste0("median_intercept_peak[", organ.ii, "]")
                median_intercept_peak <- as.numeric(fit$summary(median_intercept_column)[3])
                
                median_location_column <-paste0("median_location_peak[", organ.ii, "]") 
                median_location_peak <- as.numeric(fit$summary(median_location_column)[3])
                
                median_route_column <- paste0("median_route_peak[", route.ii, ",", organ.ii, "]") 
                median_route_peak <- as.numeric(fit$summary(median_route_column)[3])
                
                median_dose_nose_column <- paste0("median_dose_peak[1,", tissue.ii, "]") 
                median_dose_nose_peak <- as.numeric(fit$summary(median_dose_nose_column)[3])
                
                median_dose_throat_column <- paste0("median_dose_peak[2,", tissue.ii, "]") 
                median_dose_throat_peak <- as.numeric(fit$summary(median_dose_throat_column)[3])
                
                median_dose_trachea_column <- paste0("median_dose_peak[3,", tissue.ii, "]") 
                median_dose_trachea_peak <- as.numeric(fit$summary(median_dose_trachea_column)[3])
                
                median_dose_lung_column <- paste0("median_dose_peak[4,", tissue.ii, "]") 
                median_dose_lung_peak <- as.numeric(fit$summary(median_dose_lung_column)[3])
                
                median_dose_gi_column <- paste0("median_dose_peak[5,", tissue.ii, "]") 
                median_dose_gi_peak <- as.numeric(fit$summary(median_dose_gi_column)[3])
                
                median_sex_column <- paste0("median_sex_peak[", organ.ii, "]") 
                median_sex_peak <- as.numeric(fit$summary(median_sex_column)[3])
                
                median_age_column <- paste0("median_age_peak[", age.ii, ",", organ.ii, "]")
                median_age_peak <- as.numeric(fit$summary(median_age_column)[3])
                
                median_species_column <- paste0("median_species_peak[", sp.ii, ",", organ.ii,  "]")
                median_species_peak <- as.numeric(fit$summary(median_species_column)[3])
                
                median_assay_column <- paste0("median_assay_peak[", assay.ii, ",", organ.ii, "]")
                median_assay_peak <- as.numeric(fit$summary(median_assay_column)[3])
                
                median_lab_column <- paste0("median_lab_peak[", lab_standard, ",", organ.ii,  "]")
                median_lab_peak <- as.numeric(fit$summary(median_lab_column)[3])
                
                median_time_column <- paste0("median_time_peak[", inoc + 1, ",", organ.ii, "]")
                median_time_peak <- as.numeric(fit$summary(median_time_column)[3])
                
                if (lab_effect == "No") {
                  median_lab_peak <- 0
                }
                
                median_peak <- exp(-1/10 * (median_intercept_peak + 
                                              median_location_peak * loc.ii +  
                                              median_route_peak +  
                                              median_dose_nose_peak * dose_nose +  
                                              median_dose_throat_peak * dose_throat + 
                                              median_dose_trachea_peak * dose_trachea + 
                                              median_dose_lung_peak * dose_lung + 
                                              median_dose_gi_peak * dose_gi + 
                                              median_species_peak + 
                                              median_age_peak +  
                                              median_sex_peak * sex.ii +  
                                              median_assay_peak + 
                                              median_lab_peak +
                                              median_time_peak * median_first
                ))
                #print(median_peak)
                
                # Peak titer -------------------------------------------------
                
                ## Get parameters --------------------------------------------
                
                titer_intercept_column <-paste0("titer_intercept[", organ.ii, "]") 
                titer_intercept <- as.numeric(fit$summary(titer_intercept_column)[3])
                
                titer_location_column <-paste0("titer_location[", organ.ii, "]") 
                titer_location <- as.numeric(fit$summary(titer_location_column)[3])
                
                titer_inoc_column <-paste0("titer_inoc[", organ.ii, "]") 
                titer_inoc <- as.numeric(fit$summary(titer_inoc_column)[3])
                
                titer_route_column <- paste0("titer_route[", route.ii, ",", organ.ii, "]") 
                titer_route <- as.numeric(fit$summary(titer_route_column)[3])
                
                titer_dose_nose_column <- paste0("titer_dose[1,", organ.ii, "]") 
                titer_dose_nose <- as.numeric(fit$summary(titer_dose_nose_column)[3])
                
                titer_dose_throat_column <- paste0("titer_dose[2,", organ.ii, "]") 
                titer_dose_throat <- as.numeric(fit$summary(titer_dose_throat_column)[3])
                
                titer_dose_trachea_column <- paste0("titer_dose[3,", organ.ii, "]") 
                titer_dose_trachea <- as.numeric(fit$summary(titer_dose_trachea_column)[3])
                
                titer_dose_lung_column <- paste0("titer_dose[4,", organ.ii, "]") 
                titer_dose_lung <- as.numeric(fit$summary(titer_dose_lung_column)[3])
                
                titer_dose_gi_column <- paste0("titer_dose[5,", organ.ii, "]") 
                titer_dose_gi <- as.numeric(fit$summary(titer_dose_gi_column)[3])
                
                titer_sex_column <- paste0("titer_sex[", organ.ii, "]")
                titer_sex <- as.numeric(fit$summary(titer_sex_column)[3])
                
                titer_age_column <- paste0("titer_age[", age.ii, ",", organ.ii, "]")
                titer_age <- as.numeric(fit$summary(titer_age_column)[3])
                
                titer_species_column <- paste0("titer_species[", sp.ii, ",", organ.ii, "]")
                titer_species <- as.numeric(fit$summary(titer_species_column)[3])
                
                titer_assay_column <- paste0("titer_assay[", assay.ii, ",", organ.ii, "]")
                titer_assay <- as.numeric(fit$summary(titer_assay_column)[3])
                
                titer_lab_column <- paste0("titer_lab[", lab_standard, ",", organ.ii,  "]")
                titer_lab <- as.numeric(fit$summary(titer_lab_column)[3])
                
                titer_time_column <- paste0("titer_time[", inoc + 1, ",", organ.ii, "]")
                titer_time <- as.numeric(fit$summary(titer_time_column)[3])
                
                titer_sd_column <- paste0("titer_intercept_sd[", organ.ii, "]")
                titer_sd <- as.numeric(fit$summary(titer_sd_column)[3])
                
                if (lab_effect == "No") {
                  titer_lab <- 0
                }
                
                
                ## Make calculation ------------------------------------------
                
                titer_mean <- titer_intercept +
                  titer_location * loc.ii + 
                  titer_inoc * inoc + 
                  titer_route +
                  titer_dose_nose * dose_nose +
                  titer_dose_throat * dose_throat +
                  titer_dose_trachea * dose_trachea +
                  titer_dose_lung * dose_lung +
                  titer_dose_gi * dose_gi +
                  titer_sex * sex.ii +
                  titer_age +
                  titer_species +
                  titer_assay +
                  titer_lab +
                  titer_time * median_peak # median peak
                
                #print(titer_mean)
                
                # Last positive time -----------------------------------------
                
                ## Get shape parameter ------
                
                shape_last_column <- paste0("shape_intercept_last[", organ.ii, "]")
                shape_last <- as.numeric(fit$summary(shape_last_column)[3])
                
                ## Get scale parameter ----
                
                median_intercept_column <- paste0("median_intercept_last[", organ.ii, "]")
                median_intercept_last <- as.numeric(fit$summary(median_intercept_column)[3])
                
                median_location_column <-paste0("median_location_last[", organ.ii, "]") 
                median_location_last <- as.numeric(fit$summary(median_location_column)[3])
                
                median_route_column <- paste0("median_route_last[", route.ii, ",", organ.ii, "]") 
                median_route_last <- as.numeric(fit$summary(median_route_column)[3])
                
                median_dose_nose_column <- paste0("median_dose_last[1,", tissue.ii, "]") 
                median_dose_nose_last <- as.numeric(fit$summary(median_dose_nose_column)[3])
                
                median_dose_throat_column <- paste0("median_dose_last[2,", tissue.ii, "]") 
                median_dose_throat_last <- as.numeric(fit$summary(median_dose_throat_column)[3])
                
                median_dose_trachea_column <- paste0("median_dose_last[3,", tissue.ii, "]") 
                median_dose_trachea_last <- as.numeric(fit$summary(median_dose_trachea_column)[3])
                
                median_dose_lung_column <- paste0("median_dose_last[4,", tissue.ii, "]") 
                median_dose_lung_last <- as.numeric(fit$summary(median_dose_lung_column)[3])
                
                median_dose_gi_column <- paste0("median_dose_last[5,", tissue.ii, "]") 
                median_dose_gi_last <- as.numeric(fit$summary(median_dose_gi_column)[3])
                
                median_sex_column <- paste0("median_sex_last[", organ.ii, "]") 
                median_sex_last <- as.numeric(fit$summary(median_sex_column)[3])
                
                median_age_column <- paste0("median_age_last[", age.ii, ",", organ.ii, "]")
                median_age_last <- as.numeric(fit$summary(median_age_column)[3])
                
                median_species_column <- paste0("median_species_last[", sp.ii, ",", organ.ii,  "]")
                median_species_last <- as.numeric(fit$summary(median_species_column)[3])
                
                median_assay_column <- paste0("median_assay_last[", assay.ii, ",", organ.ii, "]")
                median_assay_last <- as.numeric(fit$summary(median_assay_column)[3])
                
                median_lab_column <- paste0("median_lab_last[", lab_standard, ",", organ.ii,  "]")
                median_lab_last <- as.numeric(fit$summary(median_lab_column)[3])
                
                median_time_column <- paste0("median_time_last[", inoc + 1, ",", organ.ii, "]")
                median_time_last <- as.numeric(fit$summary(median_time_column)[3])
                
                if (lab_effect == "No") {
                  median_lab_last <- 0
                }
                
                median_last <- exp(-1/10 * (median_intercept_last + 
                                              median_location_last * loc.ii +  
                                              median_route_last +  
                                              median_dose_nose_last * dose_nose +  
                                              median_dose_throat_last * dose_throat + 
                                              median_dose_trachea_last * dose_trachea + 
                                              median_dose_lung_last * dose_lung + 
                                              median_dose_gi_last * dose_gi + 
                                              median_species_last + 
                                              median_age_last +  
                                              median_sex_last * sex.ii +  
                                              median_assay_last + 
                                              median_lab_last +
                                              median_time_last * median_peak
                ))
                
                #print(median_last)
                
                # Add estimates -------------------------------------------
                
                new.obs <- data.frame(sample_num = "Median",
                                      organ_group = organ_group,
                                      organ_idx = organ.ii,
                                      dose_nose = dose_nose,
                                      dose_throat = dose_throat,
                                      dose_trachea = dose_trachea,
                                      dose_lung = dose_lung,
                                      dose_gi = dose_gi,
                                      dose_total = total_dose,
                                      route_idx = route.ii,
                                      sp_idx = sp.ii,
                                      age_idx = age.ii,
                                      sex_idx = sex.ii,
                                      assay_idx = assay.ii,
                                      location_idx = loc.ii,
                                      lab_idx = lab_standard,
                                      percent_positive = percent,
                                      first_pos_median = median_first,
                                      peak_median = median_peak,
                                      titer_mean = titer_mean,
                                      last_median = median_last)
                #print(new.obs)
                
                p.fits <- rbind(p.fits, new.obs)
                
                #print(sample_num)
                for (ii in sample_num) {
                  
                  # First positive  -------------------------
                  
                  ### Cure fraction -------
                  
                  percent_intercept_column <-paste0("percent_intercept[", organ.ii, "]") 
                  percent_intercept <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_intercept_column)])
                  
                  percent_location_column <-paste0("percent_location[", organ.ii, "]") 
                  percent_location <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_location_column)])
                  
                  percent_route_column <- paste0("percent_route[", route.ii, ",", organ.ii, "]") 
                  percent_route <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_route_column)])
                  
                  percent_dose_nose_column <- paste0("percent_dose[1,", tissue.ii, "]") 
                  percent_dose_nose <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_nose_column)])
                  
                  percent_dose_throat_column <- paste0("percent_dose[2,", tissue.ii, "]") 
                  percent_dose_throat <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_throat_column)])
                  
                  percent_dose_trachea_column <- paste0("percent_dose[3,", tissue.ii, "]") 
                  percent_dose_trachea <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_trachea_column)])
                  
                  percent_dose_lung_column <- paste0("percent_dose[4,", tissue.ii, "]") 
                  percent_dose_lung <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_lung_column)])
                  
                  percent_dose_gi_column <- paste0("percent_dose[5,", tissue.ii, "]") 
                  percent_dose_gi <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_dose_gi_column)])
                  
                  percent_sex_column <- paste0("percent_sex[", organ.ii, "]") 
                  percent_sex <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_sex_column)])
                  
                  percent_age_column <- paste0("percent_age[", age.ii, ",", organ.ii, "]")
                  percent_age <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_age_column)])
                  
                  percent_species_column <- paste0("percent_species[", sp.ii, ",", organ.ii, "]")
                  percent_species <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_species_column)])
                  
                  percent_assay_column <- paste0("percent_assay[", assay.ii, ",", organ.ii, "]")
                  percent_assay <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_assay_column)])
                  
                  percent_lab_column <- paste0("percent_lab[", lab_standard, ",", organ.ii,  "]")
                  percent_lab <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == percent_lab_column)])
                  
                  if (lab_effect == "No") {
                    percent_lab <- 0
                  }
                  
                  percent_trans <- percent_intercept + 
                    #percent_inoc * inoc + 
                    percent_location * loc.ii +
                    percent_dose_nose * dose_nose +
                    percent_dose_throat * dose_throat +
                    percent_dose_trachea * dose_trachea +
                    percent_dose_lung * dose_lung +
                    percent_dose_gi * dose_gi +
                    percent_route + 
                    percent_age +
                    percent_species +
                    percent_sex * sex.ii +
                    percent_assay + 
                    percent_lab
                  
                  percent <- exp(percent_trans) / (1 + exp(percent_trans))
                  
                  
                  ## Get shape parameter ------
                  
                  shape_first_column <- paste0("shape_intercept_first[", organ.ii, "]")
                  shape_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_first_column)])
                  
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_first[", organ.ii, "]")
                  median_intercept_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                  
                  median_location_column <-paste0("median_location_first[", organ.ii, "]") 
                  median_location_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_first[", route.ii, ",", organ.ii, "]") 
                  median_route_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
                  
                  median_dose_nose_column <- paste0("median_dose_first[1,", tissue.ii, "]") 
                  median_dose_nose_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_first[2,", tissue.ii, "]") 
                  median_dose_throat_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_first[3,", tissue.ii, "]") 
                  median_dose_trachea_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_first[4,", tissue.ii, "]") 
                  median_dose_lung_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_first[5,", tissue.ii, "]") 
                  median_dose_gi_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_first[", organ.ii, "]") 
                  median_sex_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_first[", age.ii, ",", organ.ii, "]")
                  median_age_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_first[", sp.ii, ",", organ.ii,  "]")
                  median_species_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_first[", assay.ii, ",", organ.ii, "]")
                  median_assay_first <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_first[", lab_standard, ",", organ.ii,  "]")
                  median_lab_first <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_first <- 0
                  }
                  
                  median_first <- exp(-1/10 * (median_intercept_first + 
                                                 median_location_first * loc.ii +  
                                                 median_route_first +  
                                                 median_dose_nose_first * dose_nose +  
                                                 median_dose_throat_first * dose_throat + 
                                                 median_dose_trachea_first * dose_trachea + 
                                                 median_dose_lung_first * dose_lung + 
                                                 median_dose_gi_first * dose_gi + 
                                                 median_species_first + 
                                                 median_age_first +  
                                                 median_sex_first * sex.ii +  
                                                 median_assay_first + 
                                                 median_lab_first
                  ))
                  
                  scale_first <- median_first / (log(2)^(1 / shape_first))
                
                  
                  # Peak time --------------------------------------------------
                  
                  ## Get shape parameter ------
                  
                  shape_peak_column <- paste0("shape_intercept_peak[", organ.ii, "]")
                  #shape_peak_column <- paste0("shape_intercept_peak[", assay.ii, ",", organ.ii, "]")
                  shape_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_peak_column)])
                  
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_peak[", organ.ii, "]")
                  median_intercept_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                  
                  median_location_column <-paste0("median_location_peak[", organ.ii, "]") 
                  median_location_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_peak[", route.ii, ",", organ.ii, "]") 
                  median_route_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
                  
                  median_dose_nose_column <- paste0("median_dose_peak[1,", tissue.ii, "]") 
                  median_dose_nose_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_peak[2,", tissue.ii, "]") 
                  median_dose_throat_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_peak[3,", tissue.ii, "]") 
                  median_dose_trachea_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_peak[4,", tissue.ii, "]") 
                  median_dose_lung_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_peak[5,", tissue.ii, "]") 
                  median_dose_gi_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_peak[", organ.ii, "]") 
                  median_sex_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_peak[", age.ii, ",", organ.ii, "]")
                  median_age_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_peak[", sp.ii, ",", organ.ii,  "]")
                  median_species_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_peak[", assay.ii, ",", organ.ii, "]")
                  median_assay_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_peak[", lab_standard, ",", organ.ii,  "]")
                  median_lab_peak <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  median_time_column <- paste0("median_time_peak[", inoc + 1, ",", organ.ii, "]")
                  median_time_peak <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_peak <- 0
                  }
                  
                  median_peak <- exp(-1/10 * (median_intercept_peak + 
                                                median_location_peak * loc.ii +  
                                                median_route_peak +  
                                                median_dose_nose_peak * dose_nose +  
                                                median_dose_throat_peak * dose_throat + 
                                                median_dose_trachea_peak * dose_trachea + 
                                                median_dose_lung_peak * dose_lung + 
                                                median_dose_gi_peak * dose_gi + 
                                                median_species_peak + 
                                                median_age_peak +  
                                                median_sex_peak * sex.ii +  
                                                median_assay_peak + 
                                                median_lab_peak +
                                                median_time_peak * median_first
                  ))
                
                  
                  # Peak titer -------------------------------------------------
                  
                  ## Get parameters --------------------------------------------
                  
                  titer_intercept_column <-paste0("titer_intercept[", organ.ii, "]") 
                  titer_intercept <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_intercept_column)])
                  
                  titer_location_column <-paste0("titer_location[", organ.ii, "]") 
                  titer_location <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_location_column)])
                  
                  titer_inoc_column <-paste0("titer_inoc[", organ.ii, "]") 
                  titer_inoc <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_inoc_column)])
                  
                  titer_route_column <- paste0("titer_route[", route.ii, ",", organ.ii, "]") 
                  titer_route <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_route_column)])
                  
                  titer_dose_nose_column <- paste0("titer_dose[1,", organ.ii, "]") 
                  titer_dose_nose <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_nose_column)])
                  
                  titer_dose_throat_column <- paste0("titer_dose[2,", organ.ii, "]") 
                  titer_dose_throat <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_throat_column)])
                  
                  titer_dose_trachea_column <- paste0("titer_dose[3,", organ.ii, "]") 
                  titer_dose_trachea <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_trachea_column)])
                  
                  titer_dose_lung_column <- paste0("titer_dose[4,", organ.ii, "]") 
                  titer_dose_lung <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_lung_column)])
                  
                  titer_dose_gi_column <- paste0("titer_dose[5,", organ.ii, "]") 
                  titer_dose_gi <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_dose_gi_column)])
                  
                  titer_sex_column <- paste0("titer_sex[", organ.ii, "]")
                  titer_sex <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_sex_column)])
                  
                  titer_age_column <- paste0("titer_age[", age.ii, ",", organ.ii, "]")
                  titer_age <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_age_column)])
                  
                  titer_species_column <- paste0("titer_species[", sp.ii, ",", organ.ii, "]")
                  titer_species <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_species_column)])
                  
                  titer_assay_column <- paste0("titer_assay[", assay.ii, ",", organ.ii, "]")
                  titer_assay <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_assay_column)])
                  
                  titer_lab_column <- paste0("titer_lab[", lab_standard, ",", organ.ii,  "]")
                  titer_lab <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_lab_column)])
                  
                  titer_time_column <- paste0("titer_time[", inoc + 1, ",", organ.ii, "]")
                  titer_time <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == titer_time_column)])
                  
                  titer_sd <- sample.pars$titer_intercept_sd[ii]
                  
                  if (lab_effect == "No") {
                    titer_lab <- 0
                  }
                  
                  
                  ## Make calculation ------------------------------------------
                  
                  titer_mean <- titer_intercept +
                    titer_location * loc.ii + 
                    titer_inoc * inoc + 
                    titer_route +
                    titer_dose_nose * dose_nose +
                    titer_dose_throat * dose_throat +
                    titer_dose_trachea * dose_trachea +
                    titer_dose_lung * dose_lung +
                    titer_dose_gi * dose_gi +
                    titer_sex * sex.ii +
                    titer_age +
                    titer_species +
                    titer_assay +
                    titer_lab +
                    titer_time * median_peak # median peak
                  
                  
                  # Last positive time -----------------------------------------
                  
                  ## Get shape parameter ------
                  
                  shape_last_column <- paste0("shape_intercept_last[", organ.ii, "]")
                  shape_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == shape_last_column)])
                  
                  ## Get scale parameter ----
                  
                  median_intercept_column <- paste0("median_intercept_last[", organ.ii, "]")
                  median_intercept_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_intercept_column)])
                  
                  median_location_column <-paste0("median_location_last[", organ.ii, "]") 
                  median_location_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_location_column)])
                  
                  median_route_column <- paste0("median_route_last[", route.ii, ",", organ.ii, "]") 
                  median_route_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_route_column)])
                  
                  median_dose_nose_column <- paste0("median_dose_last[1,", tissue.ii, "]") 
                  median_dose_nose_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_nose_column)])
                  
                  median_dose_throat_column <- paste0("median_dose_last[2,", tissue.ii, "]") 
                  median_dose_throat_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_throat_column)])
                  
                  median_dose_trachea_column <- paste0("median_dose_last[3,", tissue.ii, "]") 
                  median_dose_trachea_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_trachea_column)])
                  
                  median_dose_lung_column <- paste0("median_dose_last[4,", tissue.ii, "]") 
                  median_dose_lung_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_lung_column)])
                  
                  median_dose_gi_column <- paste0("median_dose_last[5,", tissue.ii, "]") 
                  median_dose_gi_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_dose_gi_column)])
                  
                  median_sex_column <- paste0("median_sex_last[", organ.ii, "]") 
                  median_sex_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_sex_column)])
                  
                  median_age_column <- paste0("median_age_last[", age.ii, ",", organ.ii, "]")
                  median_age_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_age_column)])
                  
                  median_species_column <- paste0("median_species_last[", sp.ii, ",", organ.ii,  "]")
                  median_species_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_species_column)])
                  
                  median_assay_column <- paste0("median_assay_last[", assay.ii, ",", organ.ii, "]")
                  median_assay_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_assay_column)])
                  
                  median_lab_column <- paste0("median_lab_last[", lab_standard, ",", organ.ii,  "]")
                  median_lab_last <- as.numeric(sample.pars[ii, which(colnames(sample.pars) == median_lab_column)])
                  
                  median_time_column <- paste0("median_time_last[", inoc + 1, ",", organ.ii, "]")
                  median_time_last <- as.numeric(
                    sample.pars[ii, which(colnames(sample.pars) == median_time_column)])
                  
                  if (lab_effect == "No") {
                    median_lab_last <- 0
                  }
                  
                  median_last <- exp(-1/10 * (median_intercept_last + 
                                                median_location_last * loc.ii +  
                                                median_route_last +  
                                                median_dose_nose_last * dose_nose +  
                                                median_dose_throat_last * dose_throat + 
                                                median_dose_trachea_last * dose_trachea + 
                                                median_dose_lung_last * dose_lung + 
                                                median_dose_gi_last * dose_gi + 
                                                median_species_last + 
                                                median_age_last +  
                                                median_sex_last * sex.ii +  
                                                median_assay_last + 
                                                median_lab_last +
                                                median_time_last * median_peak
                  ))
                  
                  
                  
                  # Add estimates -------------------------------------------
                  
                  new.obs <- data.frame(sample_num = ii,
                                        organ_group = organ_group,
                                        organ_idx = organ.ii,
                                        dose_nose = dose_nose,
                                        dose_throat = dose_throat,
                                        dose_trachea = dose_trachea,
                                        dose_lung = dose_lung,
                                        dose_gi = dose_gi,
                                        dose_total = total_dose,
                                        route_idx = route.ii,
                                        sp_idx = sp.ii,
                                        age_idx = age.ii,
                                        sex_idx = sex.ii,
                                        assay_idx = assay.ii,
                                        location_idx = loc.ii,
                                        lab_idx = lab_standard,
                                        percent_positive = percent,
                                        first_pos_median = median_first,
                                        peak_median = median_peak,
                                        titer_mean = titer_mean,
                                        last_median = median_last
                  )
                  
                  p.fits <- rbind(p.fits, new.obs)
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  p.fits <- p.fits[-1, ]
  return(p.fits)
}



get_weibull_curves <- function(estimates, 
                               event = "first positive", 
                               x_vals = seq(0, 20, 0.1), 
                               include_cure = FALSE) {
  
  n_rows <- nrow(estimates)
  n_x <- length(x_vals)
  
  shape_col <- switch(event,
                      "first positive" = "first_shape",
                      "peak time" = "peak_shape",
                      "last positive" = "last_shape")
  
  scale_col <- switch(event,
                      "first positive" = "first_scale",
                      "peak time" = "peak_scale",
                      "last positive" = "last_scale")
  
  shapes <- estimates[[shape_col]]
  scales <- estimates[[scale_col]]
  
  # Repeat values once
  x_rep <- rep(x_vals, times = n_rows)
  shape_rep <- rep(shapes, each = n_x)
  scale_rep <- rep(scales, each = n_x)
  
  pdf <- dweibull(x_rep, shape = shape_rep, scale = scale_rep)
  cdf <- pweibull(x_rep, shape = shape_rep, scale = scale_rep)
  
  if (include_cure) {
    cure_rep <- rep(estimates$percent_positive, each = n_x)
    cdf <- cure_rep * cdf
  }
  
  df.surv <- estimates[rep(seq_len(n_rows), each = n_x), ]
  
  df.surv$day_post_infection <- x_rep
  df.surv$pdf <- pdf
  df.surv$cdf <- cdf
  df.surv$line_num <- rep(seq_len(n_rows), each = n_x)
  
  return(df.surv)
}

# Generates colors
assign_colors <- function() {
  
  df <- data.frame(Nose = "#85649B",
                   Throat = "#3A9AB2",
                   Trachea = "#BDC881",
                   Lung = "#E3B710",
                   'Upper GI' = "#EC7A05",
                   'Lower GI' = "#F11B00")
  
  return(df)
  
}


