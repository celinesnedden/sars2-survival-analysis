functions{

  vector lower_bound_adjusted(vector wait_times) {
    
    vector[num_elements(wait_times)] adjusted_wait_times = wait_times;
    
    for (ii in 1:num_elements(wait_times)) {
      if (wait_times[ii] < 0) {
         adjusted_wait_times[ii] = -1e-100; 
      }
    }
    
    return(adjusted_wait_times);
  }
  
  real lower_bound_adjusted_real(real wait_time) {
    
    real adjusted_wait_time = wait_time;
    
    if (wait_time < 0) {
      adjusted_wait_time = 0;
    }
    
    return(adjusted_wait_time);
  }
  
  vector upper_bound_adjusted(vector wait_times) {
    
    vector[num_elements(wait_times)] adjusted_wait_times = wait_times;
    
    for (ii in 1:num_elements(wait_times)) {
      if (wait_times[ii] <= 0) {
        adjusted_wait_times[ii] = 0.01;
      }
    }
    
    return(adjusted_wait_times);
  }
  
  // Calculate the number of unique elements, by checking repeats
  int num_unique(array[] int x) {
    
    array[num_elements(x)] int res = sort_desc(x);
    
    // Count how many times there are repeats
    int counter = 0;
    for (jj in 2:num_elements(res)) {
      if (res[jj] ==  res[jj-1]) {
        counter = counter + 1;
      };
    }
    
    // Return the difference of the input length and number of repeated elements
    return num_elements(x) - counter;
  }
  
  array[] int get_unique_vec(array[] int x) {
    array[num_elements(x)] int x_sorted = sort_desc(x);
    array[num_unique(x)] int unique_vec;
    int counter = 2; 
    
    // Loop to find & assign the unique elements
    unique_vec[1] = x[1];
    for (jj in 2:num_elements(x)) {
      if (x_sorted[jj] != x_sorted[jj-1]) {
        unique_vec[counter] = x_sorted[jj];
        counter = counter + 1;
      }
    }
    return unique_vec;
  }
  
  // Check whether marginalization is necessary for a given data point
  //  for unknown sex, age, assay type, or tissue location
  real check_for_marginalizing(real sex, real age, 
                               real assay, real location) {
    if (sex != -9999 && age != -9999 && 
        assay != -9999 && location != -9999) {
      return(0.0);
    }
    else {
      return(1.0);
    }
  }
  
  // Calculate transformed percent parameter
  real calculate_percent(real intercept, 
                         real location,
                         real route,
                         real dose_nose,
                         real dose_throat,
                         real dose_trachea,
                         real dose_lung,
                         real dose_gi,
                         real sex,
                         real age,
                         real species,
                         real assay,
                         real lab) {
                           
   real percent = inv_logit(intercept + location + route + dose_nose +
                            dose_throat + dose_trachea + dose_lung + dose_gi + 
                            sex + age +species + assay + lab
                            );
    
   return(percent);          
  }
  
  
  // Calculate transformed median parameter
  real calculate_median(real intercept, 
                        real location,
                        real route,
                        real dose_nose,
                        real dose_throat,
                        real dose_trachea,
                        real dose_lung,
                        real dose_gi,
                        real sex,
                        real age,
                        real species,
                        real assay,
                        real lab,
                        real time) {
                           
   real median = exp(-1/10.0 * (intercept + location + route + dose_nose +
                      dose_throat + dose_trachea + dose_lung + 
                      dose_gi + sex + age + species + assay + lab + time
                      ));
    
   return(median);          
  }
  
  // Calculate the Weibull scale parameter from the median and shape
  real calculate_scale(real median, real shape) {
      real scale = median / (log(2)^(1 / shape));
      
      if (scale == 0) {
        scale = 1e-100;
      }
      else if (is_inf(scale)) {
        scale = 1e100;
      }
      
      return (scale);
  }
  
  // Get the options for marginizalition of a given cofactor
  array[] int get_cofactor_options(int cofactor, // -9999 if unknown
                                   array[] int options_vector) {

      if (cofactor == -9999) {
        return(options_vector);
      } 
      else {
        return(rep_array(cofactor, num_elements(options_vector)));
      }
  }
  
  // Get the log probabilities of a cofactor when marginalizing
  real get_cofactor_probs(int cofactor, real cofactor_prob) {
    if (cofactor == -9999) {
      return(log(cofactor_prob));
    }
    else {
      return(log(1));
    }
  }
  
  // Calculate log likelihood
  real calculate_log_lik(real ever_pos, 
                         real event_time,
                         real true_wait,
                         real first_lb,
                         real percent, 
                         real shape, 
                         real median) {
                           
      // Calculate scale using shape & median
      real scale = calculate_scale(median, shape);
      real log_lik = negative_infinity(); // Will be changed below
      
      // If the individual never tested positive
      if (ever_pos == 0 && event_time == 1) {
        log_lik = log_sum_exp(log1m(percent),
                              log(percent) +
                              weibull_lccdf(first_lb | shape, scale));
      }
      // If they do test positive, and the observation is interval censored
      else if (ever_pos == 1 && event_time == 1) {
        log_lik = log(percent) + 
                  weibull_lpdf(true_wait | shape, scale);
      }
      // For peak time and time to undetectability, we do not include a cure fraction 
      //    (all individuals who have a peak time or last positive time 
      //     by definition have tested postiive)
      else if (ever_pos == 1 && event_time == 2) {
        log_lik = weibull_lpdf(true_wait | shape, scale);
      }
      
      return(log_lik);
  }
  
}

data {
  
   // Data
   int N;  // total number of samples
   int L_assay; // the number of different assay categories (e.g., total RNA vs. PFU)
   int L_lab;  // the number of different labs included with event times
   int L_lab_titer;  // the number of labs with quantitative peak titer information
   vector<lower=0, upper=1>[N] ever_positive; // whether an individual ever tests positive, where 0: never positive; 1: positive
   vector<lower=0, upper=1>[N] has_peak; // whether an individual has a censored peak time, where 0: doesn't contribute, 1: does contribute
   vector<lower=0, upper=1>[N] has_titer; // whether an individual has an observed peak titer, where 0: doesn't contribute, 1: does contribute
   vector<lower=0, upper=1>[N] has_last_positive; // whether an individual has a censored last positive time, where 0: does not contribute; 1: does contribute
   vector<lower=0, upper=1>[N] last_positive_type; // the type of censored observation for the last positive, where 0: right censored; 1: interval censored
   
   // Location information
   array[N] int<lower=1, upper=3> organ_location; // the organ group for each datapoint, where: 1: URT; 2: LRT; 3: GI
   array[N] int<upper=6> tissue_location; // the tissue location for each datapoint, where: 1: nose; 2: throat; 3: trachea; 4: lung; 5: upper GI; 6: lower GI 
   array[N] int<upper=1> location; // distinguishes between tissue locations within an organ group (e.g., between nose & throat for URT)
   
   // Cofactors
   array[N] int<lower=0, upper=1> inoc;  // whether the sampled tissue was exposed for each datapoint, where 0: not exposed; 1: exposed
   array[N] int<lower=1, upper=5> route; // the exposure route category for each data point
   array[N] real<lower=0> dose_nose;  // the dose administered to the nose for the individual the datapoint comes from (range: 0-1; log10 rescaled)
   array[N] real<lower=0> dose_throat; // the dose administered to the throat for the individual the datapoint comes from (range: 0-1; log10 rescaled)
   array[N] real<lower=0> dose_trachea; // the dose administered to the trachea for the individual the datapoint comes from (range: 0-1; log10 rescaled)
   array[N] real<lower=0> dose_lung; // the dose administered to the lung for the individual the datapoint comes from (range: 0-1; log10 rescaled)
   array[N] real<lower=0> dose_gi; // the dose administered to the upper gi for the individual the datapoint comes from (range: 0-1; log10 rescaled)
   array[N] int<upper=1> sex;  // the sex of the individual each datapoint comes from, where 0: female, 1: male; -9999: unknown
   array[N] int<upper=3> age;  // the age class of the individual each datapoint comes from, where 1: juvenile, 2: adult, 3: geriatric, -9999: unknown
   array[N] int<upper=3> species; // the species of the individual each datapoint comes from, where 1: rhesus macaque, 2: cynomolgus macaque, 3: african green monkey
   array[N] int<upper=6> assay; // the assay category for each datapoint; range 1-5 integers only
   array[N] int<lower=0, upper=L_lab> lab; // the lab for each datapoint; integers only
   
   // Male/female distributions 
   array[N] real<lower=0,upper=1> prob_male; // the probability a datapoint comes from a male, if unknown, given M/F distributions in the whole study
   array[N] real<lower=0,upper=1> prob_female; // the probability a datapoint comes from a female, if unknown, given M/F distributions in the whole study
   
   // Age distributions
   array[N] real<lower=0,upper=1> prob_geriatric; // the probability a datapoint comes from a geriatric, given age distributions in the whole study
   
   // Bounds on event times
   vector<lower=0>[N] first_lower_bounds; // lower bounds on time to detectability
   vector<lower=first_lower_bounds>[N] first_upper_bounds; // upper bounds on time to detectability
   
   vector<lower=first_lower_bounds>[N] peak_lower_bounds; // lower bounds on time to peak titer
   vector<lower=peak_lower_bounds>[N] peak_upper_bounds;  // upper bounds on time to peak titer
   
   vector<lower=peak_lower_bounds>[N] last_lower_bounds; // lower bounds on time to undetectability
   vector<lower=last_lower_bounds>[N] last_upper_bounds; // upper bounds on time to undetectability
   
   // Observed peak titer values
   vector[N] peak_observed_titer;
}

parameters {
  
   // Cure fraction (i.e., percent of individuals that ever test positive)
   array[3] real percent_intercept; // intercept for each organ group, 1: URT; 2: LRT; 3: GI
   array[3] real percent_location; // the effect of location for each group, 1: URT; 2: LRT; 3: GI
   matrix[5, 6] percent_route; // the effect of exposure route, where rows: routes; columns: tissues sampled 
   matrix<lower=0>[5, 6] percent_dose; // the effect of dose, where rows: dose in tissue locations; columns: tissues sampled
   array[3] real percent_sex; // the effect of sex for each organ group, 1: URT; 2: LRT; 3: GI
   matrix[3, 3] percent_age; // the effects of age classes for each organ group, where rows: age classes; columns: organ groups
   matrix[3, 3] percent_species; // the effects of species for each organ group, where rows: species; columns: organ groups
   matrix[L_assay, 3] percent_assay; // the effects of assay type for each organ group, where rows: assays; columns: organ groups
   matrix[L_lab, 3] percent_lab; // the effects of lab for each organ group, where rows: lab groups; columns: organ groups
   
   // Shape for the time to detectabiltiy 
   array[3] real<lower=1> shape_intercept_first; // the shape parameter for each organ group, 1: URT; 2: LRT; 3: GI
   
   // Median for the time to detectability
   array[3] real median_intercept_first; /// intercept for each organ group, 1: URT; 2: LRT; 3: GI
   array[3] real median_location_first;  // the effect of location for each group, 1: URT; 2: LRT; 3: GI
   matrix[5, 6] median_route_first;      // the effect of exposure route, where rows: routes; columns: tissues sampled 
   matrix[5, 6] median_dose_first;       // the effect of dose, where rows: dose in tissue locations; columns: tissues sampled
   array[3] real median_sex_first;       // the effect of sex for each organ group, 1: URT; 2: LRT; 3: GI
   matrix[3, 3] median_age_first;        // the effects of age classes for each organ group, where rows: age classes; columns: organ groups
   matrix[3, 3] median_species_first;    // the effects of species for each organ group, where rows: species; columns: organ groups
   matrix[L_assay, 3] median_assay_first;// the effects of assay type for each organ group, where rows: assays; columns: organ groups
   matrix[L_lab, 3] median_lab_first;    // the effects of lab for each organ group, where rows: lab groups; columns: organ groups
   
   // Shape for the time to peak titer
   array[3] real<lower=1> shape_intercept_peak;
   
   // Median for the time to peak titer
   array[3] real median_intercept_peak; /// intercept for each organ group,
   array[3] real median_location_peak;  // the effect of location for each group, 
   matrix[5, 6] median_route_peak; // the effect of exposure route, where rows: routes; columns: tissues sampled 
   matrix[5, 6] median_dose_peak; // the effect of dose, where rows: dose in tissue locations; columns: tissues sampled
   array[3] real median_sex_peak; // the effect of sex for each organ group,
   matrix[3, 3] median_age_peak; // the effects of age classes for each organ group, where rows: age classes; columns: organ groups
   matrix[3, 3] median_species_peak; // the effects of species for each organ group, where rows: species; columns: organ groups
   matrix[L_assay, 3] median_assay_peak; // the effects of assay type for each organ group, where rows: assays; columns: organ groups
   matrix[L_lab, 3] median_lab_peak; // the effects of lab for each organ group, where rows: lab groups; columns: organ groups
   matrix[2, 3] median_time_peak; // the effects of the time to detectability, for rows: not exposed tissues (1) or exposed tissue (2); columns: organ groups
    
   // True peak titers
   array[3] real titer_intercept; /// intercept for each organ group, 
   array[3] real titer_location; // the effect of location for each group, 
   matrix[5, 6] titer_route; // the effect of exposure route, where rows: routes; columns: tissues sampled 
   matrix[5, 6] titer_dose; // the effect of dose, where rows: dose in tissue locations; columns: tissues sampled
   array[3] real titer_sex; // the effect of sex for each organ group, 
   matrix[3, 3] titer_age; // the effects of age classes for each organ group, where rows: age classes; columns: organ groups
   matrix[3, 3] titer_species; // the effects of species for each organ group, where rows: species; columns: organ groups
   matrix[L_assay, 3] titer_assay; // the effects of assay type for each organ group, where rows: assays; columns: organ groups
   matrix[2, 3] titer_time; // the effects of the time to peak titer, for rows: not exposed tissues (1) or exposed tissue (2); columns: organ groups
   array[3] real<lower=0> titer_true_sd; // standard deviation around the true peak titer
   
   // Observed peak titers
   matrix[L_lab_titer, 3] titer_lab; // the effects of lab on observed peak titer, where rows: lab; columns: organ groups
   array[3] real<lower=0> titer_obs_intercept_sd; // the intercept for the SD for each organ group
   matrix[L_lab_titer, 3] titer_obs_lab_sd; // lab-specific offsets to the standard deviations for each organ group
   
   // Shape for the time to undetectability
   array[3] real<lower=1> shape_intercept_last;
   
   //  Median for the time to undetectability
   array[3] real median_intercept_last; /// intercept for each organ group
   array[3] real median_location_last; // the effect of location for each group
   matrix[5, 6] median_route_last; // the effect of exposure route, where rows: routes; columns: tissues sampled 
   matrix[5, 6] median_dose_last; // the effect of dose, where rows: dose in tissue locations; columns: tissues sampled
   array[3] real median_sex_last; // the effect of sex for each organ group
   matrix[3, 3] median_age_last; // the effects of age classes for each organ group, where rows: age classes; columns: organ groups
   matrix[3, 3] median_species_last; // the effects of species for each organ group, where rows: species; columns: organ groups
   matrix[L_assay, 3] median_assay_last; // the effects of assay type for each organ group, where rows: assays; columns: organ groups
   matrix[L_lab, 3] median_lab_last; // the effects of lab for each organ group, where rows: lab groups; columns: organ groups
   matrix[2, 3] median_time_last; // the effects of the time to peak titer, for rows: not exposed tissues (1) or exposed tissue (2); columns: organ groups
   
   // True event times & peak titers for each individual, based on bounds
   vector<lower=first_lower_bounds, upper=first_upper_bounds>[N] true_first_wait; // true wait time until detectability
   vector<lower=lower_bound_adjusted(peak_lower_bounds-true_first_wait),
          upper=peak_upper_bounds-true_first_wait>[N] true_peak_wait;   // true wait time until peak titer
   vector<lower=lower_bound_adjusted(last_lower_bounds-true_peak_wait-true_first_wait),
          upper=upper_bound_adjusted(last_upper_bounds-true_peak_wait-true_first_wait)>[N] true_last_wait; // true wait time until undetectability       
   vector<lower=peak_observed_titer-1, upper=14>[N] true_peak_titer_adj; // true peak titers

}

transformed parameters {
  
  // Don't estimate true peak titers for individuals that don't contribute
  //   any information on peak time or peak titer
  vector[N] true_peak_titer = true_peak_titer_adj;
  for (jj in 1:N) {
    if (has_titer[jj] == 0 && has_peak[jj] == 0) {
      true_peak_titer[jj] = 1;
    }
  }
}

model {
  
  // Priors
  // Cure fraction (i.e., percent of individuals that ever test positive)
  percent_intercept ~ normal(0, 3);
  percent_location ~ normal(0, 0.25);
  to_vector(percent_route) ~ normal(0, 0.25);
  to_vector(percent_dose) ~ normal(0, 7*0.25);
  percent_sex ~ normal(0, 0.25);
  to_vector(percent_age) ~ normal(0, 0.25);
  to_vector(percent_species) ~ normal(0, 0.25);
  to_vector(percent_assay) ~ normal(0, 0.25);
  to_vector(percent_lab) ~ normal(0, 0.25);

  // Shape for the time to detectabiltiy 
  to_vector(shape_intercept_first) ~ normal(1, 4);
  
  // Median for the time to detectability
  median_intercept_first ~ normal(-12, 6);
  median_location_first ~ normal(0, 1);
  to_vector(median_route_first) ~ normal(0, 1);
  to_vector(median_dose_first) ~ normal(0, 7 * 0.5);
  median_sex_first ~ normal(0, 1);
  to_vector(median_age_first) ~ normal(0, 1);
  to_vector(median_species_first) ~ normal(0, 1);
  to_vector(median_assay_first) ~ normal(0, 1);
  to_vector(median_lab_first) ~ normal(0, 0.5);
  
  // Shape for the time to peak titer
  to_vector(shape_intercept_peak) ~ normal(1, 4);

  // Median for the time to peak titer
  median_intercept_peak ~ normal(-12, 6);
  median_location_peak ~ normal(0, 1);
  to_vector(median_route_peak) ~ normal(0, 1);
  to_vector(median_dose_peak) ~ normal(0, 7 * 0.5);
  median_sex_peak ~ normal(0, 1);
  to_vector(median_age_peak) ~ normal(0, 1);
  to_vector(median_species_peak) ~ normal(0, 1);
  to_vector(median_assay_peak) ~ normal(0, 1);
  to_vector(median_lab_peak) ~ normal(0, 0.5);
  to_vector(median_time_peak) ~ normal(0, 1);
  
  // True peak titer
  titer_intercept ~ normal(3, 2);
  titer_location  ~ normal(0, 0.25);
  to_vector(titer_route) ~ normal(0, 0.25);
  to_vector(titer_dose) ~ normal(0, 7 * 0.25);
  titer_sex ~ normal(0, 0.25);
  to_vector(titer_age) ~ normal(0, 0.25);
  to_vector(titer_species) ~ normal(0, 0.25);
  to_vector(titer_assay) ~ normal(0, 0.25);
  to_vector(titer_time) ~ normal(0, 0.25);
  to_vector(titer_true_sd) ~ exponential(1);
  
  // Observed vs. true peak titers
  to_vector(titer_lab) ~ normal(0, 0.25);
  titer_obs_intercept_sd ~ normal(-2, 1);
  to_vector(titer_obs_lab_sd) ~ normal(0, 0.1);

  // Shape for the time to undetectability
  to_vector(shape_intercept_last) ~ normal(1, 4);
  
  //  Median for the time to undetectability
  median_intercept_last ~ normal(-12, 6);
  median_location_last ~ normal(0, 1);
  to_vector(median_route_last) ~ normal(0, 1);
  to_vector(median_dose_last) ~ normal(0, 7 * 0.5);
  median_sex_last ~ normal(0, 1);
  to_vector(median_age_last) ~ normal(0, 1);
  to_vector(median_species_last) ~ normal(0, 1);
  to_vector(median_assay_last) ~ normal(0, 1);
  to_vector(median_lab_last) ~ normal(0, 0.5);
  to_vector(median_time_last) ~ normal(0, 1);
  
  
  // Model that loops over each datapoint
  for (ii in 1:N) {
    
    // Determine whether anything needs to be marginalized
    real marginalize = check_for_marginalizing(sex[ii], age[ii],
                                               assay[ii], location[ii]);
    
    // If no need to marginalize, proceed as normal
    if (marginalize == 0) {
      
      // Calculate the percent, shape parameters, and medians based on covariates
      real percent = calculate_percent(percent_intercept[organ_location[ii]],
                                       percent_location[organ_location[ii]] * location[ii],
                                       percent_route[route[ii], tissue_location[ii]],
                                       percent_dose[1, tissue_location[ii]] * dose_nose[ii], 
                                       percent_dose[2, tissue_location[ii]] * dose_throat[ii],
                                       percent_dose[3, tissue_location[ii]] * dose_trachea[ii],
                                       percent_dose[4, tissue_location[ii]] * dose_lung[ii],
                                       percent_dose[5, tissue_location[ii]] * dose_gi[ii],
                                       percent_sex[organ_location[ii]] * sex[ii],
                                       percent_age[age[ii], organ_location[ii]],
                                       percent_species[species[ii], organ_location[ii]],
                                       percent_assay[assay[ii], organ_location[ii]],
                                       percent_lab[lab[ii], organ_location[ii]]
                                       );
                                       
      real shape_first = shape_intercept_first[organ_location[ii]];
                                         
      real median_first = calculate_median(median_intercept_first[organ_location[ii]],
                                           median_location_first[organ_location[ii]] * location[ii], 
                                           median_route_first[route[ii], tissue_location[ii]],
                                           median_dose_first[1, tissue_location[ii]] * dose_nose[ii], 
                                           median_dose_first[2, tissue_location[ii]] * dose_throat[ii],
                                           median_dose_first[3, tissue_location[ii]] * dose_trachea[ii],
                                           median_dose_first[4, tissue_location[ii]] * dose_lung[ii],
                                           median_dose_first[5, tissue_location[ii]] * dose_gi[ii],
                                           median_sex_first[organ_location[ii]] * sex[ii],
                                           median_age_first[age[ii], organ_location[ii]],
                                           median_species_first[species[ii], organ_location[ii]],
                                           median_assay_first[assay[ii], organ_location[ii]],
                                           median_lab_first[lab[ii], organ_location[ii]],
                                           0 // Time to detectability doesn't depend on a previous event time
                                           );  
                                           
                              
      // Calculate the log likelihood
      real log_lik_first = calculate_log_lik(ever_positive[ii], 
                                             1, 
                                             true_first_wait[ii],
                                             first_lower_bounds[ii],
                                             percent, 
                                             shape_first, 
                                             median_first);
      target += log_lik_first; // Increment the log probability
      
      // Check whether this datapoint contributes information on the time to peak titers
      if (has_peak[ii] == 1) {
        
        // Calculate shape and median based on covariates
        real shape_peak = shape_intercept_peak[organ_location[ii]];
        real median_peak = calculate_median(median_intercept_peak[organ_location[ii]],
                                            median_location_peak[organ_location[ii]] * location[ii], //0,
                                            median_route_peak[route[ii], tissue_location[ii]],
                                            median_dose_peak[1, tissue_location[ii]] * dose_nose[ii], 
                                            median_dose_peak[2, tissue_location[ii]] * dose_throat[ii],
                                            median_dose_peak[3, tissue_location[ii]] * dose_trachea[ii],
                                            median_dose_peak[4, tissue_location[ii]] * dose_lung[ii],
                                            median_dose_peak[5, tissue_location[ii]] * dose_gi[ii],
                                            median_sex_peak[organ_location[ii]] * sex[ii],
                                            median_age_peak[age[ii], organ_location[ii]],
                                            median_species_peak[species[ii], organ_location[ii]],
                                            median_assay_peak[assay[ii], organ_location[ii]],
                                            median_lab_peak[lab[ii], organ_location[ii]],  
                                            median_time_peak[inoc[ii] + 1, organ_location[ii]] * true_first_wait[ii]);  
                                            
        // Compute the log likelihood      
        real log_lik_peak = calculate_log_lik(1, 
                                              2, 
                                              true_peak_wait[ii],
                                              -999999, 
                                              -999999,
                                              shape_peak, 
                                              median_peak);
        
        // Increment the log likelihood                                      
        target += log_lik_peak;
        
        // Check if the datapoint has information on the peak titer
        if (has_titer[ii] == 1) {
          
          // Increment the log likelihood based on the covariates
          target += normal_lpdf(true_peak_titer[ii] | 
                                titer_intercept[organ_location[ii]] +
                                titer_location[organ_location[ii]] * location[ii] + 
                                titer_route[route[ii], tissue_location[ii]] +
                                titer_dose[1, tissue_location[ii]] * dose_nose[ii] +
                                titer_dose[2, tissue_location[ii]] * dose_throat[ii] +
                                titer_dose[3, tissue_location[ii]] * dose_trachea[ii] +
                                titer_dose[4, tissue_location[ii]] * dose_lung[ii] +
                                titer_dose[5, tissue_location[ii]] * dose_gi[ii] +
                                titer_sex[organ_location[ii]] * sex[ii] +
                                titer_age[age[ii], organ_location[ii]] +
                                titer_species[species[ii], organ_location[ii]] +
                                titer_assay[assay[ii], organ_location[ii]] +
                                titer_time[inoc[ii] + 1, organ_location[ii]] * true_peak_wait[ii],
                                titer_true_sd[organ_location[ii]]);
          target += normal_lpdf(peak_observed_titer[ii] |
                                true_peak_titer[ii] + titer_lab[lab[ii], organ_location[ii]], 
                                exp(titer_obs_intercept_sd[organ_location[ii]] +
                                    titer_obs_lab_sd[lab[ii], organ_location[ii]]));
        }
        
        // Estimate a true peak titer for individuals with peak time information
        //    but no peak titer information
        else if (has_titer[ii] == 0) {
          true_peak_titer[ii] ~ normal(titer_intercept[organ_location[ii]] +
                                       titer_location[organ_location[ii]] * location[ii] + 
                                       titer_route[route[ii], tissue_location[ii]] +
                                       titer_dose[1, tissue_location[ii]] * dose_nose[ii] +
                                       titer_dose[2, tissue_location[ii]] * dose_throat[ii] +
                                       titer_dose[3, tissue_location[ii]] * dose_trachea[ii] +
                                       titer_dose[4, tissue_location[ii]] * dose_lung[ii] +
                                       titer_dose[5, tissue_location[ii]] * dose_gi[ii] +
                                       titer_sex[organ_location[ii]] * sex[ii] +
                                       titer_age[age[ii], organ_location[ii]] +
                                       titer_species[species[ii], organ_location[ii]] +
                                       titer_assay[assay[ii], organ_location[ii]] +
                                       titer_time[inoc[ii] + 1, organ_location[ii]] * true_peak_wait[ii], 
                                       titer_true_sd[organ_location[ii]]);
        }
      }
      
     // Check whether each datapoint contributes information on the time to undetectability
      if (has_last_positive[ii] == 1) {
       
       // Calculate the shape and median based on covariates
       real shape_last = shape_intercept_last[organ_location[ii]];
       
       real median_last = calculate_median(median_intercept_last[organ_location[ii]],
                                           median_location_last[organ_location[ii]] * location[ii],
                                           median_route_last[route[ii], tissue_location[ii]],
                                           median_dose_last[1, tissue_location[ii]] * dose_nose[ii], 
                                           median_dose_last[2, tissue_location[ii]] * dose_throat[ii],
                                           median_dose_last[3, tissue_location[ii]] * dose_trachea[ii],
                                           median_dose_last[4, tissue_location[ii]] * dose_lung[ii],
                                           median_dose_last[5, tissue_location[ii]] * dose_gi[ii],
                                           median_sex_last[organ_location[ii]] * sex[ii],
                                           median_age_last[age[ii], organ_location[ii]],
                                           median_species_last[species[ii], organ_location[ii]],
                                           median_assay_last[assay[ii], organ_location[ii]],
                                           median_lab_last[lab[ii], organ_location[ii]],
                                           median_time_last[inoc[ii] + 1, organ_location[ii]] * true_peak_titer[ii]
                                           ); 
         // Compute the log likelihood
         real log_lik_last = calculate_log_lik(1, 
                                               2, 
                                               true_last_wait[ii],
                                               -999999, 
                                               -999999, 
                                               shape_last, 
                                               median_last);
                                               
         // Incremenet the log likelihood                                      
         target += log_lik_last;                                  
       
      }
    }
    
    // If marginalizing is necessary, then must integrate over all possibilities
    else if (marginalize == 1) {
      
      // Determine which of the covariates need to be marginalized over
      array[2] int location_options = get_cofactor_options(location[ii], {0, 1});
      array[2] int sex_options = get_cofactor_options(sex[ii], {0, 1});
      array[3] int age_options;
      array[3] int assay_options = get_cofactor_options(assay[ii], {1, 2, 3});
      
      // Set the options depending on whether geriatrics is possible or not
      if (prob_geriatric[ii] == 0) {
        age_options = get_cofactor_options(age[ii], {1, 2, 2});
      }
      else {
        age_options = get_cofactor_options(age[ii], {1, 2, 3});
      }
      
      // Create vectors to store the marginalized probabilities
      int num_marginalize_options = num_elements(get_unique_vec(sex_options)) *
                                    num_elements(get_unique_vec(age_options)) * 
                                    num_elements(get_unique_vec(location_options)) *
                                    num_elements(get_unique_vec(assay_options));
      array[num_marginalize_options] real first_log_liks = rep_array(0, num_marginalize_options);
      array[num_marginalize_options] real peak_log_liks = rep_array(0, num_marginalize_options);
      array[num_marginalize_options] real last_log_liks = rep_array(0, num_marginalize_options);
      array[num_marginalize_options] real titer_log_liks = rep_array(0, num_marginalize_options);
      array[num_marginalize_options] real titer_obs_log_liks = rep_array(0, num_marginalize_options);
      
      // Loop over all possible options -- note that covariates not needing to be marginalized will have only one entry to be looped over
      int ii_marginalize = 0;
      for (assay_ii in get_unique_vec(assay_options)) {
        for (sex_ii in get_unique_vec(sex_options)) {
          for (age_ii in get_unique_vec(age_options)) {
            for (loc_ii in get_unique_vec(location_options)) {
              
              // Incremement index for marginalizing
              ii_marginalize = ii_marginalize + 1;
              int assay_ii_adjusted = assay_ii;
              
              // Figure out what tissue location integrating over
              int tissue_ii = 0;
              if (organ_location[ii] == 1 && loc_ii == 0) {
                tissue_ii = 1;
              }
              else if (organ_location[ii] == 1 && loc_ii == 1) {
                tissue_ii = 2;
              }
              else if (organ_location[ii] == 2 && loc_ii == 0) {
                tissue_ii = 3;
              }
              else if (organ_location[ii] == 2 && loc_ii == 1) {
                tissue_ii = 4;
              }
              else if (organ_location[ii] == 3 && loc_ii == 0) {
                tissue_ii = 5;
              }
              else if (organ_location[ii] == 3 && loc_ii == 1) {
                tissue_ii = 6;
              }
              
              // Get sex-specific probabilities
              real sex_prob;
              if (sex_ii == 0) {
                sex_prob = prob_female[ii];
              }
              else if (sex_ii == 1) {
                sex_prob = prob_male[ii];
              }
              
              // Set the probabilities based on the total number of options
              real log_assay_prob = get_cofactor_probs(assay[ii], 1.0 / 3); 
              real log_sex_prob = get_cofactor_probs(sex[ii], sex_prob); 
              real log_age_prob = get_cofactor_probs(age[ii], 1.0 / num_elements(get_unique_vec(age_options))); 
              real log_loc_prob = get_cofactor_probs(location[ii], 1.0 / 2); 
            
              // Calculate the percent, shape, and median for the time to detectability
              real percent = calculate_percent(percent_intercept[organ_location[ii]],
                                               percent_location[organ_location[ii]] * loc_ii,
                                               percent_route[route[ii], tissue_ii],
                                               percent_dose[1, tissue_ii] * dose_nose[ii], 
                                               percent_dose[2, tissue_ii] * dose_throat[ii],
                                               percent_dose[3, tissue_ii] * dose_trachea[ii],
                                               percent_dose[4, tissue_ii] * dose_lung[ii],
                                               percent_dose[5, tissue_ii] * dose_gi[ii],
                                               percent_sex[organ_location[ii]] * sex_ii,
                                               percent_age[age_ii, organ_location[ii]],
                                               percent_species[species[ii], organ_location[ii]],
                                               percent_assay[assay_ii, organ_location[ii]],
                                               percent_lab[lab[ii], organ_location[ii]]
                                               );
              
              real shape_first = shape_intercept_first[organ_location[ii]];
                                                 
              real median_first = calculate_median(median_intercept_first[organ_location[ii]],
                                                   median_location_first[organ_location[ii]] * loc_ii,
                                                   median_route_first[route[ii], tissue_ii],
                                                   median_dose_first[1, tissue_ii] * dose_nose[ii], 
                                                   median_dose_first[2, tissue_ii] * dose_throat[ii],
                                                   median_dose_first[3, tissue_ii] * dose_trachea[ii],
                                                   median_dose_first[4, tissue_ii] * dose_lung[ii],
                                                   median_dose_first[5, tissue_ii] * dose_gi[ii],
                                                   median_sex_first[organ_location[ii]] * sex_ii,
                                                   median_age_first[age_ii, organ_location[ii]],
                                                   median_species_first[species[ii], organ_location[ii]],
                                                   median_assay_first[assay_ii, organ_location[ii]],
                                                   median_lab_first[lab[ii], organ_location[ii]],
                                                   0
                                                   );
                                                   
              // Calculate the log likelihood
              real log_lik_first = calculate_log_lik(ever_positive[ii], 
                                                     1, 
                                                     true_first_wait[ii],
                                                     first_lower_bounds[ii],
                                                     percent, 
                                                     shape_first, 
                                                     median_first);
              
              // Combine with the probabilities for the different covariates
              first_log_liks[ii_marginalize] = log_sex_prob +
                                               log_age_prob + 
                                               log_assay_prob +
                                               log_loc_prob +
                                               log_lik_first; 
               
              // Check whether this datapoint contributes information on the time to peak titers
              if (has_peak[ii] == 1) {
                
                real shape_peak = shape_intercept_peak[organ_location[ii]]; 
                real median_peak = calculate_median(median_intercept_peak[organ_location[ii]],
                                                    median_location_peak[organ_location[ii]] * loc_ii,
                                                    median_route_peak[route[ii], tissue_ii],
                                                    median_dose_peak[1, tissue_ii] * dose_nose[ii], 
                                                    median_dose_peak[2, tissue_ii] * dose_throat[ii],
                                                    median_dose_peak[3, tissue_ii] * dose_trachea[ii],
                                                    median_dose_peak[4, tissue_ii] * dose_lung[ii],
                                                    median_dose_peak[5, tissue_ii] * dose_gi[ii],
                                                    median_sex_peak[organ_location[ii]] * sex_ii,
                                                    median_age_peak[age_ii, organ_location[ii]],
                                                    median_species_peak[species[ii], organ_location[ii]],
                                                    median_assay_peak[assay_ii, organ_location[ii]],
                                                    median_lab_peak[lab[ii], organ_location[ii]],
                                                    median_time_peak[inoc[ii] + 1, organ_location[ii]] * true_first_wait[ii]);
                
                // Calculate the log likelihood                                   
                real log_lik_peak = calculate_log_lik(1, 
                                                      2, 
                                                      true_peak_wait[ii],
                                                      -999999, 
                                                      -999999,
                                                      shape_peak, 
                                                      median_peak);
                peak_log_liks[ii_marginalize] = log_sex_prob +
                                                log_age_prob + 
                                                log_assay_prob +
                                                log_loc_prob +
                                                log_lik_peak; 
                                                
                // Check if the datapoint has information on the peak titer                                
                if (has_titer[ii] == 1) {
          
                    titer_log_liks[ii_marginalize] = 
                              normal_lpdf(true_peak_titer[ii] | 
                                          titer_intercept[organ_location[ii]] +
                                          titer_location[organ_location[ii]] * loc_ii + 
                                          titer_route[route[ii], tissue_ii] +
                                          titer_dose[1, tissue_ii] * dose_nose[ii] +
                                          titer_dose[2, tissue_ii] * dose_throat[ii] +
                                          titer_dose[3, tissue_ii] * dose_trachea[ii] +
                                          titer_dose[4, tissue_ii] * dose_lung[ii] +
                                          titer_dose[5, tissue_ii] * dose_gi[ii] +
                                          titer_sex[organ_location[ii]] * sex_ii +
                                          titer_age[age_ii, organ_location[ii]] +
                                          titer_species[species[ii], organ_location[ii]] +
                                          titer_assay[assay_ii, organ_location[ii]]  +
                                          titer_time[inoc[ii] + 1, organ_location[ii]] * true_peak_wait[ii],
                                          titer_true_sd[organ_location[ii]]);
                                          
                    titer_obs_log_liks[ii_marginalize] = 
                    normal_lpdf(peak_observed_titer[ii] |
                                          true_peak_titer[ii] +
                                          titer_lab[lab[ii], organ_location[ii]],  
                                          exp(titer_obs_intercept_sd[organ_location[ii]] +
                                              titer_obs_lab_sd[lab[ii], organ_location[ii]]));
                  }  
                  
                else if (has_titer[ii] == 0) {
                  
                  titer_obs_log_liks[ii_marginalize] = normal_lpdf(true_peak_titer[ii] | 
                                               titer_intercept[organ_location[ii]] +
                                               titer_location[organ_location[ii]] * loc_ii + 
                                               titer_route[route[ii], tissue_ii] +
                                               titer_dose[1, tissue_ii] * dose_nose[ii] +
                                               titer_dose[2, tissue_ii] * dose_throat[ii] +
                                               titer_dose[3, tissue_ii] * dose_trachea[ii] +
                                               titer_dose[4, tissue_ii] * dose_lung[ii] +
                                               titer_dose[5, tissue_ii] * dose_gi[ii] +
                                               titer_sex[organ_location[ii]] * sex_ii +
                                               titer_age[age_ii, organ_location[ii]] +
                                               titer_species[species[ii], organ_location[ii]] +
                                               titer_assay[assay_ii, organ_location[ii]] +
                                               titer_time[inoc[ii] + 1, organ_location[ii]] * true_peak_wait[ii],
                                               titer_true_sd[organ_location[ii]]);
                  
                }
              } 
              
              // Check whether this datapoint contributes information on the time to undetectability
              if (has_last_positive[ii] == 1) {
                                                
                real shape_last = shape_intercept_last[organ_location[ii]];
                
                real median_last = calculate_median(median_intercept_last[organ_location[ii]],
                                                    median_location_last[organ_location[ii]] * loc_ii,
                                                    median_route_last[route[ii], tissue_ii],
                                                    median_dose_last[1, tissue_ii] * dose_nose[ii], 
                                                    median_dose_last[2, tissue_ii] * dose_throat[ii],
                                                    median_dose_last[3, tissue_ii] * dose_trachea[ii],
                                                    median_dose_last[4, tissue_ii] * dose_lung[ii],
                                                    median_dose_last[5, tissue_ii] * dose_gi[ii],
                                                    median_sex_last[organ_location[ii]] * sex_ii,
                                                    median_age_last[age_ii, organ_location[ii]],
                                                    median_species_last[species[ii], organ_location[ii]],
                                                    median_assay_last[assay_ii, organ_location[ii]],
                                                    median_lab_last[lab[ii], organ_location[ii]],
                                                    median_time_last[inoc[ii] + 1, organ_location[ii]] * true_peak_titer[ii]
                                                    );
                
                real log_lik_last = 0;
                log_lik_last = calculate_log_lik(1, 
                                                   2, 
                                                   true_last_wait[ii],
                                                   -999999, 
                                                   -999999, 
                                                   shape_last, 
                                                   median_last);
                                                      
                last_log_liks[ii_marginalize] = log_sex_prob +
                                                log_age_prob + 
                                                log_assay_prob +
                                                log_loc_prob +
                                                log_lik_last; 
                                                
              }
            }
          }
        }
      }

     // Incremement based on all marginalized probabilities
     target += log_sum_exp(first_log_liks);
     
     if (has_peak[ii] == 1) {
       target += log_sum_exp(peak_log_liks);
     }
     if (has_titer[ii] == 1) {
       target += log_sum_exp(titer_log_liks);
       target += log_sum_exp(titer_obs_log_liks);
     }
     if (has_last_positive[ii] == 1 && last_positive_type[ii] == 1) {
       target += log_sum_exp(last_log_liks);
     }
    }
  }
}

