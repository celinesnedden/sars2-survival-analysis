# This file: - Creates the components of Figure 1 that are possible to create in R

# 1A: Data distribution --------------------------------------------------------

## Prep ------------------------------------------------------------------------

# Load the event time database
dat.event <- read.csv("./data/df-event-times.csv")
dat.event <- assign_all_names(dat.event) 


## Plot  -----------------------------------------------------------------------

# Generate each panel individually
fig1A.route <- plot_data_histograms(dat.event, "route"); fig1A.route
fig1A.dose <- plot_data_histograms(dat.event, "dose"); fig1A.dose
fig1A.sex <- plot_data_histograms(dat.event, "sex"); fig1A.sex
fig1A.age <- plot_data_histograms(dat.event, "age"); fig1A.age
fig1A.species <- plot_data_histograms(dat.event, "species"); fig1A.species
fig1A.assay <- plot_data_histograms(dat.event, "assay"); fig1A.assay
fig1A.location <- plot_data_histograms(dat.event, "location"); fig1A.location

# Combine into one figure
fig1A <- fig1A.dose +  
  fig1A.route +  
  fig1A.age +  
  fig1A.sex +  
  fig1A.species + 
  fig1A.assay + 
  fig1A.location + 
  plot_layout(nrow = 1, 
              widths = c(7, 5.2, 4.2, 3.2, 3.2, 2.2, 3.2)); fig1A


## Save  ------------------------------------------------------------------------

ggsave("./outputs/figures/fig1a-event-time-histogram.pdf",
       plot = fig1A,
       width = 6, 
       height = 1.5,
       dpi = 600)



# 1B: Sample Size Table --------------------------------------------------------

# Made in Powerpoint


# 1C: Example Censoreds --------------------------------------------------------

# Load full & event time database
dat.event <- read.csv("./data/df-event-times.csv")
dat.full <- read.csv("./data/df-all-times.csv")

# For best visualization, include only individs where the peak lower bound
#   is not the same as the first positive lower bound
dat.event <- subset(dat.event, first_lower_bound != peak_lower_bound &
                      peak_observed_time != last_lower_bound & 
                      last_lower_bound < 20 & 
                      (last_upper_bound == Inf | last_upper_bound < 20))


# Randomly sample an individual for each location
set.seed(30)
num_samps <- 5
random_indivs_urt <- sample(dat.event$indiv_sample[dat.event$location_grp == "URT"], size = num_samps)
random_indivs_lrt <- sample(dat.event$indiv_sample[dat.event$location_grp == "LRT"], size = num_samps)
random_indivs_gi <- sample(dat.event$indiv_sample[dat.event$location_grp == "GI"], size = num_samps)


# Subset event times & all times to these individuals & plottable DPI range
dat.event <- subset(dat.event, indiv_sample %in% c(random_indivs_urt, 
                                                   random_indivs_lrt,
                                                   random_indivs_gi))
dat.full <- subset(dat.full, indiv_sample %in% unique(dat.event$indiv_sample) &
                     day_post_infection <= 20)


# For best visualization, order them by bounds and set as factor
dat.event <- dat.event[order(dat.event$first_lower_bound,
                             dat.event$first_upper_bound,
                             dat.event$peak_observed_time,
                             dat.event$last_lower_bound,
                             dat.event$last_upper_bound), ]
dat.full$indiv_sample <- factor(dat.full$indiv_sample, levels = dat.event$indiv_sample)
dat.event$indiv_sample <- factor(dat.event$indiv_sample, levels = dat.event$indiv_sample)


# Classify censor type for the last positive
dat.event$last_censor_type[dat.event$last_upper_bound == Inf] <- "RIGHT"
dat.event$last_censor_type[dat.event$last_upper_bound != Inf] <- "INT"



## Plot -----------------------------------------------------------------------

fig1c <- ggplot(dat.event) + 
  
  # Colored lines connecting event times
  geom_segment(aes(x = -0.5,
                   xend = 17.5,
                   y = indiv_sample,
                   yend = indiv_sample), 
               color = "grey77",
               alpha = 1,
               linewidth = 0.03) +
  geom_segment(aes(x = as.numeric(first_lower_bound),
                   xend = as.numeric(first_upper_bound),
                   y = indiv_sample,
                   yend = indiv_sample), 
               color = "#4F609C",
               alpha = 1,
               linewidth = 0.5) +
  geom_segment(aes(x = as.numeric(peak_lower_bound),
                   xend = as.numeric(peak_upper_bound),
                   y = indiv_sample,
                   yend = indiv_sample), 
               color =  "#9E4472",
               alpha = 1,
               linewidth = 0.5) +
  geom_segment(aes(x = as.numeric(last_lower_bound),
                   xend = as.numeric(last_upper_bound),
                   y = indiv_sample,
                   yend = indiv_sample,
                   linetype = last_censor_type), 
               color = "#628D56",
               alpha = 1,
               linewidth = 0.5) +
  
  # All sample times
  geom_point(aes(x = 0,
                 y = indiv_sample,
                 fill = "Sample Day"),
             shape = 21, size = 0.5) +
  geom_point(data = dat.full, 
             aes(x = as.numeric(day_post_infection),
                 y = indiv_sample,
                 fill = "Sample Day"),
             shape = 21, size = 0.5) +
  
  # Event times
  
  ## First positive
  geom_point(aes(x = as.numeric(first_upper_bound),
                 y = indiv_sample, 
                 fill = "First Positive",
                 shape = "First Positive"),
             alpha = 1,
             size = 2.5) +
  
  ## Peak time
  geom_point(aes(x = as.numeric(peak_observed_time),
                 y = indiv_sample,
                 fill = "Peak",
                 shape = "Peak"),
             #shape = 22, 
             size = 2.5) +
  
  ## Last Positive
  geom_point(aes(x = as.numeric(last_lower_bound),
                 y = indiv_sample, 
                 fill = "Last Positive",
                 shape = "Last Positive"),
             alpha = 1,
             size = 2.5) +
  
  
  # Manual aesthetics
  scale_shape_manual(values = c("First Positive" = 23,
                                "Peak" = 24,
                                "Last Positive" = 22,
                                "Sample Day" = 21),
                     guide = "none") +
  scale_fill_manual(values = c("First Positive" = "#4F609C",
                               "Peak" = "#9E4472",
                               "Last Positive" = "#628D56",
                               "Sample Day" = "black"),
                    breaks = c("First Positive",
                               "Peak", 
                               "Last Positive",
                               "Sample Day")) +
  scale_linetype_manual(values = c("INT" = "solid",
                                   "RIGHT" = "dashed"),
                        guide = "none") +
  
  labs(y = "Individual", x = "Days post infection",
       fill = "") +
  guides(fill = guide_legend(override.aes = list(shape = c(23, 24, 22, 21),
                                                 size = c(1.5, 1.5, 1.5, 0.75 * 0.5)),
                             byrow = TRUE,
                             keywidth = 0.2,
                             keyheight = 0.15,
                             default.unit="inch")) +
  scale_x_continuous(limits = c(-0.5, 20),
                     breaks = seq(0, 20, 4),
                     expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0.05, 0.15)) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none",
        text = element_text(size = 9),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        strip.background  = element_rect(colour = "transparent", fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_rect(fill = "white",
                                         colour = "white",
                                         size = 0.2, linetype = "solid"),
        legend.margin = margin(c(5, 5, 5, 2)),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA,
                                    colour = "black",
                                    linewidth = 0.25),
        panel.grid.major = element_line(linewidth = 0.1, color = "grey88"), 
        panel.grid.minor = element_line(linewidth = 0.05, color = "grey88")); fig1c


## Save -------------------------------------------------------------------------

ggsave('./outputs/figures/fig1c-example-censored-data.pdf',
       plot = fig1c,
       width = 1.85, 
       height = 2.25,
       dpi = 600)


# 1D-F: Diagram ----------------------------------------------------------------

# This was made in BioRender


