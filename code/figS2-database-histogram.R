# This file: - Generates histograms for the data in the full database

# Prep -------------------------------------------------------------------------

# Load the full database
dat.full <- read.csv("./data/database-clean.csv")
dat.full <- assign_all_names(dat.full)


# Plot  ------------------------------------------------------------------------

# Generate each panel individually
figS2.route <- plot_data_histograms(dat.full, "route_full"); figS2.route
figS2.dose <- plot_data_histograms(dat.full, "dose"); figS2.dose
figS2.sex <- plot_data_histograms(dat.full, "sex"); figS2.sex
figS2.age <- plot_data_histograms(dat.full, "age"); figS2.age
figS2.species <- plot_data_histograms(dat.full, "species"); figS2.species
figS2.assay <- plot_data_histograms(dat.full, "assay_full"); figS2.assay
figS2.location <- plot_data_histograms(dat.full, "location"); figS2.location

# Combine into one figure
figS2 <- figS2.dose + 
  figS2.route +  
  figS2.assay +  
  figS2.age + 
  figS2.sex + 
  figS2.species + 
  figS2.location + 
  plot_layout(nrow = 1, 
              widths = c(8, 10, 5, 4.2, 3.2, 3.2, 4.2)); figS2


# Save  ------------------------------------------------------------------------

ggsave("./outputs/figures/figS2-database-histrogram.png",
       plot = figS2,
       width = 8, 
       height = 1.8,
       dpi = 600)
