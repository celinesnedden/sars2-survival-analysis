# Run this before everything else to set up the environment & load necessary functions

# Set working directory --------------------------------------------------------

setwd("./submission") # Adjust as necessary to your directory; all other paths will be relative


# Install / load packages needed throughout analysis ---------------------------

packages <- c("dplyr", "cmdstanr", "tidyr", "ggridges", "extraDistr", "crayon", 
              "ggplot2", "stringr", "wesanderson", "patchwork", "ggpubr",
              "foreach", "parallel", "doParallel", "cowplot", "png", "grid",
              "posterior", "data.table", "tidyverse")

missing <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(missing) > 0) install.packages(missing)

invisible(lapply(packages, library, character.only = TRUE))


# Load all functions -----------------------------------------------------------

# This contains the custom functions written for this analysis
source("./code/all-custom-functions.R")


# Set color palette -----------------------------------------------------------

df.color <- assign_colors()
df.colors <- df.color


# Set the max doses by location for rescaling ----------------------------------

#  Note these were computed and taken from the 03-model-fitting file
max_dose_nose <- 6.821186 
max_dose_throat <- 6.90768 
max_dose_trachea <- 7.350608
max_dose_lung <- 6.562293
max_dose_gi <- 7


