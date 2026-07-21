# This file: - Generates adjacency matrices according to single-tissue inoculations
#            - Uses the matrices to estimate a tissue connectivity structure


# Prep -------------------------------------------------------------------------

# Load the model
fit <- readRDS("./outputs/fits/fit-main.RDS")

# Get the parameters & samples to extract for all routes
sample.pars <- get_all_parameter_samples(fit)
set.seed(12378965)
sample_nums <- base::sample(1:nrow(sample.pars), size = 1000, replace = TRUE) 

# Set maximum location-specific doses for rescaling
#  (note these were computed and taken from the 03-model-fitting file)
max_dose_nose <- 6.821186 
max_dose_throat <- 6.90768 
max_dose_trachea <- 7.350608
max_dose_lung <- 6.562293
max_dose_gi <- 7

# Set covariates
dose <- 7
assay <- 1
age <- c(1, 2, 3)
sex <- c(0, 1)
species <- c(1, 2, 3)
lab <- 1 # note this is not used, but needs to be passed through the functions


# Generate predictions ---------------------------------------------------------

# Generate predictions for each tissue individually
df.nose <- generate_single_tissue_inoculations(sample.pars,
                                               sample_nums,
                                               tissue_inoculated = "Nose", 
                                               dose_total = dose / max_dose_nose, 
                                               assay_idx = assay, 
                                               age_idx = age, 
                                               sex_idx = sex, 
                                               species_idx = species, 
                                               lab_idx = lab, 
                                               lab_effect = "No")


df.throat <- generate_single_tissue_inoculations(sample.pars,
                                                 sample_nums,
                                                 tissue_inoculated = "Throat", 
                                                 dose_total = dose / max_dose_throat, 
                                                 assay_idx = assay, 
                                                 age_idx = age, 
                                                 sex_idx = sex, 
                                                 species_idx = species, 
                                                 lab_idx = lab, 
                                                 lab_effect = "No")


df.trachea <- generate_single_tissue_inoculations(sample.pars,
                                                  sample_nums,
                                                  tissue_inoculated = "Trachea", 
                                                  dose_total = dose / max_dose_trachea, 
                                                  assay_idx = assay, 
                                                  age_idx = age, 
                                                  sex_idx = sex, 
                                                  species_idx = species, 
                                                  lab_idx = lab, 
                                                  lab_effect = "No")


df.lung <- generate_single_tissue_inoculations(sample.pars,
                                               sample_nums,
                                               tissue_inoculated = "Lung", 
                                               dose_total = dose / max_dose_lung, 
                                               assay_idx = assay, 
                                               age_idx = age, 
                                               sex_idx = sex, 
                                               species_idx = species, 
                                               lab_idx = lab, 
                                               lab_effect = "No")


df.gi <- generate_single_tissue_inoculations(sample.pars,
                                             sample_nums,
                                             tissue_inoculated = "Upper GI", 
                                             dose_total = dose / max_dose_gi, 
                                             assay_idx = assay, 
                                             age_idx = age, 
                                             sex_idx = sex, 
                                             species_idx = species, 
                                             lab_idx = lab, 
                                             lab_effect = "No")

# Combine all the tissues
df.all <- rbind(df.nose, df.throat, df.trachea, df.lung, df.gi)

# Factor the tissue inoculated column
df.all$tissue_inoculated <- factor(df.all$tissue_inoculated,
                                   levels = c("Nose", "Throat", "Trachea", 
                                              "Lung", "Upper GI", "Lower GI"))


# Generate adjacency matrices --------------------------------------------------

# Take the floor of all first detectability times
df.all$time_floor <- floor(df.all$median_first)

# Set all inoculated tissues as day 0
df.all$time_floor[df.all$tissue_name == df.all$tissue_inoculated] <- 0

# Randomly sample whether a sample results in the tissue testing positive, based on 
#    the predicted probability of positivity
for (row_num in 1:nrow(df.all)) {
  ever_pos <- rbinom(1, 1, prob = df.all$percent[row_num])
  df.all$ever_positive[row_num] <- ever_pos
}

# Assign names for the tissue being sampled
df.all <- assign_tissue_names(df.all)

# Get the adjacency matrices for each tissue
df.adj <- subset(df.all, ever_positive == 1)
df.adj.nose <- generate_adjacency_matrix(df.adj, "Nose")
df.adj.throat <- generate_adjacency_matrix(df.adj, "Throat")
df.adj.trachea <- generate_adjacency_matrix(df.adj, "Trachea")
df.adj.lung <- generate_adjacency_matrix(df.adj, "Lung")
df.adj.gi <- generate_adjacency_matrix(df.adj, "Upper GI")


# Figure S25 --------------------------------------------------------------------

df.adj.nose$tissue_inoculated <- "Nose"
df.adj.throat$tissue_inoculated <- "Throat"
df.adj.trachea$tissue_inoculated <- "Trachea"
df.adj.lung$tissue_inoculated <- "Lung"
df.adj.gi$tissue_inoculated <- "GI"

df.adj.all <- rbind(df.adj.nose, df.adj.throat, 
                    df.adj.trachea, df.adj.lung, df.adj.gi)


df.adj.adj <- df.adj.all %>%
  rownames_to_column(var = "tissueY") %>%
  pivot_longer(
    cols = c(Nose, Throat, Trachea, Lung, Upper.GI, Lower.GI),
    names_to = "tissueX",
    values_to = "value") %>%
  mutate(tissueY = str_remove(tissueY, "\\d+$")) %>%
  mutate(tissueX = str_replace_all(tissueX, "\\.", " "))

df.adj.adj$tissue_inoculated[df.adj.adj$tissue_inoculated == "GI"] <- "Upper GI" 

df.adj.adj$tissue_inoculated <- factor(df.adj.adj$tissue_inoculated, 
                                       levels = c("Nose", "Throat", "Trachea", 
                                                  "Lung", "Lower GI", "Upper GI"))  
df.adj.adj$tissueX <- factor(df.adj.adj$tissueX, 
                             levels = c("Nose", "Throat", "Trachea", "Lung", 
                                        "Upper GI", "Lower GI"))  
df.adj.adj$tissueY <- factor(df.adj.adj$tissueY, 
                             levels = rev(c("Nose", "Throat", "Trachea", 
                                            "Lung", "Upper GI", "Lower GI")))

## Plot ------------------------------------------------------------------------

figS25 <- ggplot(df.adj.adj, aes(x = tissueX, y = tissueY)) + 
  geom_tile(aes(fill = value), color = "black", linewidth = 0.5) +
  geom_tile(data = subset(df.adj.adj, tissueX == tissueY), 
            fill = "black", color = "black", linewidth = 0.5) +
  geom_tile(data = subset(df.adj.adj, tissueX == tissue_inoculated), 
            fill = "black", color = "black", linewidth = 0.5) +
  geom_tile(data = subset(df.adj.adj, tissueY == tissue_inoculated &
                            tissueY != tissueX), 
             fill = "transparent", color = "red", linewidth = 2) +
  geom_tile(data = subset(df.adj.adj, tissueX != tissue_inoculated &
                            tissueX == "Nose" & tissueY != "Nose"), 
            fill = "transparent", color = "orange", linewidth = 1,
            alpha = 0.5) +
  geom_text(aes(label = value), size = 4) +
  facet_wrap(.~ tissue_inoculated, nrow = 2) +
  coord_cartesian(expand = FALSE) +
  labs(x = "Tissue To", y = "Tissue From", 
       fill = "Number of\nsubsequent\npositives") +
  scale_fill_gradient(low  = "white", high = df.color$Throat, breaks = seq(0, 16000, 8000)) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill ="grey", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black")); figS25


## Save ------------------------------------------------------------------------

ggsave("./outputs/figures/figS25-adjacency-matrices.png",
       plot = figS25,
       width = 9, 
       height = 6,
       dpi = 600)


# Figure 5A  -------------------------------------------------------------------

## Calculate outflows ----------------------------------------------------------

# Combine the rows of each matrix 
df.adj.outflows <- rbind(df.adj.nose[1, ],
                         df.adj.throat[2, ],
                         df.adj.trachea[3, ],
                         df.adj.lung[4, ],
                         df.adj.gi[5, ],
                         c(0, 0, 0, 0, 0, 0)); df.adj.outflows

df.adj.outflows$tissue_from <- c("Nose", "Throat", "Trachea", "Lung",
                                 "Upper GI", "Lower GI")

# Convert to long format for standardizing & plotting
df.adj.outflows.long <- subset(df.adj.outflows, select = -tissue_inoculated) %>%
  pivot_longer(cols = !tissue_from, values_to = "value", names_to = "tissue_to") %>%
  group_by(tissue_from) %>%
  arrange(tissue_from, value) %>%
  mutate(order = row_number(),
         rel_value = value / max(value, na.rm = TRUE))

df.adj.outflows.long$tissue_to <- gsub("\\.", " ", df.adj.outflows.long$tissue_to)


# Make some necessary adjustements for plotting purposes
df.adj.outflows.long$order[df.adj.outflows.long$tissue_from == "Lower GI"] <- 2
df.adj.outflows.long$order[df.adj.outflows.long$tissue_from == df.adj.outflows.long$tissue_to] <- 
  df.adj.outflows.long$tissue_from[df.adj.outflows.long$tissue_from == df.adj.outflows.long$tissue_to]

## Plot -------------------------------------------------------------------------

fig.outflows <- 
  ggplot(df.adj.outflows.long) +
  geom_tile(aes(x = tissue_to, y = tissue_from), 
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, rel_value >= 0.9),
            aes(x = tissue_to, y = tissue_from), 
            fill = "grey33", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, rel_value == 1),
            aes(x = tissue_to, y = tissue_from), 
            fill = "black", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_from == "Lower GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Nose" & tissue_from == "Nose"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Nose, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Throat" & tissue_from == "Throat"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Throat, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Trachea" & tissue_from == "Trachea"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Trachea, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Lung" & tissue_from == "Lung"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Lung, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Upper GI" & tissue_from == "Upper GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Upper.GI, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.outflows.long, tissue_to == "Lower GI" & tissue_from == "Lower GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Lower.GI, color = "black", linewidth = 0.4) +
  geom_text(data = subset(df.adj.outflows.long, tissue_from != tissue_to),
            aes(x = tissue_to, y = tissue_from, label = round(rel_value, 2))) +
  geom_text(data = subset(df.adj.outflows.long, order == 6 | rel_value >= 0.9),
            aes(x = tissue_to, y = tissue_from, label = round(rel_value, 2)), color = "white") +
  
  scale_fill_gradient(low = "white", high = "grey11") +
  scale_color_gradient(low = "black", high = "white") +
  scale_x_discrete(limits = c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI"),
                 expand = c(0, 0)) +
  scale_y_discrete(limits = rev(c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI")),
                   expand = c(0, 0)) +
  facet_wrap(.~ "Outflows") +
  labs(x = "Tissue To", y = "Tissue From") + 
  guides(fill = "none") +
  theme(panel.background = element_rect(fill ="grey", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()); fig.outflows


# Figure 5B  -------------------------------------------------------------------

## Calculate inflows -----------------------------------------------------------

# Get the columns from all the adjacency matrices, except for the tissue being
#    studied, and compute the average inflow across routes
nose_inflow <- rowMeans(cbind(df.adj.throat[, 1], df.adj.trachea[, 1], df.adj.lung[, 1], df.adj.gi[, 1]))
throat_inflow <- rowMeans(cbind(df.adj.nose[, 2], df.adj.trachea[, 2], df.adj.lung[, 2], df.adj.gi[, 2]))
trachea_inflow <- rowMeans(cbind(df.adj.nose[, 3], df.adj.throat[, 3], df.adj.lung[, 3], df.adj.gi[, 3]))
lung_inflow <- rowMeans(cbind(df.adj.nose[, 4], df.adj.throat[, 4], df.adj.trachea[, 4], df.adj.gi[, 4]))
upgi_inflow <- rowMeans(cbind(df.adj.nose[, 5], df.adj.throat[, 5], df.adj.trachea[, 5], df.adj.lung[, 5]))
logi_inflow <- rowMeans(cbind(df.adj.nose[, 6], df.adj.throat[, 6], df.adj.trachea[, 6], df.adj.lung[, 6], df.adj.gi[, 6]))

# Combine the averages into one dataframe
df.adj.inflows <- as.data.frame(
  cbind(nose_inflow, throat_inflow, trachea_inflow,
        lung_inflow, upgi_inflow, logi_inflow)); df.adj.inflows

# Set the column and row names
colnames(df.adj.inflows) <- c("Nose", "Throat", "Trachea", "Lung",
                              "Upper GI", "Lower GI")
df.adj.inflows$tissue_from <- c("Nose", "Throat", "Trachea", "Lung",
                                "Upper GI", "Lower GI")

# Convert to long format, and find the correct ordering & relative effects
df.adj.inflows.long <- df.adj.inflows %>%
  pivot_longer(cols = !tissue_from, values_to = "value", names_to = "tissue_to") %>%
  group_by(tissue_to) %>%
  arrange(tissue_to, value) %>%
  mutate(order = row_number(),
         rel_value = value / max(value, na.rm = TRUE))


# Make some necessary adjustments for plotting purposes
df.adj.inflows.long$order[df.adj.inflows.long$tissue_from == df.adj.inflows.long$tissue_to] <- 
  df.adj.inflows.long$tissue_from[df.adj.inflows.long$tissue_from == df.adj.inflows.long$tissue_to]


## Plot ------------------------------------------------------------------------

fig.inflows <- 
  ggplot(df.adj.inflows.long) +
  geom_tile(aes(x = tissue_to, y = tissue_from), 
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, rel_value >= 0.9),
            aes(x = tissue_to, y = tissue_from), 
            fill = "grey33", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, rel_value == 1),
            aes(x = tissue_to, y = tissue_from), 
            fill = "black", color = "black", linewidth = 0.4) +
  
  #geom_tile(aes(x = tissue_to, y = tissue_from, fill = rel_value),
  #          color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_from == "Lower GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Nose" & tissue_from == "Nose"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Nose, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Throat" & tissue_from == "Throat"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Throat, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Trachea" & tissue_from == "Trachea"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Trachea, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Lung" & tissue_from == "Lung"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Lung, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Upper GI" & tissue_from == "Upper GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Upper.GI, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Lower GI" & tissue_from == "Lower GI"),
            aes(x = tissue_to, y = tissue_from),
            fill = assign_colors()$Lower.GI, color = "black", linewidth = 0.4) +
  geom_text(data = subset(df.adj.inflows.long, tissue_from != tissue_to),
            aes(x = tissue_to, y = tissue_from, label = round(rel_value, 2))) +
  geom_text(data = subset(df.adj.inflows.long, order == 6 | rel_value >= 0.9),
            aes(x = tissue_to, y = tissue_from, label = round(rel_value, 2)), color = "white") +
  scale_fill_gradient(low = "white", high = "grey11") +
  scale_x_discrete(limits = c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI"),
                   expand = c(0, 0)) +
  scale_y_discrete(limits = rev(c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI")),
                   expand = c(0, 0)) +
  facet_wrap(.~ "Inflows") +
  labs(y = "Tissue From", x = "Tissue To") + 
  guides(fill = "none") +
  theme(panel.background = element_rect(fill ="grey", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()); fig.inflows


# Figure 5C --------------------------------------------------------------------

# Load the diagram constructed in Powerpoint & saved as a png
img <- readPNG("./outputs/figures/fig4-connectivity-schematic.png")
img_grob <- rasterGrob(img, interpolate = TRUE)

fig.diagram <- ggplot() +
  annotation_custom(img_grob,
                    xmin = -Inf, xmax = Inf,
                    ymin = -Inf, ymax = Inf) +
  theme_void()


# Figure 5D --------------------------------------------------------------------

## Compute the probability of positivity ---------------------------------------

df.percent <- df.all %>%
  group_by(tissue_name, tissue_inoculated, ever_positive) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(median_percent =  n / sum(n))

df.percent <- subset(df.percent, ever_positive == 1)

df.prob <- df.all %>%
  group_by(tissue_name, tissue_inoculated) %>%
  summarise(median_percent = median(percent, na.rm = TRUE))


## Plot ------------------------------------------------------------------------

fig.percent <- ggplot(df.prob) +
  geom_tile(data = subset(df.prob, median_percent >= 0.9), 
            aes(y = tissue_inoculated, x = tissue_name),
                fill = "black", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.prob, median_percent < 0.9), 
            aes(y = tissue_inoculated, x = tissue_name),
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Lower GI"),
            aes(y = tissue_to, x = tissue_from),
            fill = "white", color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Nose" & tissue_from == "Nose"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Nose, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Throat" & tissue_from == "Throat"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Throat, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Trachea" & tissue_from == "Trachea"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Trachea, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Lung" & tissue_from == "Lung"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Lung, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Upper GI" & tissue_from == "Upper GI"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Upper.GI, color = "black", linewidth = 0.4) +
  geom_tile(data = subset(df.adj.inflows.long, tissue_to == "Lower GI" & tissue_from == "Lower GI"),
            aes(y = tissue_to, x = tissue_from),
            fill = assign_colors()$Lower.GI, color = "black", linewidth = 0.4) +
  geom_text(aes(y = tissue_inoculated, x = tissue_name, label = round(100 * median_percent, 0)),
            color = "black") +
  geom_text(data = subset(df.prob, median_percent >= 0.9),
            aes(y = tissue_inoculated, x = tissue_name, label = round(100 * median_percent, 0)),
            color = "white") +
  scale_x_discrete(limits = c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI"),
                   expand = c(0, 0)) +
  scale_y_discrete(limits = rev(c("Nose", "Throat", "Trachea", "Lung", "Upper GI", "Lower GI")),
                   expand = c(0, 0)) +
  scale_fill_gradient(low = "white", high = "grey11") +
  facet_wrap(.~ "Probability of Positivity") +
  labs(x = "Tissue Sampled", y = "Tissue Inoculated") + 
  guides(fill = "none") +
  theme(panel.background = element_rect(fill ="grey", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank()); fig.percent



# Combine ----------------------------------------------------------------------

fig5 <- fig.outflows + labs(tag = "A") + theme(legend.position = "none") + 
  fig.inflows + labs(tag = "B") + #theme(legend.position = "bottom") + 
  fig.diagram + labs(tag = "C") +
  fig.percent + labs(tag = "D") +
  plot_layout(nrow = 1, widths = c(1, 1, 1.3, 1)); fig5



# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/fig5-tissue-connectivity-structure.pdf",
       plot = fig5,
       width = 12.5, 
       height = 3.3,
       dpi = 600)


