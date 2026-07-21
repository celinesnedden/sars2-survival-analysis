# This file: - Shows the number of datapoints available for all possible pairs of cofactors


# Prep -------------------------------------------------------------------------

# Load the event time database
dat <- read.csv("./data/df-event-times.csv")

# Assign all cofactor names
dat <- assign_all_names(dat)

# Group doses by order of magnitude
dat$dose_grouped <- floor(log10(dat$inoc_dose_total_pfu))


# Column 1: Route --------------------------------------------------------------

# Route & Route
dat.count.route.route <- dat %>%
  count(route_name)

fig.route.route <- ggplot(dat.count.route.route) +
  geom_tile(aes(x = route_name, y = factor(route_name, levels = rev(levels(dat$route_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = route_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Exposure Route") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.route


# Route & Dose
dat.count.route.dose <- dat %>%
  count(route_name, dose_grouped) %>%
  complete(route_name, dose_grouped, fill = list(n = 0))

fig.route.dose <- ggplot(dat.count.route.dose) +
  geom_tile(aes(x = route_name,  y = factor(as.character(dose_grouped), levels = rev(seq(1, 7))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name,  y = factor(as.character(dose_grouped), levels = rev(seq(1, 7))), 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Exposure Dose (log10 pfu)") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.dose

# Route & species
dat.count.route.species <- dat %>%
  count(route_name, species_name) %>%
  complete(route_name, species_name, fill = list(n = 0))

fig.route.species <- ggplot(dat.count.route.species) +
  geom_tile(aes(x = route_name, y = factor(species_name, levels = rev(levels(dat$species_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = factor(species_name, levels = rev(levels(dat$species_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Species") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.species


# Route & age
dat.count.route.age <- dat %>%
  count(route_name, age_name) %>%
  complete(route_name, age_name, fill = list(n = 0))

fig.route.age <- ggplot(dat.count.route.age) +
  geom_tile(aes(x = route_name, y = factor(age_name, levels = rev(levels(dat$age_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = age_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Age") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.age


# Route & sex
dat.count.route.sex <- dat %>%
  count(route_name, sex_name) %>%
  complete(route_name, sex_name, fill = list(n = 0))

fig.route.sex <- ggplot(dat.count.route.sex) +
  geom_tile(aes(x = route_name, y = factor(sex_name, levels = rev(levels(dat$sex_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = sex_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.sex

# Route & location
dat$tissue_name <- droplevels(dat$tissue_name)

dat.count.route.location <- dat %>%
  count(route_name, tissue_name) %>%
  complete(route_name, tissue_name, fill = list(n = 0))

fig.route.location <- ggplot(dat.count.route.location) +
  geom_tile(aes(x = route_name, y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.location


# Route & assay
dat.count.route.assay <- dat %>%
  count(route_name, assay_type) %>%
  complete(route_name, assay_type, fill = list(n = 0))

fig.route.assay <- ggplot(dat.count.route.assay) +
  geom_tile(aes(x = route_name, y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = route_name, y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Exposure Route") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.route.assay

remove_xaxis <- theme(legend.position = "none", 
                      axis.text.x = element_blank(), 
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank())

fig.route <- 
  (fig.route.route + remove_xaxis) / 
  (fig.route.dose + remove_xaxis) /
  (fig.route.species + remove_xaxis) / 
  (fig.route.age + remove_xaxis) / 
  (fig.route.sex + remove_xaxis) /
  (fig.route.location + remove_xaxis) /
  fig.route.assay + theme(legend.position = "none") +
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.route


# Column 2: Dose ---------------------------------------------------------------

# Dose 
dat.count.dose.dose <- dat %>%
  count(dose_grouped) 

fig.dose.dose <- ggplot(dat.count.dose.dose) +
  geom_tile(aes(x = as.character(dose_grouped), 
                y = factor(as.character(dose_grouped), levels = rev(seq(1, 7))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = dose_grouped, y = factor(as.character(dose_grouped), levels = rev(seq(1, 7))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Exposure Dose (log10 pfu)") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.dose

# Dose &  species
dat.count.dose.species <- dat %>%
  count(dose_grouped, species_name) %>%
  complete(dose_grouped, species_name, fill = list(n = 0))

fig.dose.species <- ggplot(dat.count.dose.species) +
  geom_tile(aes(x = as.character(dose_grouped), 
                y = factor(species_name, levels = rev(levels(dat$species_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = as.character(dose_grouped), 
                y = species_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  scale_x_discrete() +
  labs(y = "Species") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.species


# Dose &  age
dat.count.dose.age <- dat %>%
  count(dose_grouped, age_name) %>%
  complete(dose_grouped, age_name, fill = list(n = 0))

fig.dose.age <- ggplot(dat.count.dose.age) +
  geom_tile(aes(x = as.character(dose_grouped), 
                y = factor(age_name, levels = rev(levels(dat$age_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = as.character(dose_grouped), 
                y = age_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Age") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.age


# Dose &  sex
dat.count.dose.sex <- dat %>%
  count(dose_grouped, sex_name) %>%
  complete(dose_grouped, sex_name, fill = list(n = 0))

fig.dose.sex <- ggplot(dat.count.dose.sex) +
  geom_tile(aes(x = as.character(dose_grouped), 
                y = factor(sex_name, levels = rev(levels(dat$sex_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = dose_grouped, y = sex_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.sex

# Dose &  location
dat$tissue_name <- droplevels(dat$tissue_name)

dat.count.dose.location <- dat %>%
  count(dose_grouped, tissue_name) %>%
  complete(dose_grouped, tissue_name, fill = list(n = 0))

fig.dose.location <- ggplot(dat.count.dose.location) +
  geom_tile(aes(x = as.character(dose_grouped), 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = as.character(dose_grouped),
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.location


# Dose &  assay
dat.count.dose.assay <- dat %>%
  count(dose_grouped, assay_type) %>%
  complete(dose_grouped, assay_type, fill = list(n = 0))

fig.dose.assay <- ggplot(dat.count.dose.assay) +
  geom_tile(aes(x = as.character(dose_grouped), y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = as.character(dose_grouped), y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Exposure Dose (log10 pfu)") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.dose.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())

fig.dose <- plot_spacer() / 
  (fig.dose.dose + remove_axes) /
  (fig.dose.species + remove_axes) / 
  (fig.dose.age + remove_axes) / 
  (fig.dose.sex + remove_axes) /
  (fig.dose.location + remove_axes) /
  fig.dose.assay + theme(legend.position = "none", 
                         axis.text.y = element_blank(), 
                         axis.title.y = element_blank(),
                         axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.dose


# Column 3: Species ----------------------------------------------------------

# Dose & species
dat.count.species.species <- dat %>%
  count(species_name)

fig.species.species <- ggplot(dat.count.species.species) +
  geom_tile(aes(x = species_name, 
                y = factor(species_name, levels = rev(levels(dat$species_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = species_name, 
                y = factor(species_name, levels = rev(levels(dat$species_name))), 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  scale_x_discrete() +
  labs(y = "Species") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.species.species


# Dose &  age
dat.count.species.age <- dat %>%
  count(species_name, age_name) %>%
  complete(species_name, age_name, fill = list(n = 0))

fig.species.age <- ggplot(dat.count.species.age) +
  geom_tile(aes(x = species_name, 
                y = factor(age_name, levels = rev(levels(dat$age_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = species_name, 
                y = age_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Age") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.species.age


# Dose &  sex
dat.count.species.sex <- dat %>%
  count(species_name, sex_name) %>%
  complete(species_name, sex_name, fill = list(n = 0))

fig.species.sex <- ggplot(dat.count.species.sex) +
  geom_tile(aes(x = species_name, 
                y = factor(sex_name, levels = rev(levels(dat$sex_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = species_name, y = sex_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.species.sex

# Dose &  location

dat.count.species.location <- dat %>%
  count(species_name, tissue_name) %>%
  complete(species_name, tissue_name, fill = list(n = 0))

fig.species.location <- ggplot(dat.count.species.location) +
  geom_tile(aes(x = species_name, 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = species_name,
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.species.location


# Dose &  assay

dat.count.species.assay <- dat %>%
  count(species_name, assay_type) %>%
  complete(species_name, assay_type, fill = list(n = 0))

fig.species.assay <- ggplot(dat.count.species.assay) +
  geom_tile(aes(x = species_name, y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = species_name, y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Species") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.species.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())

fig.species <- plot_spacer() / 
  plot_spacer() /
  (fig.species.species + remove_axes) / 
  (fig.species.age + remove_axes) / 
  (fig.species.sex + remove_axes) /
  (fig.species.location + remove_axes) /
  fig.species.assay + theme(legend.position = "none", 
                            axis.text.y = element_blank(), 
                            axis.title.y = element_blank(),
                            axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.species


# Column 4: Age ---------------------------------------------------------------

# Dose &  age
dat.count.age.age <- dat %>%
  count(age_name)

fig.age.age <- ggplot(dat.count.age.age) +
  geom_tile(aes(x = age_name, 
                y = factor(age_name, levels = rev(levels(dat$age_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = age_name, 
                y = factor(age_name, levels = rev(levels(dat$age_name))), 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  scale_x_discrete() +
  labs(y = "age") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.age.age


# Dose &  sex
dat.count.age.sex <- dat %>%
  count(age_name, sex_name) %>%
  complete(age_name, sex_name, fill = list(n = 0))

fig.age.sex <- ggplot(dat.count.age.sex) +
  geom_tile(aes(x = age_name, 
                y = factor(sex_name, levels = rev(levels(dat$sex_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = age_name, y = sex_name, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.age.sex

# Dose &  location

dat.count.age.location <- dat %>%
  count(age_name, tissue_name) %>%
  complete(age_name, tissue_name, fill = list(n = 0))

fig.age.location <- ggplot(dat.count.age.location) +
  geom_tile(aes(x = age_name, 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = age_name,
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.age.location


# Dose &  assay

dat.count.age.assay <- dat %>%
  count(age_name, assay_type) %>%
  complete(age_name, assay_type, fill = list(n = 0))

fig.age.assay <- ggplot(dat.count.age.assay) +
  geom_tile(aes(x = age_name, y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = age_name, y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Age") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.age.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())

fig.age <- plot_spacer() / 
  plot_spacer() /
  plot_spacer() / 
  (fig.age.age + remove_axes) / 
  (fig.age.sex + remove_axes) /
  (fig.age.location + remove_axes) /
  fig.age.assay + theme(legend.position = "none", 
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(),
                        axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.age


# Column 5: Sex ----------------------------------------------------------------

# Dose &  sex
dat.count.sex.sex <- dat %>%
  count(sex_name)

fig.sex.sex <- ggplot(dat.count.sex.sex) +
  geom_tile(aes(x = sex_name, 
                y = factor(sex_name, levels = rev(levels(dat$sex_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = sex_name, 
                y = factor(sex_name, levels = rev(levels(dat$sex_name))), 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  scale_x_discrete() +
  labs(y = "sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.sex.sex

# Dose &  location

dat.count.sex.location <- dat %>%
  count(sex_name, tissue_name) %>%
  complete(sex_name, tissue_name, fill = list(n = 0))

fig.sex.location <- ggplot(dat.count.sex.location) +
  geom_tile(aes(x = sex_name, 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = sex_name,
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.sex.location


# Dose &  assay

dat.count.sex.assay <- dat %>%
  count(sex_name, assay_type) %>%
  complete(sex_name, assay_type, fill = list(n = 0))

fig.sex.assay <- ggplot(dat.count.sex.assay) +
  geom_tile(aes(x = sex_name, y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = sex_name, y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Sex") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.sex.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())

fig.sex <- plot_spacer() / 
  plot_spacer() /
  plot_spacer() / 
  plot_spacer() / 
  (fig.sex.sex + remove_axes) /
  (fig.sex.location + remove_axes) /
  fig.sex.assay + theme(legend.position = "none", 
                        axis.text.y = element_blank(), 
                        axis.title.y = element_blank(),
                        axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.sex


# Column 6: Tissue ----------------------------------------------------------

# Dose &  location
dat.count.location.location <- dat %>%
  count(tissue_name)

fig.location.location <- ggplot(dat.count.location.location) +
  geom_tile(aes(x = tissue_name, 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = tissue_name, 
                y = factor(tissue_name, levels = rev(levels(dat$tissue_name))), 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  scale_x_discrete() +
  labs(y = "location") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.location.location


# Dose &  assay

dat.count.location.assay <- dat %>%
  count(tissue_name, assay_type) %>%
  complete(tissue_name, assay_type, fill = list(n = 0))

fig.location.assay <- ggplot(dat.count.location.assay) +
  geom_tile(aes(x = tissue_name, y = assay_type, fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = tissue_name, y = assay_type, label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36)) +
  labs(y = "Assay", x = "Tissue") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.9, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "white"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.location.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())

fig.location <- plot_spacer() / 
  plot_spacer() /
  plot_spacer() / 
  plot_spacer() / 
  plot_spacer() / 
  (fig.location.location + remove_axes) /
  fig.location.assay + theme(legend.position = "none", 
                             axis.text.y = element_blank(), 
                             axis.title.y = element_blank(),
                             axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.location


# Column 7: Assay --------------------------------------------------------------

# Dose &  assay
dat.count.assay.assay <- dat %>%
  count(assay_type)

fig.assay.assay <- ggplot(dat.count.assay.assay) +
  geom_tile(aes(x = factor(assay_type, levels = c("Unknown", "totRNA", "sgRNA", "gRNA", "culture")), 
                y = assay_type, 
                fill = log10(n)), 
            linewidth = 0.4, color = "black") +
  geom_text(aes(x = factor(assay_type, levels = c("Unknown", "totRNA", "sgRNA", "gRNA", "culture")), 
                y = assay_type, 
                label = n)) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36),
                      breaks = seq(0, 3, 1), labels = c(1, 10, 100, 1000),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_discrete() +
  labs(y = "assay", x = "Assay", fill = "Sample\nSize") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, hjust = 0.5, vjust = 1),
        legend.key.size = unit(1, "line"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.key = element_rect(color = "black"),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_rect(fill = "transparent", color = "black"),
        strip.text = element_text(size = 11),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()); fig.assay.assay

remove_axes <- theme(legend.position = "none", 
                     axis.text.x = element_blank(), 
                     axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y = element_blank(), 
                     axis.title.y = element_blank(),
                     axis.ticks.y = element_blank())


fig.legend <- get_legend(fig.assay.assay)


fig.assay <- plot_spacer() / 
  plot_spacer() /
  plot_spacer() / 
  plot_spacer() / 
  plot_spacer() / 
  fig.legend /
  fig.assay.assay + theme(legend.position = "none", 
                          axis.text.y = element_blank(), 
                          axis.title.y = element_blank(),
                          axis.ticks.y = element_blank()) + 
  plot_layout(heights = c(5, 7, 3, 4, 3, 6, 5)); fig.assay



# Combine Columns --------------------------------------------------------------

fig.S3 <- (fig.route | fig.dose | fig.species | fig.age | fig.sex | fig.location | fig.assay) + 
  plot_layout(widths = c(5, 7, 3, 4, 3, 6, 5)); fig.S3



# Table for Metrics --------------------------------------------------------------

# Data passed to stan for plotting sample sizes
dat.stan <- readRDS("./data/df-passed-to-Stan.RDS")

# Categorize assays by PCR vs. Culture
dat.stan$assay_type[dat.stan$assay %in% c(1:3, -9999)] <- "Any PCR assay"
dat.stan$assay_type[dat.stan$assay %in% 4] <- "Culture"

# Set up empty dataframe 
df <- data.frame(metric = NA, route_idx = NA, assay = NA, tissue_idx = NA, n = NA)

# Fill with sample size values
for (tissue.ii in 1:6) {
  for (route.ii in 1:5){
    for (assay.ii in c("Any PCR assay", "Culture")) {
      #cat("tissue", tissue.ii, " route ", route.ii, " assay", assay.ii)
      
      num_first <- length(dat.stan$ever_positive[dat.stan$route == route.ii &
                                                   dat.stan$tissue_location == tissue.ii])
      num_peak <- sum(dat.stan$has_peak[dat.stan$route == route.ii &
                                          dat.stan$assay_type == assay.ii &
                                          dat.stan$tissue_location == tissue.ii])
      num_titer <- sum(dat.stan$has_titer[dat.stan$route == route.ii &
                                            dat.stan$assay_type == assay.ii &
                                            dat.stan$tissue_location == tissue.ii])
      num_last <- sum(dat.stan$has_last_positive[dat.stan$route == route.ii &
                                                   dat.stan$assay_type == assay.ii &
                                                   dat.stan$tissue_location == tissue.ii])
      
      print(num_first); print(num_peak); print(num_titer); print(num_last)
      
      df.add <- data.frame(metric = c("Percent positive /\nTime to detectability", 
                                      "Time to peak titer", "Peak titer", 
                                      "Time to\nundetectability"),
                           route_idx = rep(route.ii, 4),
                           assay = rep(assay.ii, 4),
                           tissue_idx = rep(tissue.ii, 4),
                           n = c(num_first, num_peak, num_titer, num_last))
      
      df <- rbind(df, df.add)
    }
  }
}

df <- assign_route_names(df)

df$tissue_name[df$tissue_idx == 1] <- "Nose"
df$tissue_name[df$tissue_idx == 2] <- "Throat"
df$tissue_name[df$tissue_idx == 3] <- "Trachea"
df$tissue_name[df$tissue_idx == 4] <- "Lung"
df$tissue_name[df$tissue_idx == 5] <- "Upper\nGI"
df$tissue_name[df$tissue_idx == 6] <- "Lower\nGI"
df$tissue_name <- factor(df$tissue_name, levels = c("Nose",
                                                    "Throat",
                                                    "Trachea",
                                                    "Lung",
                                                    "Upper\nGI",
                                                    "Lower\nGI"))


df <- df[-1, ]

tbl <- ggplot(df) +
  geom_tile(aes(x = route_name, y = assay, fill = log10(n)), color = "black") +
  geom_text(aes(x = route_name, y = assay, label = n),
            size = 3) +
  geom_tile(data = subset(df, tissue_name == "Upper GI" & metric != "Percent positive /\nTime to detectability"),
            aes(x = route_name, y = assay), fill = "grey") +
  facet_grid(tissue_name ~ factor(metric, levels = c("Percent positive /\nTime to detectability", 
                                                     "Time to peak titer", "Peak titer", 
                                                     "Time to\nundetectability"))) +
  scale_fill_gradient(low = "white", high = "#00B0F6", limits = c(0, 3.36),
                      breaks = seq(0, 3, 1), labels = c(1, 10, 100, 1000),
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "Exposure Route", y = "", tag = "b") +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35,  hjust = 1, vjust = 1),
        text = element_text(size = 10),
        legend.position = "none",
        strip.background = element_rect(fill = "white", color = "transparent"),
        strip.text = element_text(size = 9, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA_character_),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA, 
                                    colour = "black",
                                    size=0.5)); tbl



# Combine full figure ----------------------------------------------------------

fig.comb <- (wrap_elements(fig.S3) + labs(tag = "a"))  +
  inset_element(
    tbl,
    left = 0.45,
    bottom = 0.63,
    right = 1,
    top = 1
  ); test


# Save -------------------------------------------------------------------------

ggsave("./outputs/figures/figS3-data-by-cofactor-pairs.png",
       plot = fig.comb,
       width = 14, 
       height = 11.8,
       dpi = 300)

