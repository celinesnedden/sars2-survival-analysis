# Exposure route, not dose, is the primary driver of infection patterns for a respiratory virus

Celine Snedden (1), Dylan Morris (1), Thomas Friedrich (2, 3), James Lloyd-Smith (1, 4)

1. Department of Ecology and Evolutionary Biology, University of California, Los Angeles, Los Angeles, CA, USA
2. Department of Pathobiological Sciences, University of Wisconsin-Madison; Madison, Wisconsin, USA.
3. Wisconsin National Primate Research Center; Madison, Wisconsin, USA.
4. Department of Computational Medicine, University of California, Los Angeles, Los Angeles, CA, USA 


## Repository information

This repository provides the data and code to reproduce the analyses presented in the associated paper (Snedden et al., 2026). Instructions to recreate all analyses are summarized below, as is a description of the database.

Please note that the data from one study (Johnston et al., 2021; https://doi.org/10.1371/journal.pone.0246366) used in our analysis is not included in any of the data files due to a data sharing agreement. Interested readers can submit their own requests following the instructions in the source paper. Because of this, reproductions of our sample sizes calculations and analyses may yield different results due to reduced samples sizes. 

## Article abstract

Infectious disease severity, shedding, and within-host kinetics can vary substantially with exposure dose, inoculation route, and host factors. However, the relative contributions of these variables to infection heterogeneity remain poorly quantified because small sample sizes in controlled in vivo experiments limit their statistical power and scientific scope. Here, by compiling and jointly analyzing the largest published database of non-human primate challenge experiments (107 studies; 721 animals), we show that exposure route drives SARS-CoV-2 infection kinetics more strongly than dose, age, sex, or species. Aerosol inoculation yields kinetics distinct from all other respiratory exposure routes. Route and dose jointly determine which tissues become infected, and 50% infectious doses vary greatly depending on the route of exposure and on the tissue sampled. Our findings underscore the central role of exposure route in pathogenesis and transmission, and highlight the untapped potential for quantitative, cross-study syntheses to extract additional insights from costly animal experiments.


## Citation information

If you use any of the code or data associated with this manuscript, please cite our work. (_Citation information to come_). 

Note that this work was posted as a preprint before publication (https://doi.org/10.64898/2026.03.03.709384).


## File structure and naming conventions

All content is organized into the following folders:

- `code`: contains all code, including generating event times, the Stan model, model fitting, generating model predictions, and figure / table generation. 
  - primary analyses are ordered by increasing numbers (e.g., `01-analysis-prep.R`, `02-extract-event-times.R`) representing the sequence in which they were run for the investigations presented in the paper.
  - unless otherwise specified, each file should run without needing to run any other files prior 
  - files that generate certain figures or tables are preceded with `fig` or `tbl` and the corresponding number from the paper.

- `data`: contains data as `.csv` files, including the database (`database-clean.csv`) and the extracted event times (`df-event-times.csv`). 

- `outputs`: contains all output files, including model fits, model predictions, and figures / tables.
  - just as for code files, each output is marked as a figure (`fig`) or a table (`tbl`) preceding the corresponding number. In many cases, the code generating the figure/table can be found with the corresponding name in the code file. Model predictions are generated in `04-generate-predictions.R`
  - `outputs/fits`: contains the final model fit.


## Installing software dependencies

To run our analyses, you will need to install certain software, described below. 

- We use the statistical programming language R, for which you can find installation instructions at: https://cran.r-project.org/doc/manuals/r-release/R-admin.html. 

- We use CmdStanR for Bayesian model compilation. We recommend that you: 

  1. review and follow the installation instructions at https://mc-stan.org/cmdstanr/articles/cmdstanr.html
  
  2. run the example models to confirm performance is as expected before proceeding. 

- We also use various other R packages. The 01-analysis-prep.R file will prompt installation of these packages. 


## Database description

Each column is explained below in the order they appear in `database-clean.csv`. Please note the following general conventions:
- `_idx` variables (`assay_idx`, `route_idx`, `age_idx`, `sex_idx`, `sp_idx`, `location_idx`) are integer encodings used for statistical modeling, with -9999 representing an unknown where applicable.
- `prob_male`, `prob_female`, `prob_geriatric` are probabilistic assignments when only group-level information (rather than individual metadata) was available. These are used for statistical modeling. 
- `dose_` variables are log10 doses delivered to specific anatomical compartments, based on our assumptions on dose distribution and the reported inoculation routes (see Materials & Methods).
- `_rep` columns (e.g., `unit_rep`, `sample_rep`) retains the original study terminology, whereas `_grp` and `_subgrp` are standardized versions at varying levels of specificity. 


| Column | Description |
|---|---| 
| article | The citation identifier for the study. |
| indiv | The unique animal identifier. To prevent confusion arising due to repeated ID names across studies (e.g., 'RM1' is commonly used as an identifier), we have prepended all source ID names with the initials of the study's first author. For example, individual 'RM1' from (Vincent) Munster et al. 2020 is referred to as `VM_RM1` in this database. In instances where no identifier was assigned or used, we define new identifiers from the first author's initials and the species of primate ('AGM' for african green monkeys, 'CM' for cynomolgus macaques, and 'RM' for rhesus macaques). For example, the two rhesus macaques included in Wang et al. 2020 are referred to as `HW_RM1` and `HW_RM2`. When no identifier is assigned by the authors, it is usually not possible to determine which samples were obtained from the same individual, and this uncertainty is incorporated into our naming system. If it is obvious within a single data type (e.g., nasal swab) that specimens across sampling days were derived from the same individual but it is not clear which trajectories come from the same individual across sample types (e.g., nasal swabs and respiratory tissues), then our identifier is appended with a single `*`. This happens most often when lines connect data points within a figure panel, but no color/symbol scheme is used to distinguish correlated lines across panels. Alternatively, if it is unclear which specimens across sampling times or sample types were derived from the same individual, then our identifier is appended with `**`. This occurs when the data is presented as scatter plots without lines connecting sampling times and without a color/symbol scheme. Any individual tagged as "No ID" in source ID means ID names / references were not available and it was not possible to correlate samples across tissues or timepoints. |
| day_post_infection | Days after viral inoculation when the sample was collected. |
| pos_value | Indicator of whether virus was detected (1 = detected, 0 = not detected). |
| value | Reported measurement value. May contain numeric values or non-numeric entries (e.g., `< LOD`, `< LOQ`, approximate values, Ct values with uncertainty). When numeric (and not flagged as Ct values), the values are in log10 units. |
| sample_type | Whether the sample was collected invasively (tissue samples at necropsy) or non-invasively (e.g, swabs). |
| location_grp | Broad anatomical sampling location. |
| location_subgrp | Intermediate anatomical location category. |
| sample_rep | Original sample name reported by the study. |
| sample_grp | Standardized sample type. |
| sample_subgrp | More specific sample subtype or anatomical location. |
| organ_system | Organ system associated with the sample. |
| method_of_quant | Method used to quantify virus (e.g., RT-qPCR or culture). |
| rna_type | Viral target measured (genomic RNA, subgenomic RNA, total RNA, culture, etc.). |
| pcr_target_gene | Viral gene targeted by PCR assay. |
| pcr_protocol_type | PCR assay or primer/protocol used. |
| assay_idx | Numeric identifier for assay type. |
| inoc_route_grp | Standardized inoculation route category. |
| route_idx | Numeric identifier for inoculation route. |
| dose_nose | Log10 inoculum delivered to the nasal compartment, according to our assumptions (see Materials & Methods). Note that 0 indicates no virus was administered in that compartment. |
| dose_throat | Log10 inoculum delivered to the throat/oropharynx, according to our assumptions (see Materials & Methods). Note that 0 indicates no virus was administered in that compartment. |
| dose_trachea | Log10 inoculum delivered to the trachea, according to our assumptions (see Materials & Methods). Note that 0 indicates no virus was administered in that compartment. |
| dose_lung | Log10 inoculum delivered directly to the lungs, according to our assumptions (see Materials & Methods). Note that 0 indicates no virus was administered in that compartment. |
| dose_gi | Log10 inoculum delivered to the gastrointestinal tract, according to our assumptions (see Materials & Methods). Note that 0 indicates no virus was administered in that compartment. |
| age_class_grp | Standardized age category. |
| age_idx | Numeric identifier for age category. |
| prob_geriatric | Probability that an individual belongs to the geriatric age class. This depends on the reported age distribution in the study. |
| sex | Reported sex of the animal(s). |
| sex_idx | Numeric encoding of sex. |
| prob_male | Probability that an individual is male. This depends on the reported sex distribution in the study. |
| prob_female | Probability that an individual is female. This depends on the reported sex distribution in the study. |
| animal_species | Non-human primate species. |
| sp_idx | Numeric identifier for species. |
| unit_rep | Original measurement units reported by the study. |
| unit_subgrp | Measurement units grouped into similar types. |
| data_source | Location within the publication where data were obtained. |
| llod | Lower limit of detection (reported on the measurement scale). |
| lod_source | Source of the reported LOD information. |
| lloq | Lower limit of quantification (reported on the measurement scale). |
| loq_source | Source of the reported LOQ information. |
| viral_strain_rep | Original viral strain name reported by the study. |
| viral_strain_grp | Standardized viral strain designation. |
| treatment | Experimental treatment or control group reported. |
| animal_subspecies | Animal subspecies or origin, when reported. |
| age_class_rep | Original age description reported by the study. |
| age_years | Reported animal age in years or age range. |
| inoc_route_rep | Original inoculation route reported by the study. |
| inoc_route_subgrp | Standardized inoculation route description. |
| inoc_dose_units_rep | Units used in the study to report the inoculum dose. |
| inoc_dose_TCID_per_mL | Inoculum concentration in TCID50/mL. |
| inoc_dose_pfu_per_mL | Inoculum concentration in pfu/mL. |
| inoc_dose_total_mL | Total inoculum volume administered. |
| inoc_dose_total_TCID | Total inoculum dose in TCID50. |
| inoc_dose_total_pfu | Total inoculum dose in pfu. |
| inoc_dose_IT_mL | Volume administered intratracheally. |
| inoc_dose_IN_mL | Volume administered intranasally. |
| inoc_dose_OR_mL | Volume administered orally. |
| inoc_dose_OC_mL | Volume administered ocularly/conjunctivally. |
| inoc_dose_IG_mL | Volume administered intragastrically. |
| inoc_dose_IV_mL | Volume administered intravenously. |
| inoc_dose_IB_mL | Volume administered intrabronchially. |
| inoc_dose_IT_pfu | Total pfu administered intratracheally. |
| inoc_dose_IT_TCID | Total TCID50 administered intratracheally. |
| inoc_dose_IN_pfu | Total pfu administered intranasally. |
| inoc_dose_IN_TCID | Total TCID50 administered intranasally. |
| inoc_dose_OR_pfu | Total pfu administered orally. |
| inoc_dose_OR_TCID | Total TCID50 administered orally. |
| inoc_dose_OC_pfu | Total pfu administered ocularly/conjunctivally. |
| inoc_dose_OC_TCID | Total TCID50 administered ocularly/conjunctivally. |
| inoc_dose_IG_pfu | Total pfu administered intragastrically. |
| inoc_dose_IG_TCID | Total TCID50 administered intragastrically. |
| inoc_dose_IV_pfu | Total pfu administered intravenously. |
| inoc_dose_IV_TCID | Total TCID50 administered intravenously. |
| inoc_dose_IB_pfu | Total pfu administered intrabronchially. |
| inoc_dose_IB_TCID | Total TCID50 administered intrabronchially. |
| inoc_dose_AE_pfu | Total estimated inhaled pfu during aerosol exposure. |
| inoc_dose_AE_TCID | Total estimated inhaled TCID50 during aerosol exposure. |
| article_doi | DOI of the publication. |
| article_title | Full article title. |
| article_journal | Journal in which the study was published. |
| preprint | Whether the included study is a preprint. |
| indiv_rep | Original animal identifier reported by the study. |
| data_availability | How the raw data were obtained. |
| cell_line | Cell line used for viral culture assays. |
| culture_assay | Culture assay measurement type (e.g., pfu or TCID50). |
| location_idx | Numeric encoding of the study location. |
| study_location | Institution or facility where the NHP study was performed. |