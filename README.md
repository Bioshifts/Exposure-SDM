# Species range shifts often speed ahead of their modeled climatic niches

This repository contains the code required to reproduce the analyses presented in the study:

Oliveira et al. (2026). *Species range shifts often speed ahead of their modeled climatic niches.* PNAS.

The study evaluates how well climate-based species distribution models (SDMs) predict the direction and rate of observed species range shifts.

## Repository Structure

```         
├── R/                               # R scripts for data processing, modeling, and analysis

│   ├── 0_supplementary/              # Scripts used for supplementary analyses
│   │   ├── 1_get_sps_for_supp.R      # Select random 10% of terrestrial and marine species
│   │   ├── 2_get_bioclimatics_sps.R  # Extract bioclimatic variables for each species
│   │   ├── 2_job_bioclimatics_sps.R  # Submit previous script to the HPC
│   │   ├── 3_run_sdm.R               # Run SDMs for supplementary species
│   │   ├── 3_job_sdm.R               # Submit SDM script to the HPC
│   │   ├── 4_gather_sdm_CV.R         # Gather model cross-validation results
│   │   ├── 5_get_shifts_edge_method.R# Calculate species range shifts from SDMs
│   │   ├── 5_job_shifts.R            # Submit range shift script to the HPC
│   │   ├── 10_merge_all.qmd          # Merge modeled and documented range shifts
│   │   └── 10_analyses.qmd           # Analyses comparing documented and modeled shifts

│   ├── 1_select_species/             # Select species for downstream analyses
│   │   └── 1_get_species_from_SDMs.R

│   ├── 2_bioclimatic_data/           # Calculation of bioclimatic variables

│   │   ├── CHELSA/                   # Climate data for terrestrial species
│   │   │   ├── 1_CHELSA_download.R           # Download climatic variables from CHELSAcruts
│   │   │   ├── 2_CHELSA_bioclimatic_proj.R   # Calculate temperature and precipitation summaries
│   │   │   ├── 2_CHELSA_bioclimatic_proj.sh  # Submit previous script to the HPC
│   │   │   └── 2_loop_years_proj.sh          # Loop batch submission across years

│   │   └── ORAS/                     # Climate data for marine species
│   │       ├── 1_download_oras5.py        # Download ORAS5 data
│   │       ├── 2_unzip_and_save.R        # Unzip and store downloaded files
│   │       ├── 2_ORAS_job_unzip_save.sh  # Submit unzip script to the HPC
│   │       ├── 3_ORAS_bioclimatics_proj.R# Calculate SST statistics
│   │       └── 3_ORAS_job_proj.sh        # Submit previous script to the HPC

│   ├── 2_2_bioclimatic_SA/           # Bioclimatic variables for study areas
│   │   ├── get_bios_SA.R             # Calculate yearly bioclimatic variables for each study area
│   │   └── job_bios_SA.R             # Submit previous script to the HPC

│   ├── 3_gbif_data/
│   │   └── 1_Get_GBIF_data.R         # Download and filter GBIF occurrence records

│   ├── 4_bioclimatics_sps/
│   │   ├── get_bioclimatics_sps.R    # Extract bioclimatic variables for GBIF records
│   │   └── job_bioclimatics_sps.R    # Submit previous script to the HPC

│   ├── 5_sdms/
│   │   ├── 1_run_sdm.R               # Run species distribution models
│   │   ├── 1_job_sdm.R               # Submit SDM jobs to the HPC
│   │   └── 2_gather_sdm_CV.R         # Gather SDM cross-validation results

│   ├── 6_range_shifts/
│   │   ├── 1_get_shifts_edge_method.R# Calculate species range shifts from SDMs
│   │   ├── 1_job_shifts.R            # Submit range shift script to the HPC
│   │   └── 2_gather_shifts_edge_method.R # Compile range shifts across species

│   ├── 9_connectivity/
│   │   ├── 1_get_connectivity_SA.R   # Calculate connectivity for each study area
│   │   ├── 1_job_connectivity_SA.R   # Submit previous script to the HPC
│   │   ├── 2_get_connectivity_sps.R  # Calculate species-specific connectivity
│   │   └── 2_job_connectivity_sps.R  # Submit previous script to the HPC

│   └── 10_analyses/
│       ├── 1_merge_all.qmd           # Compile documented and modeled range shifts
│       ├── 2_data_exploration.qmd    # Exploratory data analysis
│       └── 3_SDM_models.qmd          # Multivariate statistical analyses


├── Data/
│   └── Output/
│       ├── Bioshifts_merge_Exposure_all.csv
│       │   # Dataset containing merged documented and modeled range shifts
│       └── Bioshifts_merge_Exposure_all_supp.csv
│           # Dataset used in supplementary analyses

├── renv/                             # Project-specific R package library
├── renv.lock                         # Reproducible R package environment
└── README.md                         # Project documentation
```

## Workflow

1.  **Select species with documented range shifts in the BioShifts and CoRE databases.**

2.  **Download species occurrence data from GBIF**

3.  **Download environmental data**

    -   Terrestrial species: CHELSAcruts climate datasets

    -   Marine species: ORAS oceanographic datasets

4.  **Fit ensemble species distribution models (SDMs)** for each species.

5.  **Estimate modeled range shifts** using SDM projections across time.

6.  **Compare modeled and documented range shifts** to assess agreement in direction and magnitude.

## Reproducibility

This project uses the **`renv` package** to ensure reproducibility of the R environment.

To recreate the computational environment:

```         
install.packages("renv") 

renv::restore()
```

This will install all package versions recorded in `renv.lock`.

## Data Availability

-   **Species occurrences:** GBIF (<https://www.gbif.org>)

-   **Documented range shifts:** BioShifts and CoRE databases

-   **Environmental data:**

    -   CHELSAcruts for terrestrial species (<http://chelsa-climate.org>)

    -   ORAS5 for marine species (<https://cds.climate.copernicus.eu/datasets/reanalysis-oras5?tab=overview>)

> ⚠️ Raw data may require downloading separately due to size or licensing.

## Citation

If you use this code or build upon this work, please cite:

Oliveira, B. F. et al. (2026). *Species range shifts often speed ahead of their modeled climatic niches.* Proceedings of the National Academy of Sciences (PNAS).
