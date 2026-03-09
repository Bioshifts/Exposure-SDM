# Species range shifts often speed ahead of their modeled climatic niches 

This repository contains the code required to reproduce the analyses presented in the study:

Oliveira et al. (2026). *Species range shifts often speed ahead of their modeled climatic niches.* PNAS.

The study evaluates how well climate-based species distribution models (SDMs) predict the direction and rate of observed species range shifts.

## Repository Structure 

```         
├── R/                # R scripts for data processing, modeling, and analysis   
├── Data/             # Raw and processed data (occurrences, environmental layers)
├── Data/Output/      # Processed dataset (merged documented and modeled range shifts)
├── renv/             # renv project library
├── renv.lock         # Reproducible R package environment
├── README.md         # Project documentation
```

## Workflow

1.  **Select species with documented range shifts in the BioShifts and CoRE databases.**

2.  **Download species occurrence data from GBIF**

3.  **Download environmental data**

    -   Terrestrial species: CHELSAcruts climate datasets

    -   Marine species: ORAS oceanographic datasets

<!-- -->

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

    -   CHELSAcruts for terrestrial species ([http://chelsa-climate.org)](http://chelsa-climate.org))

    -   ORAS for marine species (<https://www.oceandata.org>)

> ⚠️ Raw data may require downloading separately due to size or licensing.

## Citation

If you use this code or build upon this work, please cite:

Oliveira, B. F. et al. (2026). *Species range shifts often speed ahead of their modeled climatic niches.* Proceedings of the National Academy of Sciences (PNAS).
