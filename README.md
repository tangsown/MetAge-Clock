# Fig. 5 source data and analysis

This repository contains the CALERIE2 and T2D participant-level score files and the R code used for Fig. 5. `CALERIE2_Fig5_participants.csv` includes MetAge and Light-MetAge trajectories. `T2D373_Fig5_participants.csv` includes Light-MetAge trajectories but excludes individual demographic and clinical covariates and 1-year HbA1c outcomes.

`Fig5_CALERIE2_T2D_GitHub.R` documents the longitudinal mixed models, follow-up contrasts, paired plots, and the Fig. 5G logistic regression. The complete script cannot be rerun from the public CSVs alone. The regression requires a restricted T2D CSV containing `Age_baseline`, `Sex`, `BMI_baseline`, `FPG_baseline_mmol_L`, and `HbA1c_target_lt6_5_1y`. The script checks for these fields and stops with an explanatory message when they are absent. It previews plots only in an interactive R session; file exports are commented out.

## Data Availability

The derived participant-level CALERIE2 data and selected T2D Light-MetAge trajectory data underlying Fig. 5 are available in this repository. Additional individual-level T2D data, including the covariates and 1-year HbA1c outcome required to reproduce the logistic-regression analysis, are not publicly available because of participant privacy and ethical restrictions. Researchers seeking access should apply to the Ethics Committee of the First Affiliated Hospital, Sun Yat-Sen University. Following ethics approval, the data may be obtained from the corresponding author in accordance with the approved request.
