# Cross-Lagged Panel Models for Temporal Causality

In this folder there are the scripts used to perform all the simulations for the cross lagged panel models (CLPMs) in different scenarios:
- simulated linear data with 4 time points and different betas and standard errors with linear (`simulations_linear.R`) or quadratic (`simulations_quadratic.R`) dependance from time
- the same scenario with missing data, at random or not (`simulations_missing.R`)
- the same scenario adding a covariate (`simulations_cov_testing.R`, `simulations_cov_without_testing.R`, `simulations_cov_larger_testing.R`, `simulations_cov_larger_without_testing.R`)
- non-linear data with CLPMs with linear regression (`simulations_linear_on_nonLinear.R`,`simulations_linear_on_nonLinear_with_t.R`)
- non-linear data with CLPMs with Generalized Additive Models (GAMs) (`simulations_nonLinear.R`,`simulations_nonLinear_with_t.R`)


All the analysis of the simulated data were performed in the long format and in the wide format and were run using the bash file `run_all.sh`. The script called `utilities.R` contains the functions used in all the other scripts.

Finally, we generated the tables for the results related to the supplementary material of the article using the script `final_tables.R`.
