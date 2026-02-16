---
editor_options: 
  markdown: 
    wrap: 72
---

In this folder you will find scripts used for the paper "Longitudinal
plasma proteomics along a menstrual cycle highlights the regulatory
effect of sex hormones on inflammatory and cardiometabolic circuits" by
Zhernakova et al., submitted.

The structure is as follows:

-   1.data_processing_and_QC/ - this folder contains the scripts used
    for data processing and QC:

    -   1.olink_bridge_normalization.Rmd - normalization and QC of Olink
        proteomics data

    -   2.phase_reclassification.py

    -   reclassification of menstrual phases based on hormonal profiles

    -   3.calculate_PRS/ - scripts used to calculate genetic score for
        each protein to use as a covariate

-   2.main_proteomics_longitudinal.R - the main function that contains
    association analysis presented in the paper and additional
    statistical analyses.

-   3.causal_inference_models - scripts used for simulation analyses of
    cross-lagged causality models. For details check the separate readme
    file located in that folder

-   4.causal_inference_real_data.R - script to apply the cross-lagged
    temporal causality model to real protein, hormone and phenotype data

-   5.network/ - scripts and data used for plotting the association
    network

    -   1.format_associations_for_network.R - prepare the edges and
        nodes tables based on the association results

    -   2.launch_network_shiny.R - run the shiny app to plot the network

    -   network_data/ - actual association results used in the published
        network

-   utils/ - utility functions used in the main script
    2.main_proteomics_longitudinal.R
