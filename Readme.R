##  Copyright 2023 Annamaria Guolo (University of Padova) ##

Replication material for the manuscript 

Guolo A. (2022). Approximate likelihood and pseudo-likelihood inference in meta-analysis of diagnostic accuracy studies accounting for disease prevalence and study's design. Submitted.


 1. DTMAprevalence_software.R file for meta-analysis of accuracy studies in the R programming language;

 2. DTMAprevalence_auxiliary_functions.c file needed to run the exact likelihood estimation. In case you want to evaluate the exact likelihood, type "R CMD SHLIB DTMAprevalence_auxiliary_functions.c" in your terminal to compile the file, and upload the resulting .so file in your R session, "dyn.load('DTMAprevalence_auxiliary_functions.so')".

 3. DTMAprevalence_data_application.R to replicate the beta-D-glucane data analysis included in the manuscript (Karageorgopoulos et al., Clin. Infect. Dis., 2011).


Please note:
the number of nodes and weights for the Gauus-Hermite quadrature needed for the exact likelihood and the pseudo-likelihood are fixed at 21. The number can be changed at user preference.


March 2023, Padova, Italy.
