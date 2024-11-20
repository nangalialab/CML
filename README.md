### Files to reproduce figures and tables in the paper "Timing the Philadelphia chromosome and trajectory to chronic myeloid leukaemia"
This produces some figures and files that are used in downstream analysis.

### Dependencies
See install_rpackages.R for a list of required packages - also see export/sessionInfo.txt for additional information about the environment and package versions.

### Steps to populate export directory 
Starting from trees with assigned mutations (../cache/PDD.RDS) the following produces tree, selection and timing figures together with associated tables in the export directory.  The process will take a few hours to run.  
```{bash}
cd phyloanalysis
Rscript run_all.R
```

Previously generated versions of the plots and tables are in figures/v1 and export/v1

