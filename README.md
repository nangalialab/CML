### Files to reproduce figures and tables in the paper "Timing and trajectory of BCR::ABL1 driven chronic myeloid leukaemia"
This produces some figures and files that are used in downstream analysis.

### Dependencies
See install_rpackages.R for a list of required packages - also see export/sessionInfo.txt for additional information about the environment and package versions.

### Steps to populate export directory 
Starting from trees with assigned mutations (../cache/PDD.RDS) the following produces tree, selection and timing figures together with associated tables in the export directory.  The process will take a few hours to run.  
```{bash}
cd phyloanalysis
Rscript run_all.R
```

Previously generated versions of the plots, tables and markdown HTMLs are in figures/v1 and export/v1 and mixedmodels/v1

The above takes about 6 hours to run on a 4 core CPU (Intel(R) Xeon(R) Gold 6226R CPU @ 2.90GHz).  The code has been tested on mac0S Sequoia and Ubuntu Jammy.   

A manifest indicating the QC status of samples available on EGA is in the "data/EGAD00001015353.sample_manifest.n1023.csv"

