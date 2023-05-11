# danube_indicators

Analysis of microbiome data for water quality status classification using machine learning. Using microbial metabarcoding data and machine learning we identified indicators for water quality status assessement in the Danube river. This repository includes scripts to run cutadapt (https://github.com/marcelm/cutadapt), dada2 (https://github.com/benjjneb/dada2) and statistical analyses for developing classifiers for water quality ased on microbial metabarcoding data.

## pull repository

```
cd path to repositories
git clone github.com/alper1976/danube_indicators.git
```

## Authors and acknowledgment
Scripts were written by Alexander Eiler, Laurent Fontaine and Jing Wei.

## License
This Code is subject to the terms of the MIT License. 

## Project status
Results from this project have been submitted to a peer-reviewed scientifc journal.

## Folders and code
The analyses code is divided into two folders "metags" and "rRNA_amplicons" representing the code to analyze whole shotgun metagenomic and SSU rRNA gene amplicon data, respectively.

### metadata.R
This represents the R code to perform the statistical analysis and data visualization on the metadata such as nutrient and gas concentrations along the chronosequences.

## released version 1.0.0

