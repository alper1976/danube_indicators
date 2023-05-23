# danube_indicators

Analysis of microbiome data for water quality status classification using machine learning. Using microbial metabarcoding data and machine learning we identified indicators for water quality status assessement in the Danube river. This repository includes scripts to run cutadapt (https://github.com/marcelm/cutadapt), dada2 (https://github.com/benjjneb/dada2) and statistical analyses for developing classifiers for water quality based on microbial metabarcoding data.

## pull repository

```
cd path to repositories
git clone github.com/alper1976/danube_indicators.git
```

## Authors and acknowledgment
Scripts were written by Laurent Fontaine, Lorenzo Pin and Alexander Eiler.

## License
This Code is subject to the terms of the MIT License. 

## Project status
Results from this project have been submitted to a peer-reviewed scientifc journal.

## Folders and code
The analyses code is divided into two folders "raw_data" and "stats" representing the code to analyze raw sequencing data and perform statistical analysis, respectively. The "Data" folder contains the data input and output of the spatiotemporal analysis.

### metadata.tsv
This represents the data from the ICPDR that were used in this manuscript. Data was obtained from http://www.icpdr.org/wq-db

### ASV_table.tsv
This represents the Amplicon sequence variants (ASV) table obtained from the raw sequencing data using cutadapt and dada2.

### Taxonomy.tsv
This represents the taxonomic annotations for the amplicon sequence variants (ASV) obtained from the raw sequencing data using cutadapt and dada2.

## released version 1.0.0

