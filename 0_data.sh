#!/usr/bin/bash

export src=~/rds/post_qc_data/believe
export X=CAMBRIDGE-BELIEVE_Freeze_One.GxS.chrX.vcf.gz
export bfile=${src}/genotype/genomewide/plink/may_2022/CAMBRIDGE-BELIEVE_Freeze_One.GxS
export ID=${src}/phenotype/BELIEVEdata_GeneticIDMapping_P5031_20210826.csv
export dat=${src}/phenotype/BELIEVEdata_P5031_20210826.csv

Rscript -e '
   options(width=200)
   datfile <- Sys.getenv("dat")
   dat <- read.csv(datfile)
   vars <- names(dat)
   id <- c("HseIdentifier","Identifier","ethnic","religion")
   bmi <- vars[grep("BMI|bmi",vars)]
   t2d <- vars[grep("Type2",vars)]
   vars_to_use <- c(id,bmi,t2d)
   dat_to_use <- dat[vars_to_use]
   head(dat_to_use)
'
