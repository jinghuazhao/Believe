#!/usr/bin/bash

export believe=~/rds/rds-jmmh2-post_qc_data/believe
export believe=~/rds/post_qc_data/believe
export X=CAMBRIDGE-BELIEVE_Freeze_One.GxS.chrX.vcf.gz
export bfile=${believe}/genotype/genomewide/plink/may_2022/CAMBRIDGE-BELIEVE_Freeze_One.GxS

export BELIEVEdata_GeneticIDMapping_P5031_20210826=BELIEVEdata_GeneticIDMapping_P5031_20210826.csv
export BELIEVEdata_P5031_20210826=BELIEVEdata_P5031_20210826.csv

export BELURBANdata_GeneticIDMapping_P5032_20210826=BELURBANdata_GeneticIDMapping_P5032_20210826.csv
export BELURBANdata_P5032_20210826=BELURBANdata_P5032_20210826.csv

Rscript -e '
   options(width=200)
   believe <- Sys.getenv("believe")
   idfile <- file.path(believe,"phenotype",Sys.getenv("BELIEVEdata_GeneticIDMapping_P5031_20210826"))
   id <- read.csv(idfile)
   head(id)
   datfile <- file.path(believe,"phenotype",Sys.getenv("BELIEVEdata_P5031_20210826"))
   dat <- read.csv(datfile)
   vars <- names(dat)
   ids <- c("HseIdentifier","Identifier","Gender","age","ethnic","religion")
   bmi <- vars[grep("BMI|Height|Weight|bmi|height2|weight2",vars)]
   t2d <- vars[grep("Type2|AntiDiabetics",vars)]
   vars_to_use <- c(ids,bmi,t2d)
   dat_to_use <- dat[vars_to_use]
   head(dat_to_use)
   idfile_urban <- file.path(believe,"phenotype",Sys.getenv("BELURBANdata_GeneticIDMapping_P5032_20210826"))
   id_urban <- read.csv(idfile_urban)
   head(id_urban)
   names(id_urban)
   datfile_urban <- file.path(believe,"phenotype",Sys.getenv("BELURBANdata_P5032_20210826"))
   dat_urban <- read.csv(datfile_urban)
   head(dat_urban)
   vars_urban <- names(dat_urban)
'
