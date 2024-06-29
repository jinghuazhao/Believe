#!/usr/bin/bash

module load ceuadmin/R

export pre_qc=~/rds/pre_qc_data/believe
export nightingale=${pre_qc}/nightingale_nmr
export olink=${pre_qc}/olink_proteomics
export somalogic=${pre_qc}/somalogic_proteomics

Rscript -e '
library(arrow)
olink <- Sys.getenv("olink")
npx <- read_parquet("Q-08620_Di_Angelantonio_NPX_2023-12-27.parque",
                    n_threads = 4,
                    compression = "gzip")

# Q-08620_Di_Angelantonio_Analysis_Report_2023-12-27.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Explore_HT_SampleForm_v2.0.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Sample_manifest_-_plate_Explore_HT_20231120.xlsx

library(SomaDataIO)
somalogic <- Sys.getenv("somalogic")
# batch 1, SS-2218747_SQS.pdf
adat_samples <- read.delim(file.path(somalogic,"Batch1","Batch1_adat_samplenames.tsv"))
snmlSMP_samples <- read.delim(file.path(somalogic,"Batch1","Batch1_adat_samplenames_anmlSMP.tsv"))
SomaSignals <- read.delim(file.path(somalogic,"Batch1","SS-2218747_SomaSignals_20220916_v4.tsv"))
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch1",
                              "SS-2218747_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch2",
                              "SS-2218747_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

# Batch 2, SS-2225082_SQS.pdf
adat_samples <- read.delim(file.path(somalogic,"Batch2","Batch2_adat_samplenames.tsv"))
snmlSMP_samples <- read.delim(file.path(somalogic,"Batch2","Batch2_adat_samplenames_anmlSMP.tsv"))
SomaSignals <- read.delim(file.path(somalogic,"Batch12,"SS-2225082_SomaSignals_20220916_v4.tsv"))
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch2",
                              "SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch2",
                              "SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

# Batch1_Batch2_combined, SS-2218747_SS-2225082_SQS.pdf
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch1_Batch2_combined",
        "SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch1_Batch2_combined",
           "SS-2218747_SS-2225082_Combined_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

# SSM-00060 - Rev 1.0 - Data Standardization and File Specification Technical Note (002).pdf'
# Test Information Guides Plasma PAV/
# Test Information Guides Plasma SST/
'
