#!/usr/bin/bash

module load ceuadmin/R

export pre_qc=~/rds/pre_qc_data/believe
export nightingale=${pre_qc}/nightingale_nmr
export olink=${pre_qc}/olink_proteomics
export somalogic=${pre_qc}/somalogic_proteomics

Rscript -e '
# Olink
Sys.setenv(LIBARROW_MINIMAL = "false") # all optional features including gzip
Sys.setenv(ARROW_WITH_GZIP = "ON") # only gzip
library(arrow)
library(OlinkAnalyze)
olink <- Sys.getenv("olink")
npx <- read_parquet(file.path(olink,"Q-08620_Di_Angelantonio_NPX_2023-12-27.parquet"))
table(npx$SampleType)
npx[c("SampleID","OlinkID","Panel","PCNormalizedNPX","Count","NPX","UniProt","AssayQC","SampleQC")]
# Q-08620_Di_Angelantonio_Analysis_Report_2023-12-27.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Explore_HT_SampleForm_v2.0.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Sample_manifest_-_plate_Explore_HT_20231120.xlsx

if (FALSE)
{
  npx_data <- read_NPX(filename = file.path(olink,"some.xlsx"))
  olink_dist_plot(npx_data,
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c('turquoise3', 'red'))
  bridge_samples <- intersect(x = npx_data1$SampleID, y = npx_data2$SampleID)
  bridge_normalized_data <- olink_normalization(df1 = npx_data1,
                            df2 = npx_data2,
                            overlapping_samples_df1 = bridge_samples,
                            df1_project_nr = "20200001",
                            df2_project_nr = "20200002",
                            reference_project = "20200001")
  ttest_results <- olink_ttest(df = npx_data, variable = "Treatment")
  ttest_sign <- ttest_results %>%
      head(n=10) %>%
      pull(OlinkID)
  olink_volcano_plot(p.val_tbl = ttest_results,
    olinkid_list = ttest_sign) +
    scale_color_manual(values = c('turquoise3', 'red'))
}

# SomaLogic
suppressMessages(library(SomaDataIO))
somalogic <- Sys.getenv("somalogic")
# batch 1, SS-2218747_SQS.pdf
adat_samples <- read.delim(file.path(somalogic,"Batch1","Batch1_adat_samplenames.tsv"))
snmlSMP_samples <- read.delim(file.path(somalogic,"Batch1","Batch1_adat_samplenames_anmlSMP.tsv"))
SomaSignals <- read.delim(file.path(somalogic,"Batch1","SS-2218747_SomaSignals_20220916_v4.tsv"),skip=3)
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch1",
                              "SS-2218747_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch1",
                              "SS-2218747_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

# Batch 2, SS-2225082_SQS.pdf
adat_samples <- read.delim(file.path(somalogic,"Batch2","Batch2_adat_samplenames.tsv"))
snmlSMP_samples <- read.delim(file.path(somalogic,"Batch2","Batch2_adat_samplenames_anmlSMP.tsv"))
SomaSignals <- read.delim(file.path(somalogic,"Batch2","SS-2225082_SomaSignals_20220916_v4.tsv"),skip=3)
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch2",
                              "SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch2",
                              "SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

# Batch1_Batch2_combined, SS-2218747_SS-2225082_SQS.pdf
adat <- SomaDataIO::read_adat(file.path(somalogic,"Batch1_Batch2_combined",
        "SS-2218747_SS-2225082_Combined_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.adat"))
anmlSMP <- SomaDataIO::read_adat(file.path(somalogic,"Batch1_Batch2_combined",
           "SS-2218747_SS-2225082_Combined_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"))

options(width=200)
library(dplyr)
library(Biobase)
class(adat)
dim(adat)
methods(class="soma_adat")
table(adat$SampleType)

cleanData <- adat |>
  filter(SampleType == "Sample") |>
  log10()
getAnalyteInfo(cleanData) |>
  select(AptName, SeqId, Target = TargetFullName, EntrezGeneSymbol, UniProt)

adat_eset <- SomaDataIO::adat2eSet(cleanData)
clinData <- pData(adat_eset)

# SSM-00060 - Rev 1.0 - Data Standardization and File Specification Technical Note (002).pdf'
# Test Information Guides Plasma PAV/
# Test Information Guides Plasma SST/

position <- function()
{
  suppressMessages(library(SomaScan.db))
  keys(SomaScan.db)
  keytypes(SomaScan.db)
  columns(SomaScan.db)
  select(SomaScan.db, keys = "18342-2", columns = c("SYMBOL", "UNIPROT"))
  pos_sel <- select(SomaScan.db, "11138-16", columns = c("SYMBOL", "ENTREZID", "ENSEMBL"))
  keys(EnsDb.Hsapiens.v75)[1:10L]
  grep("ENSEMBL", columns(SomaScan.db), value = TRUE)
  columns(EnsDb.Hsapiens.v75)
  pos_res <- select(EnsDb.Hsapiens.v75, keys = "ENSG00000020633",
                    columns = c("GENEBIOTYPE", "SEQCOORDSYSTEM", "GENEID",
                                "PROTEINID", "PROTDOMSTART", "PROTDOMEND"))
  merge(pos_sel, pos_res, by.x = "ENSEMBL", by.y = "GENEID")
}
