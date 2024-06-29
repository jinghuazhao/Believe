#!/usr/bin/bash

module load ceuadmin/R

export pre_qc=~/rds/pre_qc_data/believe
export nightingale=${pre_qc}/nightingale_nmr
export olink=${pre_qc}/olink_proteomics
export somalogic=${pre_qc}/somalogic_proteomics

Rscript -e '
library(arrow)
olink <- Sys.getenv("olink")
npx <- read_parquet(file.path(olink,"Q-08620_Di_Angelantonio_NPX_2023-12-27.parquet"))
# Q-08620_Di_Angelantonio_Analysis_Report_2023-12-27.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Explore_HT_SampleForm_v2.0.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Sample_manifest_-_plate_Explore_HT_20231120.xlsx

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

# SSM-00060 - Rev 1.0 - Data Standardization and File Specification Technical Note (002).pdf'
# Test Information Guides Plasma PAV/
# Test Information Guides Plasma SST/

position <- function()
{
  options(width=200)
  suppressMessages(library(SomaScan.db))
  pos_sel <- select(SomaScan.db, "11138-16", columns = c("SYMBOL", "ENTREZID", "ENSEMBL"))
  keys(EnsDb.Hsapiens.v75)[1:10L]
  grep("ENSEMBL", columns(SomaScan.db), value = TRUE)
  columns(EnsDb.Hsapiens.v75)
  pos_res <- select(EnsDb.Hsapiens.v75, keys = "ENSG00000020633",
                    columns = c("GENEBIOTYPE", "SEQCOORDSYSTEM", "GENEID",
                                "PROTEINID", "PROTDOMSTART", "PROTDOMEND"))
  merge(pos_sel, pos_res, by.x = "ENSEMBL", by.y = "GENEID")
}

if (FALSE)
{
  Sys.setenv(LIBARROW_MINIMAL = "false") # all optional features including gzip
  Sys.setenv(ARROW_WITH_GZIP = "ON") # only gzip
  install.packages("arrow")
  install_arrow()
  library(OlinkAnalyze)
  npx_data <- read_NPX(filename = file.path(olink,
                                  "Q-08620_DiAngelantonio_-_Olink_-_Sample_manifest_-_plate_Explore_HT_20231120.xlsx"))
  # visualize the NPX distribution per sample per panel, example for one panel
  olink_dist_plot(npx_data %>%
    filter(Panel == 'Olink CARDIOMETABOLIC')) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = c('turquoise3', 'red'))
  olink_qc_plot(npx_data1 %>% filter(Panel == 'Olink CARDIOMETABOLIC')) +
    scale_color_manual(values = c('turquoise3', 'red'))
  # identify bridge samples
  bridge_samples <- intersect(x = npx_data1$SampleID,
                 y = npx_data2$SampleID)
  # bridge normalize
  bridge_normalized_data <- olink_normalization(df1 = npx_data1,
                            df2 = npx_data2,
                            overlapping_samples_df1 = bridge_samples,
                            df1_project_nr = "20200001",
                            df2_project_nr = "20200002",
                            reference_project = "20200001")
  # t-test npx_data
  ttest_results_NPX1 <- olink_ttest(df = npx_data, variable = "Treatment")
  # select names of the top #10 most significant proteins
  ttest_sign_NPX1 <- ttest_results_NPX1 %>%
      head(n=10) %>%
      pull(OlinkID)
  # volcano plot with annotated top #10 most significant proteins
  olink_volcano_plot(p.val_tbl = ttest_results_NPX1,
    olinkid_list = ttest_sign_NPX1) +
    scale_color_manual(values = c('turquoise3', 'red'))
}
'
