#!/usr/bin/bash

module load ceuadmin/R

export pre_qc=~/rds/pre_qc_data/believe
export nightingale=${pre_qc}/nightingale_nmr
export olink=${pre_qc}/olink_proteomics
export somalogic=${pre_qc}/somalogic_proteomics

R --no-save <<END
options(width=200)
suppressMessages(library(dplyr))
# Olink
Sys.setenv(LIBARROW_MINIMAL = "false") # all optional features including gzip
Sys.setenv(ARROW_WITH_GZIP = "ON") # only gzip
library(arrow)
library(OlinkAnalyze)
olink <- Sys.getenv("olink")
npx <- read_parquet(file.path(olink,"Q-08620_Di_Angelantonio_NPX_2023-12-27.parquet"))
table(with(npx,SampleType))
npx[c("SampleID","OlinkID","Panel","PCNormalizedNPX","Count","NPX","UniProt","AssayQC","SampleQC")] %>%
filter(NPX!=0)
# Q-08620_Di_Angelantonio_Analysis_Report_2023-12-27.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Explore_HT_SampleForm_v2.0.pdf
# Q-08620_DiAngelantonio_-_Olink_-_Sample_manifest_-_plate_Explore_HT_20231120.xlsx

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

library(Biobase)
class(adat)
dim(adat)
methods(class="soma_adat")
table(with(adat,SampleType))

cleanData <- adat |>
  filter(SampleType == "Sample") |>
  log10()
getAnalyteInfo(cleanData) |>
  select(AptName, SeqId, Target = TargetFullName, EntrezGeneSymbol, UniProt)

adat_eset <- SomaDataIO::adat2eSet(cleanData)
protData <- exprs(adat_eset)
clinData <- pData(adat_eset)

# SSM-00060 - Rev 1.0 - Data Standardization and File Specification Technical Note (002).pdf
# Test Information Guides Plasma PAV/
# Test Information Guides Plasma SST/

position <- function()
{
  options(width=200)
  suppressMessages(library(SomaScan.db))
  suppressMessages(library(EnsDb.Hsapiens.v75))
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  keys(SomaScan.db)
  keytypes(SomaScan.db)
  columns(SomaScan.db)
  select(SomaScan.db, keys = "18342-2", columns = c("ENTREZID","SYMBOL","UNIPROT"))
  keys(EnsDb.Hsapiens.v75)[1:10L]
  grep("ENSEMBL", columns(SomaScan.db), value = TRUE)
  columns(EnsDb.Hsapiens.v75)
  pos_sel <- select(SomaScan.db, "11138-16", columns = c("PMID","ENTREZID","SYMBOL","ENSEMBL"))
  pos_res <- select(EnsDb.Hsapiens.v75, keys = "ENSG00000020633",
                    columns = c("GENEID","TXNAME","SEQNAME","TXSEQSTART","TXSEQEND"))
  merge(pos_sel, pos_res, by.x = "ENSEMBL", by.y = "GENEID")
# loop-over
  results <- data.frame()
  for (key in keys(SomaScan.db)) {
    pos_sel <- select(SomaScan.db, keys = key, columns = c("PMID", "ENTREZID", "SYMBOL", "ENSEMBL"))
    if (!is.null(pos_sel$ENSEMBL) && any(!is.na(pos_sel$ENSEMBL))) {
      ensembl_ids <- pos_sel$ENSEMBL[!is.na(pos_sel$ENSEMBL)]
      for (ensembl_id in unique(ensembl_ids)) {
        pos_res <- select(EnsDb.Hsapiens.v75, keys = ensembl_id,
                          columns = c("GENEID", "TXNAME", "SEQNAME", "TXSEQSTART", "TXSEQEND"))
        merged_res <- merge(pos_sel[pos_sel$ENSEMBL == ensembl_id, ], pos_res,
                            by.x = "ENSEMBL", by.y = "GENEID", all = TRUE)
        results <- rbind(results, merged_res)
      }
    }
  }
  save(results,file="work/SomaScan.rda")
}

position()
END
