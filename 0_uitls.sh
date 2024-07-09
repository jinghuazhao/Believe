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

options(width=200)
library(AnnotationDbi)
suppressMessages(library(EnsDb.Hsapiens.v75))
suppressMessages(library(SomaScan.db))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
suppressMessages(library(dplyr))
library(parallel)
library(valr)
keys(SomaScan.db)
keytypes(SomaScan.db)
columns(SomaScan.db)
keys(EnsDb.Hsapiens.v75)[1:10L]
grep("ENSEMBL", columns(SomaScan.db), value = TRUE)
columns(EnsDb.Hsapiens.v75)
pos_sel <- SomaScan.db::select(SomaScan.db, "11138-16", columns = c("PMID","ENTREZID","SYMBOL","ENSEMBL"))
pos_col <- c("GENEID","GENESEQEND","GENESEQSTART","TXNAME","SEQNAME","TXSEQSTART","TXSEQEND","UNIPROTID")
pos_res <- SomaScan.db::select(EnsDb.Hsapiens.v75, keys = "ENSG00000020633",columns = pos_col)
m <- merge(pos_sel, pos_res, by.x = "ENSEMBL", by.y = "GENEID")
mm <- dplyr::rename(m,chrom=SEQNAME,start=TXSEQSTART,end=TXSEQEND) %>%
      valr::bed_merge()
numCores <- 8
cl <- makeCluster(numCores)
clusterEvalQ(cl, library(SomaScan.db))
clusterEvalQ(cl, library(EnsDb.Hsapiens.v75))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(valr))
results <- data.frame()
results <- mclapply(keys(SomaScan.db), function(key) {
  pos_sel <- SomaScan.db::select(SomaScan.db, keys = key,
                             columns = c("PMID", "ENTREZID", "SYMBOL", "ENSEMBL"))
  if (!is.null(pos_sel$ENSEMBL) & any(!is.na(pos_sel$ENSEMBL))) {
    ensembl_ids <- pos_sel$ENSEMBL[!is.na(pos_sel$ENSEMBL)]
    sapply(unique(ensembl_ids), function(ensembl_id) {
      pos_res <- SomaScan.db::select(EnsDb.Hsapiens.v75, keys = ensembl_id,
                                 columns = c("GENEID", "SEQNAME", "TXSEQSTART", "TXSEQEND"))
      merged_res <- merge(pos_sel[pos_sel$ENSEMBL == ensembl_id, ], pos_res,
                          by.x = "ENSEMBL", by.y = "GENEID", all = TRUE)
      genenames <- unique(merged_res$SYMBOL)
      regions <- merged_res %>%
                 dplyr::rename(chrom=SEQNAME,start=TXSEQSTART,end=TXSEQEND) %>%
                 valr::bed_merge()
      data.frame(genename=genenames,region=regions)
    })
  } else {
    NULL
  }
}, mc.cores = numCores)
results <- do.call(rbind, results)
stopCluster(cl)
save(results, file="work/SomaScan.RData")
END
