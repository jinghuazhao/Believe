# Believe analysis

Web: <https://jinghuazhao.github.io/Believe>

## CSD3 location

- `~/rds/pre_qc_data/believe`
- `~/rds/rds-jmmh2-post_qc_data/believe`
- `~/rds/rds-jmmh2-post_qc_data/believe/reference_files/genetic/pLOF/pLOF_annotation_may_2022`

- TOPMed imputation data by Regeneron for the full ~72k participants. `genotype/imputed/*Freeze_Two.GxS.TOPMED_dosages*` (.pgen-GT files, .pgen-HDS files and .bgen files pending).
- PCs (counterpart to the sequenced variants, `genotype/genomewide/plink/aug_2023`??).

A growing list of scripts.

```
0_utils.sh
1_desc.sh
2_gcta.sb
```

## 0. Utilities

`0_utils.sh` reads information such as Olink Parquet and SomaLogic adat files.

## 1. Descriptive analysis

`1_desc.sh` is to extract data and to protoype some `canonical options` for handling relatedness in the study.

## 2. GCTA work

`2_gcta.sb` contains various approaches for the analysis.

## A. URLs

- Olink, <https://github.com/Olink-Proteomics>, <https://cran.r-project.org/web/packages/OlinkAnalyze/index.html>
- SomaLogic, <https://github.com/somalogic>
- BELIEVE study, <https://www.phpc.cam.ac.uk/ceu/international-vascular-health/believe/>
- BLR_GxE_Analysis, <https://github.com/angellicamarta/BLR_GxE_Analysis>
