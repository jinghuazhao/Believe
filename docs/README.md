# Believe analysis

Web: <https://jinghuazhao.github.io/Believe>

CSD3 location: `~/rds/[rds-jmmh2-]post_qc_data/believe`

- TOPMed imputation data by Regeneron for the full ~72k participants. `genotype/imputed/*` with `Freeze_Two` in the name (`*_Freeze_Two.GxS.TOPMED_dosages*`, including the .pgen file for the GT version, with the .pgen file for the HDS version and the .bgen file pending.) instead of `Freeze_One`.
- PCs in place (counterpart to the sequenced variants, `genotype/genomewide/plink/aug_2023`??).

A growing list of scripts.

```
1_desc.sh
```

## 1. Descriptive analysis

`1_desc.sh` is to extract data and to protoype some `canonical options` for handling relatedness in the study.
