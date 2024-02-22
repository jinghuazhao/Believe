# Believe analysis

- Web: [https://jinghuazhao.github.io/Believe](https://jinghuazhao.github.io/Believe).
- CSD3 location: ~/rds/[rds-jmmh2-]post_qc_data/believe
    * the BELIEVE TOPMed imputation data by Regeneron for the full ~72,000 participants.
~/rds/[rds-jmmh2-post]_qc_data/believe/genotype/imputed where the new files have
‘Freeze_Two’ in the name instead of ‘Freeze_One’.
    * PCs for the full dataset are in place (counterpart to the sequenced variants,
~/rds/[rds-jmmh2-]post_qc_data/believe/genotype/genomewide/plink/aug_2023??).

A growing list of scripts.

```
1_desc.sh
```

## 1. Descriptive analysis

`1_desc.sh` is to extract data and to protoype some `canonical options` for handling relatedness in the study.
