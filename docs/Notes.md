# Notes

Examples and resources which can be seen from this page and <https://github.com/jinghuazhao/Believe/tree/master/docs/T2D>.

## BGLR/

This folder contains data and R script for Crossa et al. (2020), see below.

## brms/

At this stage, the demonstrate how to run a SLURM job on CSD3.

Change account and $USER appropriately.

- Script, **csd3.sb**
- Error report: **_csd3.e**
- Output report: **_csd3.o**
- Plot: **brms-fit3.pdf**

## ldsc

This implements LD-score regression, which is available from CSD3 by `module load ceuadmin/ldsc` as indicated here, <https://cambridge-ceu.github.io/csd3/Python/ldsc.html>.

## Mixed effects Cox models

- **kinship**, in [here](https://github.com/jinghuazhao/jinghuazhao.github.io/tree/master/docs/software), CRAN [archive](https://cran.r-project.org/src/contrib/Archive/kinship/) and [GitHub](https://github.com/cran/kinship).
- **coxme**, <https://cran.r-project.org/web/packages/coxme/index.html>, which includes updates on **kinship**.
- **coxmeg**, <https://github.com/lhe17/coxmeg>, which is an update to CRAN [archive](https://cran.r-project.org/src/contrib/Archive/coxmeg/).

All are R packages available under `ceuadmin/R` modules.

## References

Crossa J, Montesinos-López OA, Pérez-Rodríguez P, Costa-Neto G, Fritsche-Neto R, Ortiz R, Martini JWR, Lillemo M, Montesinos-López A, Jarquin D, Breseghello F, Cuevas J, Rincent R. Genome and Environment Based Prediction Models and Methods of Complex Traits Incorporating Genotype × Environment Interaction. Methods Mol Biol. 2022;2467:245-283. doi: 10.1007/978-1-0716-2205-6_9. PMID: 35451779,  <https://link.springer.com/protocol/10.1007/978-1-0716-2205-6_9>.

He L, Kulminski AM. Fast Algorithms for Conducting Large-Scale GWAS of Age-at-Onset Traits Using Cox Mixed-Effects Models. Genetics. 2020 May;215(1):41-58. doi: 10.1534/genetics.119.302940. Epub 2020 Mar 4. Erratum in: Genetics. 2020 Aug;215(4):1191. PMID: 32132097; PMCID: PMC7198273, <https://academic.oup.com/genetics/article/215/1/41/5930490>.
