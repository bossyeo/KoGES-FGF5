# KoGES-FGF5-mQTL

R code for the manuscript

**Replication of a cis-methylation QTL at the *FGF5* blood-pressure locus in two Korean population cohorts**

The five scripts reproduce every result, table and figure in the manuscript from the KoGES source data.
No individual-level data are included in this repository.

## Analysis flow

| Script | Analysis | Manuscript output |
|---|---|---|
| `01_GWAS_36GenePanel.R` | Panel-restricted association analysis: markers within ±20 kb of 282 genes (36 symbols / prefixes) tested against 8 chronic diseases by cohort and sex (logistic regression) | Methods counts; association results; Supplementary Tables S1–S2; Supplementary Data 1 |
| `02_mQTL_AX83223998_FGF5.R` | Genome-wide mQTL of the lead marker in CITY (discovery) and ASAS (replication) with MatrixEQTL; fixed-effect inverse-variance meta-analysis | mQTL results; Table 2 |
| `03_FGF5_3CpG_validation.R` | Validation of the meta-analysis CpGs: probe SNPs, sensitivity analyses, cell-type (EpiDISH) and genetic-PC adjustment, partial R², CpG–hypertension association | Table 1; Supplementary Tables S3–S4; Supplementary Figure S1 |
| `04_GoDMC_lookup.R` | Look-up of the CpGs in the GoDMC blood mQTL database (requires internet) | GoDMC comparison in the Discussion |
| `05_Figures.R` | Figures 1–3 from the summary outputs of scripts 01 and 02 | Figures 1–3 |

Each script reads the results of the previous one, so the analysis targets are selected by the code itself:
script 02 takes the marker with the smallest P among tests with FDR < 0.05 from script 01,
and script 03 takes the CpGs passing Bonferroni correction in the meta-analysis from script 02.
Run the scripts in order: 01 → 02 → 03 → 04 → 05.

## Data

The data are from the Korean Genome and Epidemiology Study (KoGES; Ansan–Ansung and CITY/HEXA cohorts),
distributed by the National Biobank of Korea, Korea Disease Control and Prevention Agency.
The data cannot be shared by the authors; they are available on application to the National Biobank of Korea
(https://biobank.nih.go.kr).

The scripts expect the following files (names as distributed):

| Location | Files |
|---|---|
| `proj_dir` | `ASAS01to10.txt`, `CITY01to02.txt.gz` (phenotypes); `ASAS0400mx.txt.gz`, `CITY0822mx.txt.gz` (methylation M-values) |
| `chip_path` / `kch_dir` | K-chip and Exome chip genotypes as bigsnpr files (`*_v2.rds` with its `*_v2.bk`) |

The bigsnpr files can be created from the distributed PLINK files, for example:

```r
library(bigsnpr)
snp_readBed2("path/to/K-chip_folder/<K-chip>.bed", backingfile = "path/to/K-chip_folder/<K-chip>_v2")
```

The K-chip `.bk` file is about 32 GB. Keep it on a local disk; memory-mapping fails on cloud-synced drives.

## Setup

1. Set the path lines at the top of each script (`proj_dir`, `chip_path` in 01, `kch_dir` in 02 and 03).
2. Install the packages:
   - CRAN: `data.table`, `dplyr`, `bigsnpr`, `MatrixEQTL`, `httr`, `jsonlite`
   - Bioconductor: `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `GenomicFeatures`,
     `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `EpiDISH`
3. Run each script, e.g. `source("01_GWAS_36GenePanel.R", encoding = "UTF-8")`.

Versions used for the manuscript: R 4.5.2, bigsnpr 1.12.21, MatrixEQTL 2.4, org.Hs.eg.db 3.22.0,
TxDb.Hsapiens.UCSC.hg19.knownGene 3.22.1, IlluminaHumanMethylation450kanno.ilmn12.hg19 0.6.1, EpiDISH 2.26.0.

Each script writes a log file ending with `sessionInfo()`.
Lines marked `## 기대값:` ("expected values") give the values reported in the manuscript for comparison;
they are comments and are not used in any calculation.

## Outputs

| Folder (under `proj_dir`) | Script |
|---|---|
| `GWAS_results_36GenePanel/` | 01 |
| `mQTL_results/` | 02 |
| `validation_results/` | 03 |
| `Figures/` (PDF and 300 dpi PNG) | 05 |

All outputs are summary statistics. Approximate run time on a Windows 11 desktop: 01 several hours, 02 about 5 min, 03 about 11 min, 04 and 05 under 1 min.

## Notes

- Code comments are in Korean.
- The scripts use base R and standard package functions only (no user-defined functions).
