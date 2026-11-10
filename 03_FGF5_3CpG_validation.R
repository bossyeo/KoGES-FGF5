################################################################################################
## 03_FGF5_3CpG_validation.R
## KW25 논문: 메타분석 유의 CpG 3개 (cg12635937, cg10752158, cg12337525) 의 검증 / 민감도 분석
## 메틸화 파일 전체를 메모리에 올리지 않음 - gz 파일을 10만 줄씩 읽으며 필요한 행만 골라냄
##
## 하는 일 / 논문 대응
##   Part A. 450K manifest 의 probe SNP 확인                         -> Results "probe checks", Table 2 ‡
##   Part B. 필요한 CpG 행만 추출 (3개 CpG + EpiDISH 기준 CpG 333개)
##   Part C. 분석 자료 구성 + Table 1
##   Part D. 주 모형 재현, 흡연 제외 / complete-case 민감도, 부분 R2   -> Results, "Genotype explained 9.6% / 8.0%"
##   Part E. 유전형별 분포 그림                                       -> Supplementary Figure S1
##   Part F. 유전형-고혈압, CpG-고혈압                                -> Results, Supplementary Table S3
##   Part G. 범주형 공변량 + 혈구 구성(EpiDISH) + 유전 PC 5개, 메타분석 -> Supplementary Table S4
##
## 입력 (proj_dir): CITY0822mx.txt.gz, ASAS0400mx.txt.gz, CITY01to02.txt.gz, ASAS01to10.txt
##        (kch_dir): K-chip bigsnpr .rds/.bk (로컬 디스크)
## 출력 (proj_dir/validation_results/): 각 Part 의 csv, SuppFigS1 (pdf/png), 03_validation_log.txt
##
## 필요 패키지: data.table, dplyr, bigsnpr, EpiDISH, IlluminaHumanMethylation450kanno.ilmn12.hg19
## 개인 단위 자료는 저장하지 않음 - 요약 통계만 저장
## 사용자 정의 함수 없음 - base R 및 표준 패키지 함수만 사용
################################################################################################

## ---- 경로 (다른 컴퓨터에서는 이 두 줄만 바꾸면 됨) ---------------------------------------------------
proj_dir <- "path/to/project_folder"   # 표현형·메틸화 파일이 있고 결과가 저장될 폴더
t_start <- Sys.time()   # 전체 소요 시간 측정 시작
kch_dir  <- "path/to/K-chip_folder"   # K-chip bigsnpr .rds/.bk (약 30 GB - 로컬 디스크에 둘 것)

setwd(proj_dir)
out_dir <- file.path(proj_dir, "validation_results")
if (!dir.exists(out_dir)) dir.create(out_dir)
sink(file.path(out_dir, "03_validation_log.txt"), split = TRUE)
cat("실행 시각:", format(Sys.time()), "\n")

library(data.table); library(dplyr); library(bigsnpr)
library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data("centDHSbloodDMC.m", package = "EpiDISH")

## ---- 분석 대상: 02 결과에서 자동으로 읽음 -------------------------------------------------------------
## 마커 = 02 에서 mQTL 을 돌린 lead 마커 / CpG = 메타분석 Bonferroni 통과 CpG (P 순)
tab2 <- fread(file.path(proj_dir, "mQTL_results", "Table2_meta_Bonferroni_CpGs.csv"))
TARGET_MARKER <- tab2$SNP[1]
TARGET_CPG    <- tab2[order(p_M)]$gene
cat(">>> 02 결과에서 읽은 대상:", TARGET_MARKER, "|", paste(TARGET_CPG, collapse = ", "), "\n")
## 기대값: AX-83223998 | cg12635937, cg10752158, cg12337525

## ---- Part A. probe annotation ------------------------------------------------------------------------
## Probe_rs: probe 서열 안 SNP / CpG_rs: CpG 자리 SNP / SBE_rs: single-base extension 자리 SNP
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probe_tab <- as.data.frame(ann450k[TARGET_CPG, c("chr", "pos", "strand", "Type", "UCSC_RefGene_Name",
                                                 "Probe_rs", "Probe_maf", "CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf")])
print(probe_tab)
write.csv(probe_tab, file.path(out_dir, "PartA_probe_annotation.csv"))
## 기대값: cg10752158 만 Probe_rs = rs59133384 (MAF 0.13), 나머지 둘은 NA, RefGene 은 셋 다 비어 있음

## ---- 표현형 ------------------------------------------------------------------------------------------
CITY01to02 <- as.data.frame(fread("CITY01to02.txt.gz", header = TRUE, sep = "\t"))
f_asas <- "ASAS01to10.txt"
magic <- readBin(f_asas, "raw", n = 2)
con <- if (identical(as.integer(magic), c(31L, 139L))) gzfile(f_asas, "rb") else file(f_asas, "rb")
raw_asas <- readLines(con, warn = FALSE, encoding = "unknown"); close(con)
ASAS01to10 <- as.data.frame(fread(text = iconv(raw_asas, from = "CP949", to = "UTF-8", sub = "")))

## ---- Part B. 필요한 CpG 행만 추출 (3개 CpG + EpiDISH 기준 CpG) -------------------------------------
wanted   <- c(TARGET_CPG, rownames(centDHSbloodDMC.m))
met_file <- c(CITY = "CITY0822mx.txt.gz", ASAS = "ASAS0400mx.txt.gz")
M_list   <- list()
for (coh in c("CITY", "ASAS")) {
  con <- gzfile(met_file[coh], "r")
  hdr <- gsub('"', "", strsplit(readLines(con, n = 1), "\t")[[1]])
  hit <- character(0)
  repeat {
    x <- readLines(con, n = 100000)
    if (length(x) == 0) break
    hit <- c(hit, x[gsub('"', "", sub("\t.*", "", x)) %in% wanted])   # 첫 열(CpG ID)로 선택
  }
  close(con)
  sp <- do.call(rbind, strsplit(hit, "\t"))
  M  <- matrix(as.numeric(sp[, -1]), nrow = nrow(sp))
  rownames(M) <- gsub('"', "", sp[, 1])
  colnames(M) <- if (length(hdr) == ncol(sp) - 1) hdr else hdr[-1]
  stopifnot(all(TARGET_CPG %in% rownames(M)))
  cat(">>> [", coh, "] 추출:", nrow(M), "행 | EpiDISH 기준 CpG:",
      sum(rownames(M) %in% rownames(centDHSbloodDMC.m)), "/", nrow(centDHSbloodDMC.m), "\n")   # 논문 로그: 310 / 330
  M_list[[coh]] <- M
}

## ---- Part C. 분석 자료 (02 스크립트와 같은 공변량 처리 + 결측 표시) + Table 1 ------------------------
bk_path <- kch_dir
rds_candidates <- list.files(bk_path, pattern = "\\.rds$", full.names = TRUE)
v2 <- rds_candidates[grepl("_v2\\.rds$", rds_candidates)]
rds_file <- if (length(v2) > 0) v2[1] else rds_candidates[1]
KoGES_Kch <- snp_attach(rds_file)
INDEX_TARGET <- which(KoGES_Kch$map$marker.ID == TARGET_MARKER)

D_list <- list()
for (coh in c("CITY", "ASAS")) {
  if (coh == "CITY") {
    D <- dplyr::select(CITY01to02, DIST_ID, CT01_SEX, CT01_AGE, CT01_DRINK, CT01_SMOKE, CT01_HTN)
  } else {
    D <- dplyr::select(ASAS01to10, DIST_ID, AS01_SEX, AS01_AGE, AS01_DRINK, AS01_SMOKEA, AS01_PDHT)
  }
  colnames(D) <- c("DIST_ID", "SEX", "AGE", "DRKF", "SMKF", "HT")
  D$MISS_DRK <- is.na(D$DRKF); D$MISS_SMK <- is.na(D$SMKF)
  D$DRK_C <- factor(ifelse(D$MISS_DRK, "missing", D$DRKF))     # Part G 용 범주형 (결측 = 별도 범주)
  D$SMK_C <- factor(ifelse(D$MISS_SMK, "missing", D$SMKF))
  D$DRKF[D$MISS_DRK] <- 1; D$SMKF[D$MISS_SMK] <- 1             # 주 분석: 결측 -> 1
  D <- D[complete.cases(D[, c("SEX", "AGE", "DRKF", "SMKF")]), ]
  M <- M_list[[coh]]
  ID <- Reduce(intersect, list(D$DIST_ID, KoGES_Kch$fam$sample.ID, colnames(M)))
  D <- D[match(ID, D$DIST_ID), ]
  D$GENO <- KoGES_Kch$genotypes[match(ID, KoGES_Kch$fam$sample.ID), INDEX_TARGET]
  D <- cbind(D, t(M[TARGET_CPG, ID]))
  if (all(na.omit(D$HT) %in% c(1, 2))) D$HT <- D$HT - 1        # 1 = 없음, 2 = 있음 -> 0/1
  D_list[[coh]] <- D
  cat("\n>>> [", coh, "] n =", nrow(D), "| genotype 결측:", sum(is.na(D$GENO)), "\n")   # 360 / 280, 각 1명
  print(table(DRINK = D$DRK_C, useNA = "ifany")); print(table(SMOKE = D$SMK_C, useNA = "ifany"))
}

## Table 1
T1 <- data.frame()
for (coh in c("CITY", "ASAS")) {
  D <- D_list[[coh]]
  drk <- D$DRK_C[D$DRK_C != "missing"]
  T1 <- rbind(T1, data.frame(
    cohort = coh, n = nrow(D),
    men = sprintf("%d (%.1f)", sum(D$SEX == 1), 100 * mean(D$SEX == 1)),
    age = sprintf("%.1f (%.1f)", mean(D$AGE), sd(D$AGE)),
    drinking_1_2_3 = paste(sprintf("%d (%.1f)", as.vector(table(droplevels(drk))),
                                   100 * as.vector(prop.table(table(droplevels(drk))))), collapse = " / "),
    drinking_missing = sum(D$MISS_DRK),
    hypertension = sprintf("%d (%.1f)", sum(D$HT == 1, na.rm = TRUE), 100 * mean(D$HT == 1, na.rm = TRUE)),
    coded_allele_freq = round(mean(D$GENO, na.rm = TRUE) / 2, 2),
    smoking_missing = sum(D$MISS_SMK)))
}
print(T1)
write.csv(T1, file.path(out_dir, "Table1_participants.csv"), row.names = FALSE)
## 기대값: CITY 360, 남 261 (72.5), 51.1 (5.9), 93 (25.8) / 15 (4.2) / 252 (70.0), 고혈압 46 (12.8), 0.34, 흡연 결측 327
##         ASAS 280, 남 140 (50.0), 52.6 (8.2), 123 (44.4) / 21 (7.6) / 133 (48.0) [결측 3], 고혈압 43 (15.4), 0.36

## ---- Part D. 주 모형 재현 + 민감도 + 부분 R2 ----------------------------------------------------------
## lm 은 genotype 결측 1명을 제외하므로 n = 359 / 279 (MatrixEQTL 은 평균 대체로 360 / 280)
RES_D <- data.frame()
for (cpg in TARGET_CPG) {
  for (coh in c("CITY", "ASAS")) {
    D <- D_list[[coh]]
    D <- D[!is.na(D$GENO), ]
    D$Y <- D[[cpg]]
    f_main <- lm(Y ~ GENO + SEX + AGE + DRKF + SMKF, data = D)
    f_nosm <- lm(Y ~ GENO + SEX + AGE + DRKF, data = D)                                     # 흡연 제외
    f_cc   <- lm(Y ~ GENO + SEX + AGE + DRKF + SMKF, data = D[!D$MISS_DRK & !D$MISS_SMK, ])  # complete case
    f_red  <- lm(Y ~ SEX + AGE + DRKF + SMKF, data = D)                                     # 부분 R2 기준 모형
    c_main <- summary(f_main)$coefficients["GENO", ]
    c_nosm <- summary(f_nosm)$coefficients["GENO", ]
    c_cc   <- summary(f_cc)$coefficients["GENO", ]
    RES_D <- rbind(RES_D, data.frame(
      CpG = cpg, cohort = coh,
      n = nobs(f_main), beta = c_main[1], se = c_main[2], p = c_main[4],
      beta_noSMK = c_nosm[1], p_noSMK = c_nosm[4],
      n_complete = nobs(f_cc), beta_complete = c_cc[1], p_complete = c_cc[4],
      partial_R2 = (sum(resid(f_red)^2) - sum(resid(f_main)^2)) / sum(resid(f_red)^2)))
  }
}
rownames(RES_D) <- NULL
print(RES_D, digits = 3)
write.csv(RES_D, file.path(out_dir, "PartD_main_sensitivity_partialR2.csv"), row.names = FALSE)
## 기대값: cg12635937 CITY beta -0.120 (P 2.3e-9), 흡연 제외 -0.120 / ASAS complete case n 273, -0.103 (P 3.2e-6)
##         부분 R2: cg12635937 CITY 0.096 / ASAS 0.080, 나머지 두 CpG 0.041-0.070
##         (CITY 는 흡연 결측이 많아 complete case 가 33명뿐이라 해석하지 않음)

## ---- Part E. Supplementary Figure S1 (유전형별 분포) ---------------------------------------------------
for (dev_type in c("pdf", "png")) {
  if (dev_type == "pdf") pdf(file.path(out_dir, "SuppFigS1_3CpG_by_genotype.pdf"), width = 9, height = 6)
  if (dev_type == "png") png(file.path(out_dir, "SuppFigS1_3CpG_by_genotype.png"), width = 9, height = 6, units = "in", res = 300)
  par(mfrow = c(2, 3))
  for (coh in c("CITY", "ASAS")) {
    D <- D_list[[coh]]
    D <- D[!is.na(D$GENO), ]
    for (cpg in TARGET_CPG) {
      boxplot(D[[cpg]] ~ D$GENO, xlab = paste0(TARGET_MARKER, " genotype (0/1/2)"), ylab = "M-value",
              main = paste0(cpg, " [", coh, "]"), outline = FALSE, col = "grey90")
      stripchart(D[[cpg]] ~ D$GENO, vertical = TRUE, method = "jitter", pch = 16, cex = 0.4,
                 col = "#2a78d680", add = TRUE)
    }
  }
  dev.off()
}
for (coh in c("CITY", "ASAS")) {
  D <- D_list[[coh]]
  cat(coh, "cg12635937 유전형별 중앙값:", round(tapply(D$cg12635937, D$GENO, median), 2), "\n")
}
## 기대값: CITY 0.27 / 0.22 / 0.06

## ---- Part F. 고혈압 연관 (Supplementary Table S3) ------------------------------------------------------
RES_F <- data.frame()
for (coh in c("CITY", "ASAS")) {
  D <- D_list[[coh]]
  g0 <- summary(glm(HT ~ GENO + SEX + AGE + DRKF + SMKF, data = D, family = binomial))$coefficients["GENO", ]
  RES_F <- rbind(RES_F, data.frame(cohort = coh, exposure = "GENO", adj_GENO = NA,
                                   n = sum(!is.na(D$HT)), n_HT = sum(D$HT == 1, na.rm = TRUE),
                                   OR = exp(g0[1]), L95 = exp(g0[1] - 1.96 * g0[2]), U95 = exp(g0[1] + 1.96 * g0[2]), p = g0[4]))
  for (cpg in TARGET_CPG) {
    D$X <- D[[cpg]]
    g1 <- summary(glm(HT ~ X + SEX + AGE + DRKF + SMKF, data = D, family = binomial))$coefficients["X", ]
    g2 <- summary(glm(HT ~ X + GENO + SEX + AGE + DRKF + SMKF, data = D, family = binomial))$coefficients["X", ]
    RES_F <- rbind(RES_F,
      data.frame(cohort = coh, exposure = cpg, adj_GENO = FALSE, n = sum(!is.na(D$HT)), n_HT = sum(D$HT == 1, na.rm = TRUE),
                 OR = exp(g1[1]), L95 = exp(g1[1] - 1.96 * g1[2]), U95 = exp(g1[1] + 1.96 * g1[2]), p = g1[4]),
      data.frame(cohort = coh, exposure = cpg, adj_GENO = TRUE, n = sum(!is.na(D$HT)), n_HT = sum(D$HT == 1, na.rm = TRUE),
                 OR = exp(g2[1]), L95 = exp(g2[1] - 1.96 * g2[2]), U95 = exp(g2[1] + 1.96 * g2[2]), p = g2[4]))
  }
}
rownames(RES_F) <- NULL
print(RES_F, digits = 3)
write.csv(RES_F, file.path(out_dir, "SuppTableS3_hypertension.csv"), row.names = FALSE)
## 기대값: 유전형-고혈압 CITY OR 1.40 (0.89-2.20) P 0.15 / ASAS 1.73 (1.04-2.89) P 0.035
##         cg12635937 CITY 2.16 (0.65-7.14) P 0.21, 유전형 보정 3.38 (0.92-12.43) P 0.067 ... (Supp Table S3)

## ---- Part G. 범주형 공변량 + 혈구 구성 + 유전 PC (Supplementary Table S4) ------------------------------
## G-1. 혈구 구성: M-value -> beta = 2^M / (1 + 2^M), EpiDISH RPC, centDHSbloodDMC 기준
## G-2. 유전 PC 5개: 각 코호트 mQTL 표본에서, 상염색체 K-chip SNP 를 25개마다 1개 선택 후 MAF >= 0.05,
##      결측은 SNP 평균으로 대체, 표준화 후 prcomp
cell_cols <- c("Neutro", "CD4T", "CD8T", "B", "NK", "Mono")   # Eosino 는 제외 (합이 1 이라 하나는 빼야 함)
auto_idx  <- which(KoGES_Kch$map$chromosome %in% 1:22)
pca_idx   <- auto_idx[seq(1, length(auto_idx), by = 25)]

for (coh in c("CITY", "ASAS")) {
  D <- D_list[[coh]]
  M <- M_list[[coh]]
  ref_rows <- intersect(rownames(M), rownames(centDHSbloodDMC.m))
  B <- 2^M[ref_rows, D$DIST_ID] / (1 + 2^M[ref_rows, D$DIST_ID])
  cell <- epidish(beta.m = B, ref.m = centDHSbloodDMC.m, method = "RPC")$estF
  cat(">>> [", coh, "] 평균 혈구 비율:", paste(colnames(cell), round(colMeans(cell), 3), collapse = ", "), "\n")
  D <- cbind(D, cell[D$DIST_ID, cell_cols])

  Gs  <- KoGES_Kch$genotypes[match(D$DIST_ID, KoGES_Kch$fam$sample.ID), pca_idx]
  maf <- colMeans(Gs, na.rm = TRUE) / 2
  maf <- pmin(maf, 1 - maf)
  Gs  <- Gs[, !is.na(maf) & maf >= 0.05]
  na_pos <- which(is.na(Gs), arr.ind = TRUE)
  Gs[na_pos] <- colMeans(Gs, na.rm = TRUE)[na_pos[, 2]]
  cat(">>> [", coh, "] PCA 에 쓴 SNP 수:", ncol(Gs), "\n")       # 논문 로그: CITY 12,001 / ASAS 11,988
  pcs <- prcomp(Gs, center = TRUE, scale. = TRUE)$x[, 1:5]
  colnames(pcs) <- paste0("PC", 1:5)
  D_list[[coh]] <- cbind(D, pcs)
}

model_rhs <- c(
  numeric_as_main     = "GENO + SEX + AGE + DRKF + SMKF",
  categorical         = "GENO + SEX + AGE + DRK_C + SMK_C",
  categorical_noSMK   = "GENO + SEX + AGE + DRK_C",
  categorical_cell    = "GENO + SEX + AGE + DRK_C + SMK_C + Neutro + CD4T + CD8T + B + NK + Mono",
  categorical_cell_PC = "GENO + SEX + AGE + DRK_C + SMK_C + Neutro + CD4T + CD8T + B + NK + Mono + PC1 + PC2 + PC3 + PC4 + PC5")

RES_G <- data.frame()
for (cpg in TARGET_CPG) {
  for (mdl in names(model_rhs)) {
    for (coh in c("CITY", "ASAS")) {
      D <- D_list[[coh]]
      D <- D[!is.na(D$GENO), ]
      fit <- lm(as.formula(paste(cpg, "~", model_rhs[mdl])), data = D)
      cf  <- summary(fit)$coefficients["GENO", ]
      RES_G <- rbind(RES_G, data.frame(CpG = cpg, model = mdl, cohort = coh, n = nobs(fit),
                                       beta = cf[1], se = cf[2], p = cf[4]))
    }
    ## 두 코호트 고정효과 역분산가중 메타분석
    r <- RES_G[RES_G$CpG == cpg & RES_G$model == mdl, ]
    w <- 1 / r$se^2
    b_meta  <- sum(w * r$beta) / sum(w)
    se_meta <- 1 / sqrt(sum(w))
    RES_G <- rbind(RES_G, data.frame(CpG = cpg, model = mdl, cohort = "META", n = sum(r$n),
                                     beta = b_meta, se = se_meta, p = 2 * pnorm(-abs(b_meta / se_meta))))
  }
}
rownames(RES_G) <- NULL
print(RES_G, digits = 3)
write.csv(RES_G, file.path(out_dir, "SuppTableS4_sensitivity_cell_PC.csv"), row.names = FALSE)
## 기대값 (cg12635937, 분석 당시 FGF5_categorical_check.csv)
##   numeric_as_main     CITY -0.1201 (2.26e-9)  ASAS -0.1040 (1.76e-6)  META -0.1127 (5.08e-15)
##   categorical         CITY -0.1194 (2.45e-9)  ASAS -0.1073 (1.92e-6)  META -0.1140 (5.59e-15)
##   categorical_noSMK   CITY -0.1190 (2.97e-9)  ASAS -0.1053 (1.82e-6)  META -0.1128 (6.80e-15)
##   categorical_cell    CITY -0.1150 (8.05e-9)  ASAS -0.1010 (1.06e-5)  META -0.1090 (1.25e-13)
##   categorical_cell_PC CITY -0.1158 (9.07e-9)  ASAS -0.1011 (1.04e-5)  META -0.1095 (1.38e-13)
##   cg10752158 cell_PC META -0.1625 (1.29e-10) / cg12337525 cell_PC META -0.1111 (4.27e-8)

cat("\n=== 03 완료 ===\n")
cat("\n>>> 전체 소요 시간:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
sessionInfo()
sink()

