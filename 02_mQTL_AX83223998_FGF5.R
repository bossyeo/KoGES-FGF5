################################################################################################
## 02_mQTL_AX83223998_FGF5.R
## KW25 논문: rs16998073 (K-chip AX-83223998) 전장 mQTL + 두 코호트 메타분석
## (Methods "DNA methylation and mQTL analysis", "Discovery, replication and meta-analysis")
##
## 하는 일
##   1. MatrixEQTL 선형모형  CpG M-value ~ GENO + SEX + AGE + DRKF + SMKF
##      (공변량은 숫자 코드 그대로, 음주/흡연 결측은 1 로 채움; genotype 결측은 MatrixEQTL 이 평균으로 대체)
##      CITY (discovery, EPIC) / ASAS (replication, 450K) 각각 실행
##   2. 두 코호트 공통 CpG 에 대해 고정효과 역분산가중 메타분석 (SE = beta / t), Cochran's Q
##
## 논문에서 나오는 값
##   Results : CITY cg12635937 t -6.14, beta -0.120, P 2.2e-9, FDR 0.0017 / 2위 cg10752158 P 4.2e-7, FDR 0.17
##             ASAS cg12635937 1위, t -4.89, beta -0.104, P 1.7e-6 / cg10752158 3위
##             메타: 공통 413,632 CpG, Bonferroni 1.2e-7 통과 3개
##   Table 2 : 3개 CpG 의 위치, 거리, CITY/ASAS beta(P), meta beta(SE), P, P_het
##   Figure 3: 05_Figures.R 가 이 스크립트의 출력 파일을 읽어서 그림
##
## 입력 (proj_dir): CITY0822mx.txt.gz, ASAS0400mx.txt.gz, CITY01to02.txt.gz, ASAS01to10.txt
##        (kch_dir): K-chip bigsnpr .rds/.bk (로컬 디스크)
## 출력 (proj_dir/mQTL_results/)
##   meQTL_AX83223998_FGF5_CITY_sorted.txt / _ASAS_sorted.txt   전체 CpG 결과 (P 순)
##   mQTL_meta_CITY_ASAS_all.txt                             공통 CpG 전체 메타분석 결과
##   Table2_meta_Bonferroni_CpGs.csv                            Table 2
##   02_mQTL_log.txt                                            콘솔 기록 + 패키지 버전
##
## 필요 패키지: data.table, dplyr, bigsnpr, MatrixEQTL, IlluminaHumanMethylation450kanno.ilmn12.hg19
## 메모리: 메틸화 파일 전체를 읽으므로 RAM 이 충분해야 함 (CITY 를 끝내고 지운 뒤 ASAS 를 읽음)
## 사용자 정의 함수 없음 - base R 및 표준 패키지 함수만 사용
################################################################################################

## ---- Module 00: 경로 (다른 컴퓨터에서는 이 두 줄만 바꾸면 됨) -------------------------------------
proj_dir <- "path/to/project_folder"   # 표현형·메틸화 파일이 있고 결과가 저장될 폴더
t_start <- Sys.time()   # 전체 소요 시간 측정 시작
kch_dir  <- "path/to/K-chip_folder"   # K-chip bigsnpr .rds/.bk (약 30 GB - 로컬 디스크에 둘 것)

setwd(proj_dir)
out_dir <- file.path(proj_dir, "mQTL_results")
if (!dir.exists(out_dir)) dir.create(out_dir)
sink(file.path(out_dir, "02_mQTL_log.txt"), split = TRUE)
cat("실행 시각:", format(Sys.time()), "\n")

library(data.table); library(dplyr)
library(bigsnpr)
library(MatrixEQTL)

## ---- mQTL 로 넘길 마커: 01 결과에서 자동 선택 ------------------------------------------------------
## 규칙 (Methods "the lead marker"): 패널 연관분석에서 FDR < 0.05 인 검정 중 P 가 가장 작은 마커
## 메틸화 자료가 있는 사람은 K-chip 으로 유전형이 측정되어 있으므로, 선택된 마커가 K-chip 마커인지 확인
s2_gwas <- fread(file.path(proj_dir, "GWAS_results_36GenePanel", "SuppTableS2_FDRsignif.csv"))
lead    <- s2_gwas[order(P)][1]
TARGET_MARKER <- lead$MARKER
cat(">>> 01 결과에서 선택된 lead 마커:", TARGET_MARKER, "|", lead$DISEASE, "|", lead$COHORT, "| P =", signif(lead$P, 3),
    "| FDR =", signif(lead$P_FDR, 3), "\n")                     # 기대값: AX-83223998, PDHT, CITY_M, P 5.8e-11
if (lead$CHIP != "Kch") stop("lead 마커가 K-chip 마커가 아님 - mQTL 표본과 유전형 칩이 맞지 않음")

## ---- Module 01: K-chip genotype -------------------------------------------------------------------
bk_path <- kch_dir
rds_candidates <- list.files(bk_path, pattern = "\\.rds$", full.names = TRUE)
v2 <- rds_candidates[grepl("_v2\\.rds$", rds_candidates)]
rds_file <- if (length(v2) > 0) v2[1] else rds_candidates[1]
cat(">>> 사용 파일:", rds_file, "\n")
KoGES_Kch <- snp_attach(rds_file)
INDEX_TARGET <- which(KoGES_Kch$map$marker.ID == TARGET_MARKER)
if (length(INDEX_TARGET) == 0) stop("타겟 마커를 K-chip map 에서 찾을 수 없음")
print(KoGES_Kch$map[INDEX_TARGET, ])   # chr4:81,184,341 / allele1 = T (genotype 0/1/2 = T 개수)

## ---- Module 02: 공변량 ---------------------------------------------------------------------------------
CITY01to02 <- as.data.frame(fread("CITY01to02.txt.gz", header = TRUE, sep = "\t"))
CITY_COV <- dplyr::select(CITY01to02, DIST_ID, CT01_SEX, CT01_AGE, CT01_DRINK, CT01_SMOKE)
colnames(CITY_COV) <- c("DIST_ID", "SEX", "AGE", "DRKF", "SMKF")
CITY_COV[, c("DRKF", "SMKF")][is.na(CITY_COV[, c("DRKF", "SMKF")])] <- 1
CITY_COV <- CITY_COV[complete.cases(CITY_COV), ]

f_asas <- "ASAS01to10.txt"
magic <- readBin(f_asas, "raw", n = 2)
con <- if (identical(as.integer(magic), c(31L, 139L))) gzfile(f_asas, "rb") else file(f_asas, "rb")
raw_asas <- readLines(con, warn = FALSE, encoding = "unknown"); close(con)
ASAS01to10 <- as.data.frame(fread(text = iconv(raw_asas, from = "CP949", to = "UTF-8", sub = "")))
ASAS_COV <- dplyr::select(ASAS01to10, DIST_ID, AS01_SEX, AS01_AGE, AS01_DRINK, AS01_SMOKEA)
colnames(ASAS_COV) <- c("DIST_ID", "SEX", "AGE", "DRKF", "SMKF")
ASAS_COV[, c("DRKF", "SMKF")][is.na(ASAS_COV[, c("DRKF", "SMKF")])] <- 1
ASAS_COV <- ASAS_COV[complete.cases(ASAS_COV), ]

## ---- Module 03: MatrixEQTL - CITY 와 ASAS 를 같은 코드로 차례대로 실행 ---------------------------------
met_file <- c(CITY = "CITY0822mx.txt.gz", ASAS = "ASAS0400mx.txt.gz")
res_all  <- list()

for (coh in c("CITY", "ASAS")) {

  COV <- if (coh == "CITY") CITY_COV else ASAS_COV
  Met <- as.data.frame(fread(met_file[coh], header = TRUE, sep = "\t"))
  rownames(Met) <- Met$V1
  Met <- Met[, -1]
  cat("\n>>> [", coh, "] 메틸화:", nrow(Met), "CpG x", ncol(Met), "명\n")      # 논문: CITY 803,393 / ASAS 451,245

  ID <- Reduce(intersect, list(COV$DIST_ID, KoGES_Kch$fam$sample.ID, colnames(Met)))
  cat(">>> [", coh, "] 분석 인원:", length(ID), "\n")                           # 논문: CITY 360 / ASAS 280

  snp_matrix <- matrix(KoGES_Kch$genotypes[match(ID, KoGES_Kch$fam$sample.ID), INDEX_TARGET], nrow = 1,
                       dimnames = list(TARGET_MARKER, ID))
  cat(">>> [", coh, "] genotype 결측:", sum(is.na(snp_matrix)), "명 (MatrixEQTL 이 평균으로 대체)\n")
  meth_matrix <- as.matrix(Met[, ID])
  cov_matrix  <- t(as.matrix(COV[match(ID, COV$DIST_ID), c("SEX", "AGE", "DRKF", "SMKF")]))
  colnames(cov_matrix) <- ID
  rm(Met); invisible(gc())

  sd_snp <- SlicedData$new(); sd_snp$CreateFromMatrix(snp_matrix);  sd_snp$ResliceCombined(sliceSize = 2000)
  sd_met <- SlicedData$new(); sd_met$CreateFromMatrix(meth_matrix); sd_met$ResliceCombined(sliceSize = 2000)
  sd_cov <- SlicedData$new(); sd_cov$CreateFromMatrix(cov_matrix);  sd_cov$ResliceCombined(sliceSize = 2000)

  me <- Matrix_eQTL_engine(snps = sd_snp, gene = sd_met, cvrt = sd_cov,
                           output_file_name = file.path(out_dir, paste0("meQTL_AX83223998_FGF5_", coh, "_result.txt")),
                           pvOutputThreshold = 1, useModel = modelLINEAR, verbose = TRUE, pvalue.hist = FALSE)
  cat(">>> [", coh, "] 잔차 자유도:", me$param$dfFull, "(n - 6 이어야 함)\n")

  res <- as.data.table(me$all$eqtls)[order(pvalue)]
  res[, `:=`(snps = as.character(snps), gene = as.character(gene))]   # factor -> 문자 (코호트 간 병합용)
  fwrite(res, file.path(out_dir, paste0("meQTL_AX83223998_FGF5_", coh, "_sorted.txt")), sep = "\t")
  cat(">>> [", coh, "] FDR < 0.05 CpG 수:", sum(res$FDR < 0.05), "\n")
  print(head(res, 5))
  res_all[[coh]] <- res
  rm(meth_matrix, sd_met, me); invisible(gc())
}
## 기대값: CITY 1위 cg12635937 t -6.14, beta -0.120, P 2.2e-9, FDR 0.0017 (FDR < 0.05 는 1개)
##         ASAS 1위 cg12635937 t -4.89, beta -0.104, P 1.7e-6 / cg10752158 3위

## ---- Module 04: 고정효과 역분산가중 메타분석 (공통 CpG) -------------------------------------------------
m <- merge(res_all$CITY[, .(gene, beta_C = beta, t_C = statistic, p_C = pvalue)],
           res_all$ASAS[, .(gene, beta_A = beta, t_A = statistic, p_A = pvalue)], by = "gene")
m[, `:=`(se_C = beta_C / t_C, se_A = beta_A / t_A)]
m[, `:=`(w_C = 1 / se_C^2, w_A = 1 / se_A^2)]
m[, beta_M := (w_C * beta_C + w_A * beta_A) / (w_C + w_A)]
m[, se_M   := 1 / sqrt(w_C + w_A)]
m[, p_M    := 2 * pnorm(-abs(beta_M / se_M))]
m[, Q      := w_C * (beta_C - beta_M)^2 + w_A * (beta_A - beta_M)^2]
m[, p_het  := pchisq(Q, df = 1, lower.tail = FALSE)]
m[, FDR_M  := p.adjust(p_M, method = "BH")]
setorder(m, p_M)
bonf_meta <- 0.05 / nrow(m)
cat("\n>>> 공통 CpG 수:", nrow(m), "| Bonferroni:", signif(bonf_meta, 2), "| 통과:", sum(m$p_M < bonf_meta), "\n")
fwrite(m, file.path(out_dir, "mQTL_meta_CITY_ASAS_all.txt"), sep = "\t")
## 기대값: 413,632 CpG, Bonferroni 1.2e-7, 통과 3개 (cg12635937, cg10752158, cg12337525)

## ---- Module 05: Table 2 (Bonferroni 통과 CpG + 450K 위치) ---------------------------------------------
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
tab2 <- m[p_M < bonf_meta]
tab2[, `:=`(chr = ann450k[gene, "chr"], pos = ann450k[gene, "pos"],
            RefGene = ann450k[gene, "UCSC_RefGene_Name"], Probe_rs = ann450k[gene, "Probe_rs"])]
tab2[, `:=`(SNP = TARGET_MARKER, SNP_POS = KoGES_Kch$map$physical.pos[INDEX_TARGET])]   # 03, 05 가 읽어 감
tab2[, dist_from_SNP := pos - KoGES_Kch$map$physical.pos[INDEX_TARGET]]
setorder(tab2, pos)
print(tab2[, .(gene, pos, dist_from_SNP, beta_C, p_C, beta_A, p_A, beta_M, se_M, p_M, p_het, Probe_rs)], digits = 3)
fwrite(tab2, file.path(out_dir, "Table2_meta_Bonferroni_CpGs.csv"))
## 기대값 (Table 2)
##   cg12635937 81,184,925  584 bp  -0.120 (2.2e-9)  -0.104 (1.7e-6)  -0.113 (0.014)  4.8e-15  P_het 0.58
##   cg10752158 81,185,022  681 bp  -0.170 (4.2e-7)  -0.162 (1.3e-5)  -0.167 (0.024)  1.0e-11  P_het 0.88
##   cg12337525 81,185,228  887 bp  -0.099 (1.2e-4)  -0.143 (1.4e-5)  -0.116 (0.020)  7.0e-9   P_het 0.28

cat("\n=== 02 완료 ===\n")
cat("\n>>> 전체 소요 시간:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
sessionInfo()
sink()

