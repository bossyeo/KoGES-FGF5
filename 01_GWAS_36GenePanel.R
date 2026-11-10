################################################################################################
## 01_GWAS_36GenePanel.R
## KW25 논문: 후보 유전자 패널 연관분석 (Methods "Panel-restricted association analysis")
##
## 하는 일
##   36개 패널 항목 -> 실제 유전자 심볼로 확장(계열 접두어 포함, org.Hs.eg.db) -> 유전자 영역 +/-20 kb
##   칩 x 그룹(ASAS/CITY x 남/여) x 8개 만성질환 로지스틱 회귀:  질환 ~ GENO + AGE + ALC + SMK
##   P_FDR = 칩 x 그룹 x 질환 조합 안에서 BH 보정 / 전체 Bonferroni = 0.05 / 전체 검정 수
##
## 논문에서 나오는 값
##   Methods     : 282 유전자, K-chip 3,603 / Exome 580 마커, 전체 133,836 검정, Bonferroni 3.7e-7
##   Results     : 12개 FDR < 0.05 검정, AX-83223998 OR / 95% CI / P / FDR (4개 그룹), exm-rs16998073
##   Supp Table S1 / Supplementary Data 1 : 확장된 유전자 목록
##   Supp Table S2 : FDR < 0.05 검정 표
##   Figure 2      : 05_Figures.R 가 이 스크립트의 출력 파일을 읽어서 그림
##
## 입력: proj_dir 의 ASAS01to10.txt, CITY01to02.txt.gz / chip_path 의 칩별 bigsnpr .rds/.bk
##       (K-chip .bk 는 약 30 GB 라 구글 드라이브에서는 "Error when mapping file" 에러가 남
##        -> K-chip 폴더는 로컬 디스크에 두고 읽음. Exome 은 작아서 구글 드라이브에서도 됨)
## 출력: proj_dir/GWAS_results_36GenePanel/ (파일 목록은 맨 아래 Module 10 참고)
##
## 필요 패키지: data.table, bigsnpr, org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg19.knownGene, GenomicFeatures
## 사용자 정의 함수 없음 - base R 및 표준 패키지 함수만 사용
################################################################################################

## ---- Module 00: 경로 (다른 컴퓨터에서는 이 세 줄만 바꾸면 됨) -------------------------------------
proj_dir <- "path/to/project_folder"   # 표현형·메틸화 파일이 있고 결과가 저장될 폴더
chip_path <- c(Exo = "path/to/Exome_chip_folder",    # 칩별 bigsnpr .rds/.bk 폴더 (.rds 와 .bk 가 같이 있는 곳)
               Kch = "path/to/K-chip_folder")         # K-chip .bk 는 약 30 GB - 로컬 디스크에 둘 것

## 실행 전 확인: genotype 파일이 로컬에 있는지, .rds 와 .bk 가 짝으로 있는지
chip_list <- c("Exo", "Kch")      # 논문 분석에서 결과가 나온 칩 (Affymetrix 칩은 파일이 없어 분석 당시에도 스킵됨)
for (chip in chip_list) {
  rds_i <- list.files(chip_path[chip], pattern = "_v2\\.rds$", full.names = TRUE)
  if (length(rds_i) == 0) stop(chip, ": ", chip_path[chip], " 에 _v2.rds 가 없음")
  if (!file.exists(sub("\\.rds$", ".bk", rds_i[1]))) stop(chip, ": .rds 와 같은 폴더에 .bk 파일이 없음")
  cat(">>>", chip, "파일 확인:", rds_i[1], "\n")
}

path_out <- file.path(proj_dir, "GWAS_results_36GenePanel")
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)
sink(file.path(path_out, "01_GWAS_36GenePanel_log.txt"), split = TRUE)
cat("실행 시각:", format(Sys.time()), "\n")

## ---- Module 01: 패키지 -----------------------------------------------------------------------------
library(data.table)
library(bigsnpr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)

## ---- Module 02: 후보 유전자 패널 36개 항목 ------------------------------------------------------------
## 완전한 유전자 심볼이면 그대로, 아니면 그 접두어로 시작하는 유전자 전부로 확장 (예: FGF -> FGF1...FGF23, FGFR1...)
gene_panel_raw <- c(
  "TRPC", "KL", "ORAI", "STIM", "NPHS", "APOL", "P2RX", "PLCE", "RHOA", "ROCK",
  "LIMK", "CFL", "ACTN", "SYMPO", "CD2AP", "CACNA1C", "WT1", "FOXO1", "NFE2L2", "SOD",
  "NOX", "PPARGC1A", "CTCF", "FN1", "MMP", "TIMP", "ACTA", "VEGFA", "ANGPTL", "BCL",
  "HIF1A", "KIM", "LCN", "PAX", "FGF", "SIRT"
)
flank_bp <- 20000

all_symbol_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
token_report <- data.table()
token_gene   <- data.table()
for (token in gene_panel_raw) {
  if (token %in% all_symbol_keys) {
    genes_i <- token; mode_i <- "exact"
  } else {
    genes_i <- all_symbol_keys[grepl(paste0("^", token), all_symbol_keys)]; mode_i <- "prefix"
  }
  token_report <- rbind(token_report, data.table(TOKEN = token, MODE = mode_i, N_GENES = length(genes_i)))
  if (length(genes_i) > 0) token_gene <- rbind(token_gene, data.table(TOKEN = token, SYMBOL = genes_i))
}
gene_panel <- unique(token_gene$SYMBOL)
print(token_report)
cat(">>> 확장 후 유전자 심볼 수:", length(gene_panel), "(논문: 282, KIM / SYMPO 는 0개)\n")
cat(">>> org.Hs.eg.db 버전:", as.character(packageVersion("org.Hs.eg.db")), "(논문: 3.22.0)\n")

## ---- Module 03: 유전자 좌표 (hg19, UCSC knownGene) ---------------------------------------------------
sym2entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_panel, keytype = "SYMBOL", columns = "ENTREZID")
all_gene_ranges <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
match_idx <- match(sym2entrez$ENTREZID, names(all_gene_ranges))
keep <- which(!is.na(match_idx))
gene_coords <- data.table(
  SYMBOL = sym2entrez$SYMBOL[keep],
  CHR    = sub("^chr", "", as.character(seqnames(all_gene_ranges[match_idx[keep]]))),
  START  = start(all_gene_ranges[match_idx[keep]]) - flank_bp,
  END    = end(all_gene_ranges[match_idx[keep]]) + flank_bp
)
cat(">>> 좌표가 있는 유전자 수:", nrow(gene_coords), "/", length(gene_panel), "\n")

## Supplementary Data 1: 항목별 유전자 목록 + 좌표
supp_data1 <- merge(token_gene, gene_coords, by = "SYMBOL", all.x = TRUE, sort = FALSE)
setcolorder(supp_data1, c("TOKEN", "SYMBOL", "CHR", "START", "END"))
fwrite(supp_data1, file.path(path_out, "SuppData1_panel_gene_list.csv"))
fwrite(token_report, file.path(path_out, "SuppTableS1_panel_token_report.csv"))
setkey(gene_coords, CHR, START, END)

## ---- Module 04: 표현형 (CP949 -> UTF-8, gzip 자동 감지) ------------------------------------------
## fread 에 경로를 직접 주면 한글+공백 경로에서 압축 해제가 실패할 수 있어 R 커넥션으로 읽음
f_asas <- file.path(proj_dir, "ASAS01to10.txt")
magic <- readBin(f_asas, "raw", n = 2)
con <- if (identical(as.integer(magic), c(31L, 139L))) gzfile(f_asas, "rb") else file(f_asas, "rb")
raw_asas <- readLines(con, warn = FALSE, encoding = "unknown"); close(con)
pheno_ASAS <- fread(text = iconv(raw_asas, from = "CP949", to = "UTF-8", sub = ""))

f_city <- file.path(proj_dir, "CITY01to02.txt.gz")
magic <- readBin(f_city, "raw", n = 2)
con <- if (identical(as.integer(magic), c(31L, 139L))) gzfile(f_city, "rb") else file(f_city, "rb")
raw_city <- readLines(con, warn = FALSE, encoding = "unknown"); close(con)
pheno_CITY <- fread(text = iconv(raw_city, from = "CP949", to = "UTF-8", sub = ""))

## ---- Module 05: 실행 목록 (칩 목록은 Module 00 에 있음) ----------------------------------------------
cohort_list   <- c("ASAS_M", "ASAS_F", "CITY_M", "CITY_F")
disease_codes <- c("PDTH", "PDHT", "PDDM", "PDLP", "PDUL", "PDAS", "PDAL", "PDKD")
city_disease_map <- c(PDTH = "THY", PDHT = "HTN", PDDM = "DM", PDLP = "LIP",      # CITY 질환 열 이름
                      PDUL = "GASULCER", PDAS = "ASTH", PDAL = "ALLER", PDKD = "ANEPH")

gwas_results_all <- list()
panel_marker_n   <- data.table()

## ---- Module 06~08: 칩 x 그룹 x 질환 ----------------------------------------------------------------
for (chip in chip_list) {

  cat("\n=== 칩:", chip, "===\n")
  rds_file <- list.files(chip_path[chip], pattern = "_v2\\.rds$", full.names = TRUE)[1]   # 재생성한 _v2 파일
  cat(">>> 사용 파일:", rds_file, "\n")

  obj_bigsnp <- snp_attach(rds_file)
  G <- obj_bigsnp$genotypes; map <- obj_bigsnp$map; fam <- obj_bigsnp$fam

  map_dt <- data.table(IDX = seq_len(nrow(map)), CHR = as.character(map$chromosome),
                       START = map$physical.pos, END = map$physical.pos)
  hits <- foverlaps(map_dt, gene_coords, by.x = c("CHR", "START", "END"), type = "within", nomatch = 0L)
  panel_idx <- sort(unique(hits$IDX))
  cat(">>>", chip, "패널 영역 마커 수:", length(panel_idx), "/", nrow(map), "(논문: Kch 3,603 / Exo 580)\n")
  if (length(panel_idx) == 0) next
  panel_marker_n <- rbind(panel_marker_n, data.table(CHIP = chip, N_PANEL_MARKERS = length(panel_idx)))

  for (cohort in cohort_list) {

    cohort_base <- sub("_[MF]$", "", cohort)
    sex_code    <- ifelse(grepl("_M$", cohort), 1, 2)
    pheno       <- if (cohort_base == "ASAS") pheno_ASAS else pheno_CITY

    round_prefix <- if (cohort_base == "ASAS") "AS01_" else "CT01_"
    id_col <- if ("DIST_ID" %in% names(pheno)) "DIST_ID" else "ID"
    common_id <- intersect(fam$sample.ID, pheno[[id_col]])
    if (length(common_id) < 30) { cat(">>>", chip, cohort, ": 공유 표본", length(common_id), "명, 스킵\n"); next }

    idx_geno  <- match(common_id, fam$sample.ID)
    pheno_sub <- pheno[match(common_id, pheno[[id_col]]), ]
    sex_col   <- paste0(round_prefix, "SEX")
    keep_sex  <- which(pheno_sub[[sex_col]] == sex_code)
    if (length(keep_sex) < 30) { cat(">>>", chip, cohort, ": 성별 필터 후 표본 부족, 스킵\n"); next }
    idx_geno_f  <- idx_geno[keep_sex]
    pheno_sub_f <- pheno_sub[keep_sex, ]

    ## 공변량 열: 이름 패턴에 처음 걸리는 열 (분석 당시 방식 그대로). 실제로 쓰인 열을 로그에 기록
    age_col <- paste0(round_prefix, "AGE")
    alc_col <- grep(paste0("^", round_prefix, ".*(ALC|DRINK)"), names(pheno_sub_f), value = TRUE)[1]
    smk_col <- grep(paste0("^", round_prefix, ".*(SMK|SMOKE)"), names(pheno_sub_f), value = TRUE)[1]
    cat(">>>", chip, cohort, "공변량 열: AGE =", age_col, "| ALC =", alc_col, "| SMK =", smk_col, "\n")
    covar <- as.matrix(data.frame(AGE = pheno_sub_f[[age_col]], ALC = pheno_sub_f[[alc_col]],
                                  SMK = pheno_sub_f[[smk_col]]))

    for (dcode in disease_codes) {

      dcol <- paste0(round_prefix, if (cohort_base == "CITY") city_disease_map[dcode] else dcode)
      if (!(dcol %in% names(pheno_sub_f))) { cat(">>>", chip, cohort, dcode, ":", dcol, "없음, 스킵\n"); next }
      y <- as.numeric(pheno_sub_f[[dcol]])
      if (all(na.omit(unique(y)) %in% c(1, 2))) y <- y - 1           # 1 = 없음, 2 = 있음 -> 0/1
      valid <- which(!is.na(y) & y %in% c(0, 1))
      if (length(valid) < 30) next

      sub_mat  <- G[idx_geno_f[valid], panel_idx]                      # 읽기 전용
      res_list <- vector("list", length(panel_idx))

      for (cc in seq_along(panel_idx)) {
        col_vals <- sub_mat[, cc]
        if (anyNA(col_vals)) {                                         # genotype 결측 -> 최빈값
          tab <- table(col_vals)
          col_vals[is.na(col_vals)] <- as.numeric(names(tab)[which.max(tab)])
        }
        if (length(unique(col_vals)) < 2) next                         # 단형 마커 제외
        df <- data.frame(Y = y[valid], GENO = col_vals, AGE = covar[valid, "AGE"],
                         ALC = covar[valid, "ALC"], SMK = covar[valid, "SMK"])
        fit <- suppressWarnings(glm(Y ~ GENO + AGE + ALC + SMK, data = df, family = binomial()))  # 공변량 결측 제외
        if (!fit$converged) next
        ct <- summary(fit)$coefficients
        if (!("GENO" %in% rownames(ct))) next
        mk <- panel_idx[cc]
        res_list[[cc]] <- data.table(
          MARKER = map$marker.ID[mk], CHR = map$chromosome[mk], POS = map$physical.pos[mk],
          ALLELE1 = map$allele1[mk], ALLELE2 = map$allele2[mk],
          BETA = ct["GENO", "Estimate"], SE = ct["GENO", "Std. Error"], P = ct["GENO", "Pr(>|z|)"],
          N = nobs(fit), N_CASE = sum(fit$y == 1))
      }

      res <- rbindlist(res_list)
      if (nrow(res) == 0) next
      res[, P_BONF := pmin(P * .N, 1)]                                 # 조합 안 Bonferroni
      res[, P_FDR  := p.adjust(P, method = "BH")]                      # 조합 안 BH-FDR
      res[, `:=`(CHIP = chip, COHORT = cohort, DISEASE = dcode)]
      setcolorder(res, c("CHIP", "COHORT", "DISEASE"))
      fwrite(res, file.path(path_out, paste0("GWAS_36GenePanel_", chip, "_", cohort, "_", dcode, ".txt")), sep = "\t")
      gwas_results_all[[paste(chip, cohort, dcode, sep = "_")]] <- res
    }
  }
}

## ---- Module 09: 전체 결과 통합 --------------------------------------------------------------------
gwas_combined <- rbindlist(gwas_results_all, use.names = TRUE, fill = TRUE)
setorder(gwas_combined, P)
gwas_combined[, COHORT_BASE := sub("_[MF]$", "", COHORT)]
gwas_combined[, SEX := sub(".*_", "", COHORT)]
fwrite(gwas_combined, file.path(path_out, "GWAS_36GenePanel_ALL_combined.txt"), sep = "\t")

## ---- Module 10: 논문 값 ------------------------------------------------------------------------------
n_tests  <- nrow(gwas_combined)
bonf_thr <- 0.05 / n_tests
cat("\n>>> 전체 검정 수:", n_tests, "| Bonferroni 기준:", signif(bonf_thr, 2), "(논문: 133,836 / 3.7e-7)\n")
print(gwas_combined[, .N, by = CHIP])                                  # 논문: Kch 115,296 / Exo 18,540
fwrite(merge(panel_marker_n, gwas_combined[, .(N_TESTS = .N), by = CHIP], by = "CHIP"),
       file.path(path_out, "panel_marker_and_test_counts.csv"))

## (1) Supp Table S2: FDR < 0.05 검정 (Figure 2a 입력)
s2 <- gwas_combined[P_FDR < 0.05]
s2[, `:=`(OR = exp(BETA), OR_L95 = exp(BETA - 1.96 * SE), OR_U95 = exp(BETA + 1.96 * SE),
          BONFERRONI = P < bonf_thr)]
setorder(s2, P_FDR)
print(s2[, .(MARKER, DISEASE, COHORT, OR, OR_L95, OR_U95, P, P_FDR, BONFERRONI)], digits = 3)
fwrite(s2, file.path(path_out, "SuppTableS2_FDRsignif.csv"))
## 기대값: 12개 (AX-83223998 CITY_M OR 1.34 (1.23-1.47), P 5.8e-11, FDR 2.1e-7 ... exm1692187 CITY_M)

## (2) rs16998073 (AX-83223998 / exm-rs16998073) - 고혈압, 4개 그룹 (본문, Figure 2b 입력)
fgf5 <- gwas_combined[MARKER %in% c("AX-83223998", "exm-rs16998073") & DISEASE == "PDHT"]
fgf5[, `:=`(OR = exp(BETA), OR_L95 = exp(BETA - 1.96 * SE), OR_U95 = exp(BETA + 1.96 * SE))]
setorder(fgf5, MARKER, COHORT)
print(fgf5[, .(CHIP, MARKER, ALLELE1, COHORT, N, N_CASE, OR, OR_L95, OR_U95, P, P_FDR)], digits = 3)
fwrite(fgf5, file.path(path_out, "FGF5_rs16998073_PDHT_byGroup.csv"))
## 기대값 (AX-83223998): CITY_M 1.34 (1.23-1.47) P 5.8e-11 / CITY_F 1.20 (1.12-1.30) P 1.7e-6
##                      ASAS_M 1.30 (1.10-1.54) P 1.9e-3 FDR 0.77 / ASAS_F 1.16 (0.99-1.36) P 0.062
##         (exm-rs16998073): CITY_F 1.52 (1.27-1.83) FDR 0.0039

cat("\n=== 01 완료 ===\n")
sessionInfo()
sink()



