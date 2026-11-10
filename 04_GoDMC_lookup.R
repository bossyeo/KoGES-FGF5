################################################################################################
## 04_GoDMC_lookup.R
## KW25 논문: 3개 CpG 를 GoDMC 혈액 mQTL 데이터베이스에서 조회 (Discussion 의 GoDMC 비교)
##
## 하는 일
##   GoDMC 조회 API (http://api.godmc.org.uk/v0.1/query) 에 3개 CpG 를 POST 로 보내고
##   (모든 연관, clumping 없음) 결과를 저장한 뒤 rs16998073 행을 뽑아 우리 결과와 대립유전자 방향을 맞춤
##
## 논문에서 나오는 값 (Discussion)
##   rs16998073 - cg12635937: beta -0.212, P 3.7e-103, n 24,467 (allele1 = T, T 빈도 0.30)
##   cg12337525: P 2.5e-54 / cg10752158: P 5.9e-44
##   각 CpG 의 가장 강한 연관 SNP (예: rs12054632, rs1379705)
##
## 참고
##   - 인터넷 연결이 필요함. 분석 당시 조회일: 2026-10-13
##   - GoDMC 서버는 GET 요청에 "hello" 만 돌려주므로 POST 를 써야 함
##   - 서버가 응답하지 않으면, 분석 당시 저장한 GoDMC_FGF5_CpG_lookup.csv 를 읽어서 Module 02 부터 실행
##
## 필요 패키지: httr, jsonlite, data.table
## 사용자 정의 함수 없음 - base R 및 표준 패키지 함수만 사용
################################################################################################

proj_dir <- "path/to/project_folder"   # 표현형·메틸화 파일이 있고 결과가 저장될 폴더
t_start <- Sys.time()   # 전체 소요 시간 측정 시작
setwd(proj_dir)
library(httr); library(jsonlite); library(data.table)

TARGET_CPG <- c("cg12635937", "cg10752158", "cg12337525")
out_file   <- "GoDMC_FGF5_CpG_lookup.csv"

## ---- Module 01: API 조회 (모든 연관, clumped = 0) --------------------------------------------------
resp <- try(POST("http://api.godmc.org.uk/v0.1/query",
                 body = list(cpgs = TARGET_CPG, clumped = 0), encode = "json", timeout(120)), silent = TRUE)

if (!inherits(resp, "try-error") && status_code(resp) == 200) {
  godmc <- as.data.table(fromJSON(content(resp, as = "text", encoding = "UTF-8")))
  fwrite(godmc, out_file)
  cat(">>> GoDMC 조회 완료:", nrow(godmc), "행 저장 ->", out_file, "\n")
} else {
  cat(">>> GoDMC 서버 응답 없음 - 분석 당시 저장한", out_file, "을 사용\n")
  godmc <- fread(out_file)
}

## ---- Module 02: rs16998073 결과와 각 CpG 의 최상위 SNP ----------------------------------------------
print(godmc[rsid == "rs16998073", .(cpg, rsid, a1, a2, freq_a1, beta_a1, se, pval, samplesize, hetpval)])
## 기대값: cg12635937 a1 = T, freq_a1 0.2967, beta_a1 -0.212, pval 3.71e-103, samplesize 24,467
##         cg12337525 pval 2.47e-54 / cg10752158 pval 5.85e-44
## 우리 자료도 T 를 센 genotype (AX-83223998 allele1 = T, 빈도 0.34-0.36) 이므로 같은 방향 (beta < 0)

top_snp <- godmc[order(pval)][, .SD[1], by = cpg]
print(top_snp[, .(cpg, rsid, a1, freq_a1, beta_a1, pval, samplesize)])
## 기대값: 최상위 SNP 는 rs16998073 이 아니라 인근 SNP
##         cg12635937 rs12054632 (빈도 0.10, P 3.8e-303) / cg12337525 rs1379705 (0.10, P 2.4e-210)
##         cg10752158 rs11099098 (0.30, P 5.4e-47)

cat("\n>>> 전체 소요 시간:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
