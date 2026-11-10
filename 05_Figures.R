################################################################################################
## 05_Figures.R
## KW25 논문: Figure 1-3 (base R 그래픽). Supplementary Figure S1 은 03 스크립트에서 만들어짐
##
## 입력 (01, 02 스크립트의 출력 - 개인 단위 자료 없음, 요약 결과만 사용)
##   GWAS_results_36GenePanel/SuppTableS2_FDRsignif.csv            Figure 1, 2a
##   GWAS_results_36GenePanel/FGF5_rs16998073_PDHT_byGroup.csv     Figure 1, 2b
##   GWAS_results_36GenePanel/panel_marker_and_test_counts.csv     Figure 1
##   mQTL_results/meQTL_AX83223998_FGF5_CITY_sorted.txt            Figure 1, 3a
##   mQTL_results/meQTL_AX83223998_FGF5_ASAS_sorted.txt            Figure 1
##   mQTL_results/mQTL_meta_CITY_ASAS_all.txt                   Figure 1
##   mQTL_results/Table2_meta_Bonferroni_CpGs.csv                  Figure 1, 3a, 3b, 3c
## 출력 (proj_dir/Figures/): Fig1_workflow, Fig2_GWAS_panel, Fig3_mQTL_FGF5_locus (.pdf + 300 dpi .png)
##
## 필요 패키지: data.table (그림은 base R graphics, cairo 장치)
## 사용자 정의 함수 없음 - base R 및 표준 패키지 함수만 사용
################################################################################################

proj_dir <- "path/to/project_folder"   # 표현형·메틸화 파일이 있고 결과가 저장될 폴더
t_start <- Sys.time()   # 전체 소요 시간 측정 시작
gwas_dir <- file.path(proj_dir, "GWAS_results_36GenePanel")
mqtl_dir <- file.path(proj_dir, "mQTL_results")
fig_dir  <- file.path(proj_dir, "Figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir)
library(data.table)

## ---- 공통 설정 ----------------------------------------------------------------------------------------
FONT <- "Arial"
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID <- "#e1e0d9"; AXIS <- "#c3c2b7"
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; BLUE_T <- "#e8f0fb"; PANEL <- "#f4f3ef"
MINUS <- "\u2212"
FGF5_TSS  <- 81187742      # FGF5 전사 시작 (RefSeq NM_004464, + strand, hg19)
disease_name <- c(PDHT = "Hypertension", PDDM = "Diabetes mellitus", PDAL = "Allergic disease",
                  PDLP = "Dyslipidemia", PDKD = "Kidney disease", PDAS = "Asthma",
                  PDUL = "Gastric ulcer", PDTH = "Thyroid disease")

## ---- 입력 -----------------------------------------------------------------------------------------------
s2     <- fread(file.path(gwas_dir, "SuppTableS2_FDRsignif.csv"))
fgf5   <- fread(file.path(gwas_dir, "FGF5_rs16998073_PDHT_byGroup.csv"))
counts <- fread(file.path(gwas_dir, "panel_marker_and_test_counts.csv"))
cityr  <- fread(file.path(mqtl_dir, "meQTL_AX83223998_FGF5_CITY_sorted.txt"))
asasr  <- fread(file.path(mqtl_dir, "meQTL_AX83223998_FGF5_ASAS_sorted.txt"))
n_meta <- nrow(fread(file.path(mqtl_dir, "mQTL_meta_CITY_ASAS_all.txt"), select = "gene"))
tab2   <- fread(file.path(mqtl_dir, "Table2_meta_Bonferroni_CpGs.csv"))
setorder(cityr, pvalue); setorder(asasr, pvalue); setorder(tab2, pos)
LEAD    <- tab2$SNP[1]          # 02 에서 mQTL 을 돌린 lead 마커 (AX-83223998)
SNP_POS <- tab2$SNP_POS[1]      # 그 위치 (chr4:81,184,341, hg19)

ax  <- fgf5[MARKER == LEAD]
lead <- ax[COHORT == "CITY_M"]
top_city <- cityr[1]; top_asas <- asasr[gene == top_city$gene]
probe_snp_cpg <- tab2$gene[!is.na(tab2$Probe_rs) & tab2$Probe_rs != ""]   # probe 안 SNP 가 있는 CpG (‡ 표시)

################################################################################################
## Figure 1. Study design
################################################################################################
W <- 7.2; H <- 3.3
for (dev_type in c("pdf", "png")) {
  if (dev_type == "pdf") cairo_pdf(file.path(fig_dir, "Fig1_workflow.pdf"), width = W, height = H, family = FONT, pointsize = 7)
  if (dev_type == "png") png(file.path(fig_dir, "Fig1_workflow.png"), width = W, height = H, units = "in", res = 300,
                             type = "cairo", family = FONT, pointsize = 7)
  par(mai = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")

  text(0.15, 3.15, "Association analysis", adj = c(0, 1), font = 2, cex = 7.8 / 7, col = INK)
  text(3.75, 3.15, "Methylation QTL analysis", adj = c(0, 1), font = 2, cex = 7.8 / 7, col = INK)
  rect(0.10, 0.12, 3.55, 2.92, col = PANEL, border = NA)
  rect(3.70, 0.12, 7.10, 2.92, col = PANEL, border = NA)

  ## 상자 위치 (x, y, 너비, 높이, 강조 여부)
  bx <- data.frame(x = c(0.25, 1.90, 0.25, 3.85, 5.50, 3.85), y = c(1.62, 1.62, 0.27, 1.62, 1.62, 0.27),
                   w = c(1.55, 1.55, 3.20, 1.55, 1.55, 3.20), h = 1.15, hi = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE))
  ## 상자 내용 (숫자는 결과 파일에서 가져옴)
  s_lead <- formatC(lead$P, format = "e", digits = 1)
  s_city <- formatC(top_city$pvalue, format = "e", digits = 1)
  s_asas <- formatC(top_asas$pvalue, format = "e", digits = 1)
  s_meta <- formatC(tab2[gene == top_city$gene]$p_M, format = "e", digits = 1)
  s_bonf <- formatC(0.05 / n_meta, format = "e", digits = 1)
  n_dir  <- sum(ax$BETA > 0)
  box_title <- c("Candidate panel", paste0(length(disease_name), " chronic diseases"),
                 paste0(nrow(s2), " tests at FDR < 0.05 (", uniqueN(s2$DISEASE), " diseases)"),
                 "Discovery: CITY", "Replication: ASAS",
                 paste0("Meta-analysis (", format(n_meta, big.mark = ","), " shared CpGs)"))
  box_body <- list(
    list("36 gene symbols /", "family prefixes (incl. FGF)",
         paste0(format(counts[CHIP == "Kch"]$N_PANEL_MARKERS, big.mark = ","), " K-chip and"),
         paste0(counts[CHIP == "Exo"]$N_PANEL_MARKERS, " Exome markers")),
    list("ASAS and CITY,", "by sex and array;", "logistic regression;", "FDR within strata"),
    list(paste0(sum(s2$BONFERRONI), " pass study-wide Bonferroni."),
         "Lead: rs16998073 (AX-83223998, near FGF5)", "with hypertension, CITY men,",
         bquote("OR" ~ .(sprintf("%.2f", lead$OR)) * "," ~ P == .(sub("e.*", "", s_lead)) %*% 10^.(as.integer(sub(".*e", "", s_lead))) *
                  "; same direction in" ~ .(paste0(n_dir, "/", nrow(ax))) ~ "groups")),
    list(paste0("n = 360; ", format(nrow(cityr), big.mark = ","), " CpGs"),
         paste0(top_city$gene, " (", tab2[gene == top_city$gene]$dist_from_SNP, " bp)"),
         bquote(beta == .(if (top_city$beta < 0) quote(-"") else quote("")) * .(sprintf("%.3f", abs(top_city$beta))) * "," ~
                  P == .(sub("e.*", "", s_city)) %*% 10^.(as.integer(sub(".*e", "", s_city)))),
         if (sum(cityr$FDR < 0.05) == 1) "only CpG at FDR < 0.05" else paste(sum(cityr$FDR < 0.05), "CpGs at FDR < 0.05")),
    list(paste0("n = 280; ", format(nrow(asasr), big.mark = ","), " CpGs"),
         paste0(top_city$gene, if (asasr$gene[1] == top_city$gene) " top CpG" else " (not top CpG)"),
         bquote(beta == .(if (top_asas$beta < 0) quote(-"") else quote("")) * .(sprintf("%.3f", abs(top_asas$beta))) * "," ~
                  P == .(sub("e.*", "", s_asas)) %*% 10^.(as.integer(sub(".*e", "", s_asas))))),
    list(bquote(.(nrow(tab2)) ~ "CpGs pass Bonferroni (" * P < .(sub("e.*", "", s_bonf)) %*% 10^.(as.integer(sub(".*e", "", s_bonf))) * "),"),
         paste0("all within ", max(tab2$pos) - min(tab2$pos), " bp (chr4:", format(min(tab2$pos), big.mark = ","),
                "\u2013", format(max(tab2$pos), big.mark = ","), ");"),
         bquote(.(top_city$gene) ~ "meta" ~ P == .(sub("e.*", "", s_meta)) %*% 10^.(as.integer(sub(".*e", "", s_meta))) *
                  ";" ~ .(if (min(tab2$p_het) > 0.05) "no evidence of heterogeneity" else "heterogeneity present"))))

  for (k in seq_len(nrow(bx))) {
    rect(bx$x[k], bx$y[k], bx$x[k] + bx$w[k], bx$y[k] + bx$h[k],
         col = if (bx$hi[k]) BLUE_T else "white", border = if (bx$hi[k]) BLUE else AXIS, lwd = if (bx$hi[k]) 1.6 else 1.2)
    text(bx$x[k] + 0.1, bx$y[k] + bx$h[k] - 0.1, box_title[k], adj = c(0, 1), font = 2, cex = 7.4 / 7, col = INK)
    for (j in seq_along(box_body[[k]])) {
      text(bx$x[k] + 0.1, bx$y[k] + bx$h[k] - 0.40 - (j - 1) * 0.125, box_body[[k]][[j]], adj = c(0, 0.5), cex = 6.6 / 7, col = INK2)
    }
  }
  arrows(c(1.80, 2.675, 5.40, 5.45, 3.45), c(2.20, 1.62, 2.20, 1.62, 0.85),
         c(1.90, 2.675, 5.50, 5.45, 3.85), c(2.20, 1.42, 2.20, 1.42, 2.05), length = 0.05, col = MUTED, lwd = 1)
  dev.off()
}

################################################################################################
## Figure 2. Panel-restricted association analysis
################################################################################################
g <- copy(s2)
setorder(g, P_FDR)
g[, mlog10FDR := -log10(P_FDR)]
g[, col := ifelse(DISEASE == "PDHT", BLUE, ifelse(DISEASE == "PDDM", ORANGE, MUTED))]
g[, is_fgf5 := MARKER %in% c("AX-83223998", "exm-rs16998073")]
grp <- ax[match(c("ASAS_M", "ASAS_F", "CITY_M", "CITY_F"), COHORT)]

W <- 8.0; H <- 4.3; n <- nrow(g); top <- 3.72; rowh <- 0.262
for (dev_type in c("pdf", "png")) {
  if (dev_type == "pdf") cairo_pdf(file.path(fig_dir, "Fig2_GWAS_panel.pdf"), width = W, height = H, family = FONT, pointsize = 7)
  if (dev_type == "png") png(file.path(fig_dir, "Fig2_GWAS_panel.png"), width = W, height = H, units = "in", res = 300,
                             type = "cairo", family = FONT, pointsize = 7)

  ## (a) lollipop: -log10(FDR)
  l <- 3.0; w <- 2.3; b <- top - n * rowh; h <- n * rowh
  par(mai = c(b, l, H - b - h, W - l - w), xpd = FALSE)
  plot(NA, xlim = c(0, 7.3), ylim = c(n + 0.5, 0.5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  abline(v = c(2, 4, 6), col = GRID, lwd = 0.8)
  abline(v = -log10(0.05), col = MUTED, lty = 2, lwd = 0.9)
  segments(0, seq_len(n), g$mlog10FDR, seq_len(n), col = g$col, lwd = 1.4)
  points(g$mlog10FDR, seq_len(n), pch = 21, bg = g$col, col = "white", cex = 1.5, lwd = 0.8)
  axis(1, at = c(0, 2, 4, 6), col = AXIS, col.axis = INK2, cex.axis = 1, lwd = 0.8, lwd.ticks = 0.8, tcl = -0.3, mgp = c(3, 0.5, 0))
  mtext(expression(-log[10](FDR)), side = 1, line = 2.1, cex = 8 / 7, col = INK)
  par(xpd = NA)
  text(-log10(0.05) + 0.08, 0.5, "FDR = 0.05", adj = c(0, -0.6), cex = 6.5 / 7, col = MUTED)
  s_fdr <- formatC(g$P_FDR[1], format = "e", digits = 1)
  text(g$mlog10FDR[1] - 0.1, 1.5, bquote(FDR == .(sub("e.*", "", s_fdr)) %*% 10^.(as.integer(sub(".*e", "", s_fdr)))),
       adj = c(1, 0.5), cex = 6.8 / 7, col = INK2)

  ## (b) 4개 그룹 -log10(P)
  l2 <- 6.3; w2 <- 1.1; h2 <- 1.15; b2 <- top - h2
  par(new = TRUE, mai = c(b2, l2, H - b2 - h2, W - l2 - w2), xpd = FALSE)
  plot(NA, xlim = c(0, 16), ylim = c(4.5, 0.5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  abline(v = c(5, 10), col = GRID, lwd = 0.8)
  rect(0, 1:4 - 0.23, -log10(grp$P), 1:4 + 0.23, col = ifelse(grp$P < 0.05, BLUE, BLUE_T), border = BLUE, lwd = 0.9)
  abline(v = -log10(0.05), col = MUTED, lty = 2, lwd = 0.9)
  axis(1, at = c(0, 5, 10), col = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.3, mgp = c(3, 0.5, 0))
  mtext(expression(-log[10](P)), side = 1, line = 2.1, cex = 8 / 7, col = INK)
  par(xpd = NA)
  for (i in 1:4) {
    s_p <- formatC(grp$P[i], format = "e", digits = 1)
    lab <- if (grp$P[i] < 0.01) bquote(P == .(sub("e.*", "", s_p)) %*% 10^.(as.integer(sub(".*e", "", s_p)))) else
      paste0("P = ", sprintf("%.3f", grp$P[i]), if (grp$P[i] >= 0.05) " (n.s.)" else "")
    text(-log10(grp$P[i]) + 0.4, i, lab, adj = c(0, 0.5), cex = 6.6 / 7, col = INK2)
    text(-0.7, i, paste(sub("_.*", "", grp$COHORT[i]), "\u00b7", sub(".*_", "", grp$COHORT[i])), adj = c(1, 0.5), col = INK2)
  }

  ## 바깥 글자 (inch 좌표)
  par(new = TRUE, mai = c(0, 0, 0, 0)); plot.new(); plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")
  text(0.03, H - 0.15, "a", adj = c(0, 1), font = 2, cex = 11 / 7, col = INK)
  text(0.38, H - 0.19, paste0("Panel-restricted association analysis: ", n, " marker\u2013disease\u2013group tests at FDR < 0.05"),
       adj = c(0, 1), cex = 8 / 7, col = INK2)
  text(c(0.1, 1.3, 2.35), top + 0.06, c("Marker", "Disease", "Group"), adj = c(0, 0), cex = 6.8 / 7, col = MUTED)
  yy <- top - (seq_len(n) - 0.5) * rowh
  text(0.1, yy, g$MARKER, adj = c(0, 0.5), font = ifelse(g$is_fgf5, 2, 1), col = INK)
  text(1.3, yy, disease_name[g$DISEASE], adj = c(0, 0.5), col = INK2)
  text(2.35, yy, paste(g$COHORT_BASE, "\u00b7", g$SEX), adj = c(0, 0.5), col = INK2)
  points(rep(5.55, 3), c(0.86, 0.64, 0.43), pch = 21, bg = c(BLUE, ORANGE, MUTED), col = "white", cex = 1.4)
  text(5.68, c(0.86, 0.64, 0.43), c("Hypertension", "Diabetes mellitus", "Other diseases"), adj = c(0, 0.5), cex = 6.8 / 7, col = INK2)
  text(0.1, 0.11, "Bold: rs16998073 (FGF5) \u2014 AX-83223998 on K-chip, exm-rs16998073 on Exome chip", adj = c(0, 0), cex = 6.5 / 7, col = MUTED)
  text(5.55, H - 0.15, "b", adj = c(0, 1), font = 2, cex = 11 / 7, col = INK)
  text(5.75, H - 0.19, "AX-83223998 (rs16998073)\nby cohort and sex, hypertension", adj = c(0, 1), cex = 8 / 7, col = INK2)
  text(5.55, b2 - 0.55, paste0("Dashed line: P = 0.05.\nEffect direction concordant\n(\u03b2 > 0) in ",
                               if (all(grp$BETA > 0)) "all four groups;" else paste0(sum(grp$BETA > 0), " of 4 groups;"),
                               "\nP < 0.05 in ", sum(grp$P < 0.05), " of 4."), adj = c(0, 1), cex = 7 / 7, col = INK2)
  dev.off()
}

################################################################################################
## Figure 3. Discovery, replication and location of the FGF5 mQTL
################################################################################################
mq <- cityr[1:20]
mq[, rank := 1:20]
mq[, mlog10P := -log10(pvalue)]
mq[, hi := gene %in% tab2$gene]
fdr_line <- -log10(0.05 / nrow(cityr))       # FDR = 0.05 에 해당하는 1위 P (BH)

W <- 7.2; H <- 7.8
for (dev_type in c("pdf", "png")) {
  if (dev_type == "pdf") cairo_pdf(file.path(fig_dir, "Fig3_mQTL_FGF5_locus.pdf"), width = W, height = H, family = FONT, pointsize = 7)
  if (dev_type == "png") png(file.path(fig_dir, "Fig3_mQTL_FGF5_locus.png"), width = W, height = H, units = "in", res = 300,
                             type = "cairo", family = FONT, pointsize = 7)

  ## (a) CITY 상위 20 CpG
  l <- 0.7; b <- 5.5; w <- 6.3; h <- 1.75
  par(mai = c(b, l, H - b - h, W - l - w), xpd = FALSE)
  plot(NA, xlim = c(0.3, 20.7), ylim = c(0, 10), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  abline(h = c(2, 4, 6, 8), col = GRID, lwd = 0.8)
  abline(h = fdr_line, col = MUTED, lty = 2, lwd = 0.9)
  segments(mq$rank, 0, mq$rank, mq$mlog10P, col = ifelse(mq$hi, BLUE, MUTED), lwd = 1.4)
  points(mq$rank, mq$mlog10P, pch = ifelse(mq$beta < 0, 25, 24), bg = ifelse(mq$hi, BLUE, MUTED), col = "white",
         cex = ifelse(mq$hi, 1.5, 1.2), lwd = 0.6)
  axis(1, at = 1:20, col = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.3, mgp = c(3, 0.5, 0))
  axis(2, at = c(0, 2, 4, 6, 8), col = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.3, mgp = c(3, 0.6, 0), las = 1)
  mtext("CpG rank by P-value in CITY (top 20 shown)", side = 1, line = 2.1, cex = 8 / 7, col = INK)
  mtext(expression(-log[10](P)), side = 2, line = 2.1, cex = 8 / 7, col = INK)
  par(xpd = NA)
  text(20.6, fdr_line - 0.15, paste0("FDR = 0.05 (", format(nrow(cityr), big.mark = ","), " CpGs)"), adj = c(1, 1), cex = 6.5 / 7, col = MUTED)
  s_p <- formatC(mq$pvalue[1], format = "e", digits = 1)
  text(1.7, 8.9, mq$gene[1], adj = c(0, 0.5), font = 2, cex = 8 / 7, col = BLUE)
  text(1.7, 8.1, bquote(P == .(sub("e.*", "", s_p)) %*% 10^.(as.integer(sub(".*e", "", s_p))) * "," ~~ FDR == .(sprintf("%.4f", mq$FDR[1]))),
       adj = c(0, 0.5), cex = 7 / 7, col = INK2)
  for (i in 2:20) {
    if (mq$hi[i]) text(i + 0.3, mq$mlog10P[i] + 0.15, paste0(mq$gene[i], if (mq$gene[i] %in% probe_snp_cpg) " \u2021" else ""),
                       adj = c(0, 0.5), cex = 6.8 / 7, col = BLUE)
  }
  points(c(14.6, 16.4, 18.2), rep(9.4, 3), pch = c(25, 24, 25), bg = c(MUTED, MUTED, BLUE), col = "white", cex = 1.2)
  text(c(14.6, 16.4, 18.2) + 0.35, 9.4, c("\u03b2 < 0", "\u03b2 > 0", "meta-analysis sig."), adj = c(0, 0.5), cex = 6.8 / 7, col = INK2)

  ## (b) forest plot: CITY / ASAS / Meta
  l2 <- 1.75; b2 <- 2.75; w2 <- 3.1; h2 <- 1.75
  par(new = TRUE, mai = c(b2, l2, H - b2 - h2, W - l2 - w2), xpd = FALSE)
  plot(NA, xlim = c(-0.26, 0), ylim = c(10.5, -0.5), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  abline(v = c(-0.25, -0.20, -0.15, -0.10, -0.05), col = GRID, lwd = 0.8)
  abline(v = 0, col = AXIS, lwd = 0.8)
  fr <- data.frame()
  for (k in seq_len(nrow(tab2))) {
    fr <- rbind(fr, data.frame(y = (k - 1) * 4 + 0:2, lab = c("CITY", "ASAS", "Meta"),
                               beta = c(tab2$beta_C[k], tab2$beta_A[k], tab2$beta_M[k]),
                               se = c(tab2$se_C[k], tab2$se_A[k], tab2$se_M[k]),
                               p = c(tab2$p_C[k], tab2$p_A[k], tab2$p_M[k])))
  }
  fr$meta <- fr$lab == "Meta"
  segments(fr$beta - 1.96 * fr$se, fr$y, fr$beta + 1.96 * fr$se, fr$y, col = ifelse(fr$meta, BLUE, INK2), lwd = ifelse(fr$meta, 1.8, 1.2))
  points(fr$beta, fr$y, pch = ifelse(fr$meta, 23, ifelse(fr$lab == "CITY", 21, 22)), bg = ifelse(fr$meta, BLUE, INK2),
         col = "white", cex = ifelse(fr$meta, 1.8, 1.3), lwd = 0.6)
  axis(1, at = c(-0.25, -0.20, -0.15, -0.10, -0.05, 0), labels = c(paste0(MINUS, c("0.25", "0.20", "0.15", "0.10", "0.05")), "0"),
       col = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.3, mgp = c(3, 0.5, 0))
  mtext("\u03b2 per allele (methylation M-value), 95% CI", side = 1, line = 2.1, cex = 8 / 7, col = INK)
  par(xpd = NA)
  for (r in seq_len(nrow(fr))) {
    s_p <- formatC(fr$p[r], format = "e", digits = 1)
    text(0.035, fr$y[r], bquote(P == .(sub("e.*", "", s_p)) %*% 10^.(as.integer(sub(".*e", "", s_p)))),
         adj = c(0, 0.5), cex = 6.5 / 7, col = if (fr$meta[r]) INK else INK2, font = if (fr$meta[r]) 2 else 1)
    text(-0.272, fr$y[r], fr$lab[r], adj = c(1, 0.5), cex = 6.8 / 7, col = INK2)
  }
  for (k in seq_len(nrow(tab2))) {
    text(-0.38, (k - 1) * 4 + 1, paste0(tab2$gene[k], if (tab2$gene[k] %in% probe_snp_cpg) " \u2021" else ""),
         adj = c(0, 0.5), font = 2, col = INK)
  }

  ## (c) 유전자좌 지도
  l3 <- 0.7; b3 <- 0.55; w3 <- 6.3; h3 <- 1.5
  x0 <- 81183900; x1 <- 81189200
  par(new = TRUE, mai = c(b3, l3, H - b3 - h3, W - l3 - w3), xpd = FALSE)
  plot(NA, xlim = c(x0, x1), ylim = c(-1.7, 1.9), xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
  segments(x0, 0, FGF5_TSS, 0, col = AXIS, lwd = 1.2)
  polygon(c(FGF5_TSS, x1 - 150, x1, x1 - 150, FGF5_TSS), c(-0.42, -0.42, 0, 0.42, 0.42), col = "#dedcd4", border = MUTED, lwd = 0.8)
  text((FGF5_TSS + x1) / 2, 0, "FGF5", font = 4, cex = 8 / 7, col = INK2)
  rect(min(tab2$pos) - 30, -0.22, max(tab2$pos) + 30, 0.22, col = BLUE_T, border = NA)
  points(SNP_POS, 0, pch = 21, bg = BLUE, col = "white", cex = 1.8, lwd = 0.8)
  points(tab2$pos, rep(0, nrow(tab2)), pch = 22, bg = "white", col = BLUE, cex = 1.6, lwd = 1.4)
  axis(1, at = seq(81184000, 81189000, by = 1000), labels = format(seq(81184000, 81189000, by = 1000), big.mark = ","),
       col = AXIS, col.axis = INK2, lwd = 0.8, tcl = -0.3, mgp = c(3, 0.5, 0))
  mtext("chr4 (hg19), bp", side = 1, line = 2.0, cex = 8 / 7, col = INK)
  par(xpd = NA)
  segments(SNP_POS, 0.12, SNP_POS, 1.15, col = AXIS, lwd = 0.6)
  text(SNP_POS, 1.2, paste0(LEAD, " (rs16998073)\nchr4:", format(SNP_POS, big.mark = ",")), adj = c(0.5, 0), font = 2, cex = 7 / 7, col = INK)
  segments(min(tab2$pos), 0.55, max(tab2$pos) + 70, 0.85, col = AXIS, lwd = 0.6)
  text(max(tab2$pos) + 80, 0.95, paste0(nrow(tab2), " meta-analysis CpGs\nchr4:", format(min(tab2$pos), big.mark = ","), "\u2013",
                                        format(max(tab2$pos), big.mark = ","), "\n(", max(tab2$pos) - min(tab2$pos), " bp window)"),
       adj = c(0, 0.5), font = 2, cex = 7 / 7, col = INK)
  arrows(SNP_POS, -0.75, max(tab2$pos), -0.75, code = 3, length = 0.04, col = INK2, lwd = 0.8)
  text((SNP_POS + max(tab2$pos)) / 2, -0.9, paste0(min(tab2$dist_from_SNP), "\u2013", max(tab2$dist_from_SNP), " bp from SNP"),
       adj = c(0.5, 1), cex = 7 / 7, col = INK2)
  arrows(max(tab2$pos), -0.75, FGF5_TSS, -0.75, code = 3, length = 0.04, col = MUTED, lwd = 0.7)
  text((max(tab2$pos) + FGF5_TSS) / 2, -0.9, paste0("\u2248 ", sprintf("%.1f", (FGF5_TSS - max(tab2$pos)) / 1000), " kb to FGF5 start"),
       adj = c(0.5, 1), cex = 7 / 7, col = MUTED)

  ## 바깥 글자 (inch 좌표)
  par(new = TRUE, mai = c(0, 0, 0, 0)); plot.new(); plot.window(xlim = c(0, W), ylim = c(0, H), xaxs = "i", yaxs = "i")
  text(0.03, H - 0.10, "a", adj = c(0, 1), font = 2, cex = 11 / 7, col = INK)
  text(0.30, H - 0.14, "Discovery: mQTL of AX-83223998 (rs16998073) against genome-wide CpGs in CITY", adj = c(0, 1), cex = 8 / 7, col = INK2)
  text(0.03, 4.98, "b", adj = c(0, 1), font = 2, cex = 11 / 7, col = INK)
  text(0.30, 4.94, bquote("Replication in ASAS (" * .(format(nrow(asasr), big.mark = ",")) ~ "CpGs) and fixed-effect meta-analysis (" *
                            .(format(n_meta, big.mark = ",")) ~ "shared CpGs;" ~ P[het] >= .(sprintf("%.2f", min(tab2$p_het))) * ")"),
       adj = c(0, 1), cex = 8 / 7, col = INK2)
  text(0.03, 2.38, "c", adj = c(0, 1), font = 2, cex = 11 / 7, col = INK)
  text(0.30, 2.34, paste0("Genomic context (hg19; FGF5 on + strand, start ", format(FGF5_TSS, big.mark = ","), "; distances to scale)"),
       adj = c(0, 1), cex = 8 / 7, col = INK2)
  text(0.30, 0.04, "\u2021 Probe contains a known SNP (rs59133384); interpret with caution.", adj = c(0, 0), cex = 6.5 / 7, col = MUTED)
  dev.off()
}

cat(">>> 그림 저장 위치:", fig_dir, "\n")
cat("\n>>> 전체 소요 시간:", format(round(difftime(Sys.time(), t_start, units = "mins"), 1)), "\n")
sessionInfo()
