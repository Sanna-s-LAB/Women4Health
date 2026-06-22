setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")

# Read W4H protein vs phase results
d <- read.delim(paste0(out_basedir, "/prot_vs_tp_gam.txt"), as.is = T, check.names = F, sep = "\t")
d$signif <- d$gam_BH_pval < 0.05

# read UKB results
ukb <- read.delim("../UKB_riishede_results.txt", as.is = T, check.names = F, sep = "\t")

overlap_prots <- intersect(ukb$Protein, d$prot)

ukb <- ukb[ukb$Protein %in% overlap_prots,]
d <- d[d$prot %in% overlap_prots,]

merged <- full_join(d, ukb, by = c("prot" = "Protein"), suffix = c("_w4h", "_ukb"))

nrow(d[d$signif == T,])
nrow(merged[merged$signif ==T & merged$P < 0.05,])

nrow(ukb[ukb$FDR < 0.05,])
nrow(merged[merged$FDR < 0.05 & merged$gam_pval < 0.05,])


ggplot(merged, aes(x = -log10(gam_pval), y = -log10(P))) + geom_point() + ylim(0,10) + xlim(0,10)
