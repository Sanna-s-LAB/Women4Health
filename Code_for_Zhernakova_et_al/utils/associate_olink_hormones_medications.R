library(future) 
library(future.apply)

set.seed(123)

####################
### Remove 6 women taking thyroid hormones and re-run the protein - hormone associations
####################
ids_rm <- c("S002", "S023", "X008", "X067", "X111", "X082")
d_wide_2 <- d_wide[! d_wide$ID %in% ids_rm,]

all_hormones <- c("PROG", "FSH", "17BES", "LH", "PRL")
all_prots <- colnames(d_wide_2)[! colnames(d_wide_2) %in% c("SampleID", "ID", "phase", "TP")]

# Create a grid of all combinations to process
tasks <- expand.grid(prot = all_prots, ph = all_hormones, stringsAsFactors = FALSE)

# Setup parallel backend (use most available cores)
plan(multisession, workers = parallel::detectCores() - 1)

# run calculations in parallel
results_list <- future_lapply(1:nrow(tasks), function(idx) {
  prot <- tasks$prot[idx]
  ph <- tasks$ph[idx]
  
  res_gam <- gam_prot_pheno_adj_covar(d_wide_2, pheno, prot, ph, covariates, 
                                      scale = T, adjust_timepoint = 'spline', 
                                      anova_pval = F, longitudinal = T)

  return(c(prot = prot, pheno = ph, unlist(res_gam))) 
}, future.seed = TRUE) # Important for reproducibility/random processes in GAMs

# Combine results
gam_res_no_thyr <- as.data.frame(do.call(rbind, results_list))
colnames(gam_res_no_thyr) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples")

# Clean up
gam_res_no_thyr <- gam_res_no_thyr %>%
  mutate(across(c(pval, estimate, SE, n, n_samples), as.numeric)) %>% 
  na.omit() %>%
  arrange(pval) %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH'))


write.table(gam_res_no_thyr, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.no_thyroid_med.txt"), 
            quote = F, sep = "\t", row.names = FALSE)


m <- full_join(gam_res_hormones, gam_res_no_thyr, by = c("pheno", "prot"), suffix = c(".base", ".no_thyroid_med"))
m$sign <- "not significant any"
m[m$BH_pval.base < 0.05 & m$BH_pval.no_thyroid_med < 0.05,"sign"] <- "significant both"
m[m$BH_pval.base < 0.05 & m$BH_pval.no_thyroid_med > 0.05,"sign"] <- "significant base"
m[m$BH_pval.base > 0.05 & m$BH_pval.no_thyroid_med < 0.05,"sign"] <- "significant no thyroid med"
m$sign <- as.factor(m$sign)

write.table(m, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.CMP.thyroid_med.txt"), 
            quote = F, sep = "\t", row.names = FALSE)


cor <- cor(m$estimate.base, m$estimate.no_thyroid_med)
# [1] 0.9879248

p1 <- ggplot(m, aes(x = estimate.base, y = estimate.no_thyroid_med)) + 
  geom_abline(slope = 1, color = 'grey', lty = 2) +
  geom_point(alpha = 0.5) + 
  theme_minimal() + ggtitle(paste("Pearson r = ", formatC(cor, digits = 3))) +
  ylab("estimate when removing individuals taking thyroid hormone medication") + xlab("estimate base model") +
  ggtitle("Hormone - protein associations\ncomparison of effect estimates")

p2 <- ggplot(m, aes(x = -log10(pval.base), y = -log10(pval.no_thyroid_med))) + 
  #geom_rect(aes(xmin = 0, xmax = 20, ymin = 0, ymax = 20), fill = NA, color = "grey", linetype = "dashed")  +
  geom_abline(slope = 1, color = 'grey', lty = 2) +
  geom_point(alpha = 0.5) + 
  theme_minimal() +
  ylab("-log10(P) when removing individuals taking thyroid hormone medication") + xlab("-log10(P) base model") +
  ggtitle("Hormone - protein associations\ncomparison of -log10(P)")

p3 <- ggplot(m, aes(x = -log10(pval.base), y = -log10(pval.no_thyroid_med), color = sign)) + 
  geom_abline(slope = 1, color = 'grey', lty = 2) +
  geom_point(alpha = 0.5) + 
  theme_minimal() + ylim(c(0,20)) + xlim(c(0,20)) +
  scale_colour_brewer(palette = "Dark2") +
  guides(color=guide_legend(title="BH significance")) +
  ylab("-log10(P) when removing individuals taking thyroid hormone medication") + xlab("-log10(P) base model") +
  ggtitle("Zoom into the comparison of -log10(P) ")

pdf(paste0(out_basedir, "plots/prot_vs_hormones.cmp_thyroid_med.pdf"), height = 12, width = 12)
(p2 + p1) / (p3 + plot_spacer())
dev.off()


######
# Adjusting for medication use in the association analysis
#####

med <- read.delim("../../phenotypes/batch12/medications_per_visit.txt", as.is = T, check.names = F, sep = "\t")
med$Supp_check = med$Med_check = med$ID = med$Med_or_supp_check = med$integratori_names = NULL
med$farmaci_names = med$pathology= NULL
med$paracetamol <- grepl("Paracetamol", med$drug_category, fixed = T)
med$NSAID <- grepl("NSAID", med$drug_category, fixed = T)
med$Antihistamine <- grepl("Antihistamine", med$drug_category, fixed = T)
med$drug_category = NULL

covar_tmp <- covariates
covar_tmp$SampleID <- paste0(covar_tmp$ID, "_", as.numeric(covar_tmp$phase))
covar_tmp$TP <- as.numeric(covar_tmp$phase)
covar_tmp$phase <- NULL
covar_med <- inner_join(covar_tmp, med, by = c("SampleID" = "Code"))

tasks <- expand.grid(prot = all_prots, ph = all_hormones, stringsAsFactors = FALSE)

#  Setup parallel backend (use most available cores)
plan(multisession, workers = parallel::detectCores() - 1)

# run calculations in parallel
results_list <- future_lapply(1:nrow(tasks), function(idx) {
  prot <- tasks$prot[idx]
  ph <- tasks$ph[idx]

  res_gam <- gam_prot_pheno_adj_covar(d_wide, pheno, prot, ph, covar_med, 
                                      scale = T, adjust_timepoint = 'spline', 
                                      anova_pval = F, longitudinal = T)
  
  return(c(prot = prot, pheno = ph, unlist(res_gam))) 
}, future.seed = TRUE) # Important for reproducibility/random processes in GAMs

# Combine results
gam_res_med <- as.data.frame(do.call(rbind, results_list))
colnames(gam_res_med) <- c("prot", "pheno", "pval", "estimate", "SE", "n", "n_samples")

# Clean up
gam_res_med <- gam_res_med %>%
  mutate(across(c(pval, estimate, SE, n, n_samples), as.numeric)) %>% 
  na.omit() %>%
  arrange(pval) %>%
  mutate(BH_pval = p.adjust(pval, method = 'BH'))

write.table(gam_res_med, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.MED.txt"), 
            quote = F, sep = "\t", row.names = FALSE)

m <- full_join(gam_res_hormones, gam_res_med, by = c("pheno", "prot"), suffix = c(".base", ".adj_med"))
m$sign <- "not significant any"
m[m$BH_pval.base < 0.05 & m$BH_pval.adj_med < 0.05,"sign"] <- "significant both"
m[m$BH_pval.base < 0.05 & m$BH_pval.adj_med > 0.05,"sign"] <- "significant base"
m[m$BH_pval.base > 0.05 & m$BH_pval.adj_med < 0.05,"sign"] <- "significant no thyroid med"
m$sign <- as.factor(m$sign)

write.table(m, file = paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.CMP.MED.txt"), 
            quote = F, sep = "\t", row.names = FALSE)


cor <- cor(m$estimate.base, m$estimate.adj_med)
# [1] 0.9662154
p1 <- ggplot(m, aes(x = estimate.base, y = estimate.adj_med)) + 
  geom_abline(slope = 1, color = 'grey', lty = 2) +
  geom_point(alpha = 0.5) + 
  theme_minimal() + ggtitle(paste("Pearson r = ", formatC(cor, digits = 3))) +
  ylab("estimate when adjusting for medications") + xlab("estimate base model") +
  ggtitle("Hormone - protein associations\ncomparison of effect estimates")

p2 <- ggplot(m, aes(x = -log10(pval.base), y = -log10(pval.adj_med))) + 
  #geom_rect(aes(xmin = 0, xmax = 20, ymin = 0, ymax = 20), fill = NA, color = "grey", linetype = "dashed")  +
  geom_abline(slope = 1, color = 'grey', lty = 2) +
  geom_point(alpha = 0.5) + 
  theme_minimal() +
  ylab("-log10(P) when when adjusting for medications") + xlab("-log10(P) base model") +
  ggtitle("Hormone - protein associations\ncomparison of -log10(P)")


pdf(paste0(out_basedir, "plots/prot_vs_hormones.cmp_adj_med.pdf"), height = 7, width = 12)
(p2 + p1) 
dev.off()


