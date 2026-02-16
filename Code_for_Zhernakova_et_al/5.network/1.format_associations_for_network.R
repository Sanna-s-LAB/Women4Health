library(dplyr)
set.seed(123)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
out_basedir <- "results12/intensity_shared_prots_261125/"

gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.shared_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_phenos<- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.shared_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- read.delim(paste0(out_basedir, "all_pheno_vs_all_pheno_spline_gam.shared_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- gam_res_hormones_pheno[gam_res_hormones_pheno$prot != gam_res_hormones_pheno$pheno,]

causal_links <- read.delim(paste0(out_basedir, "causality/causal_links.v2.txt"), as.is = T, sep = "\t", check.names  = FALSE)


network_data <- rbind(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_phenos[gam_res_phenos$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_hormones_pheno[gam_res_hormones_pheno$BH_pval < 0.05, c("prot", "pheno", "estimate")])

causal_forward <- paste(causal_links$cause, causal_links$consequence)

# Fix the order of the nodes to enable arrows later
network_data_swapped <- network_data %>%
  mutate(
    # check if there is direcion 
    is_forward = paste(prot, pheno) %in% causal_forward,
    is_reverse = paste(pheno, prot) %in% causal_forward,
    
    # check if association has defined direction
    has_direction = is_forward | is_reverse,
    
    # if bi-directional we need 2 rows, if not - 1 row
    n_rows = ifelse(is_forward & is_reverse, 2, 1)
  ) %>%
  # add a duplicate row for bi-directional
  uncount(n_rows, .id = "row_id") %>%
  mutate(
    # Swap the nodes if:
    # reverse relationship (Reverse == T, Forward == F)
    # or it's the second copy of a bi-directional relationship (row_id == 2)
    should_swap = (is_reverse & !is_forward) | (row_id == 2),
    
    prot_new = ifelse(should_swap, pheno, prot),
    pheno_new = ifelse(should_swap, prot, pheno)
  ) %>%
  select(prot = prot_new, pheno = pheno_new, estimate, has_direction) %>%
  distinct()


nodes_data <- unique(rbind(data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"pheno"], type = "hormone"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"pheno"], type = "phenotype")))


write.table(network_data_swapped, file = paste0(out_basedir, "network_data/network.edges.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(nodes_data, file = paste0(out_basedir, "network_data/network.nodes.txt"), quote = F, sep = "\t", row.names = FALSE)
