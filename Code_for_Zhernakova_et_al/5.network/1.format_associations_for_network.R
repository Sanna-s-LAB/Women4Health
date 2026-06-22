library(dplyr)
set.seed(123)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12")
out_basedir <- "results12/intensity_all_prots_220526/"

gam_res_hormones <- read.delim(paste0(out_basedir, "prot_vs_hormones.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_phenos<- read.delim(paste0(out_basedir, "prot_vs_phenotypes.gam.spline.all_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- read.delim(paste0(out_basedir, "all_pheno_vs_all_pheno_spline_gam.shared_prots.txt"), as.is = T, sep = "\t", check.names  = FALSE)
gam_res_hormones_pheno <- gam_res_hormones_pheno[gam_res_hormones_pheno$prot != gam_res_hormones_pheno$pheno,]
#colnames(gam_res_hormones_pheno)[c(1,2)] <-c("prot", "pheno")

gam_res_hormones$prot <- gsub("PRL", "PRL_olink", gam_res_hormones$prot)
gam_res_phenos$prot <- gsub("PRL", "PRL_olink", gam_res_phenos$prot)  

causal_1 <- read.delim(paste0(out_basedir, "causality/hormone_vs_phenotypes.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)
causal_2 <- read.delim(paste0(out_basedir, "causality/prot_vs_hormones.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)
causal_3 <- read.delim(paste0(out_basedir, "causality/prot_vs_phenotypes.causality_wide.linear.txt"), as.is = T, sep = "\t", check.names  = FALSE)

# rm hormone - pheno pairs whose association is not nominally significant:
sign_pairs <- c(
  paste(gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$prot, gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$pheno, sep = "_"),
  paste(gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$pheno, gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05, ]$prot, sep = "_")
)
causal_1 <- causal_1[paste(causal_1$phenotype, causal_1$protein, sep = "_") %in% sign_pairs,]

causal <- rbind(causal_1, causal_2, causal_3)

causal_links <- causal[,1:3] %>%
  separate_rows(summary, sep = ";") %>%               # one row per timepoint
  mutate(summary = str_remove(summary, "^\\d+:")) %>% # remove "1:", "2:", etc.
  filter(summary != "nothing significant") %>%        # drop empty timepoints
  rowwise() %>%
  mutate(
    pairs = if (summary == "both directions significant") {
      list(c(paste(phenotype, protein, sep = " → "),
             paste(protein, phenotype, sep = " → ")))
    } else {
      list(summary)
    }
  ) %>%
  unnest(pairs) %>%
  separate(pairs, into = c("cause", "consequence"), sep = " → ") %>%
  select(cause, consequence) %>%
  unique()

write.table(causal_links, file = paste0(out_basedir, "causality/causal_links.v2.txt"), quote = F, sep = "\t", row.names = F)



network_data <- rbind(gam_res_hormones[gam_res_hormones$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_phenos[gam_res_phenos$BH_pval < 0.05, c("prot", "pheno", "estimate")],
                      gam_res_hormones_pheno[gam_res_hormones_pheno$pval < 0.05 , c("prot", "pheno", "estimate")])

# add hormone-pheno links that do not have a significant association, but do have a significant causal link
additional_hormones_pheno <- rbind(
  inner_join(gam_res_hormones_pheno[c("prot", "pheno", "estimate")], causal_links, by = c("prot" = "cause", "pheno" = "consequence")),
  inner_join(gam_res_hormones_pheno[c("prot", "pheno", "estimate")], causal_links, by = c("prot" = "consequence", "pheno" = "cause")))

network_data <- rbind(network_data, additional_hormones_pheno)


causal_forward <- paste(causal_links$cause, causal_links$consequence)

network_data_swapped <- network_data %>%
  mutate(
    # Check existence of both directions
    is_forward = paste(prot, pheno) %in% causal_forward,
    is_reverse = paste(pheno, prot) %in% causal_forward,
    
    # Determine if this edge has a defined direction
    has_direction = is_forward | is_reverse,
    
    # Determine how many rows we need.
    # If it is bi-directional (both forward and reverse exist), we need 2 rows.
    # Otherwise, we keep 1 row.
    n_rows = ifelse(is_forward & is_reverse, 2, 1)
  ) %>%
  # Expand the dataframe: this duplicates the rows where n_rows is 2.
  # .id = "row_id" creates a column (1 or 2) to distinguish the copies.
  uncount(n_rows, .id = "row_id") %>%
  mutate(
    # Logic to determine if we swap the nodes:
    # Swap if:
    # 1. It is strictly a reverse relationship (Reverse is T, Forward is F)
    # 2. OR It is the second copy of a bi-directional relationship (row_id == 2)
    should_swap = (is_reverse & !is_forward) | (row_id == 2),
    
    prot_new = ifelse(should_swap, pheno, prot),
    pheno_new = ifelse(should_swap, prot, pheno)
  ) %>%
  # Select and rename final columns
  dplyr::select(prot = prot_new, pheno = pheno_new, estimate, has_direction) %>%
  # Ensure no exact duplicates 
  distinct()


nodes_data <- unique(rbind(data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_hormones[gam_res_hormones$BH_pval < 0.05,"pheno"], type = "hormone"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"prot"], type = "protein"),
                           data.frame(feature = gam_res_phenos[gam_res_phenos$BH_pval < 0.05,"pheno"], type = "phenotype")))

nodes_data$feature <- rename_p4_e2(nodes_data$feature)
network_data_swapped$prot <- rename_p4_e2(network_data_swapped$prot)
network_data_swapped$pheno <- rename_p4_e2(network_data_swapped$pheno)


write.table(network_data_swapped, file = paste0(out_basedir, "network/network.spline.edges.causality2.with_pheno-pheno.txt"), quote = F, sep = "\t", row.names = FALSE)
write.table(nodes_data, file = paste0(out_basedir, "network/network.spline.nodes.txt"), quote = F, sep = "\t", row.names = FALSE)

