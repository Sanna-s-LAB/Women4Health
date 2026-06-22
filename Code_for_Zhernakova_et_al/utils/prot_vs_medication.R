#get data
med <- read.csv(
  "/mnt/sannaLAB-Temp/dasha/olink/batch2/data/medications_per_visit.txt",
  sep = "\t"
)

prot <- read.csv(
  "/mnt/sannaLAB-Temp/dasha/olink/batch2/data/olink_batch12.intensity.bridged_all_proteins_lod150_wide_rm_outliers_4sd.adj_batch_storage.txt",
  sep = "\t"
)

rownames(prot) <- prot$SampleID
prot <- prot[med$Code, ]


# add missing ids
prob_id <- c()
for(id in c("S002", "S023", "X008", "X067", "X111", "X082")){
  if(any(!grepl('Thyroid',med[med$ID == id,'drug_category']))){
    prob_id <- c(prob_id,id)
  }
}
problem <- med[med$ID %in% prob_id,'drug_category']


med[med$ID %in% prob_id,'drug_category'] <- ifelse(nchar(problem) == 0,paste0('Thyroid hormone ',problem),paste0('Thyroid hormone, ',problem))
med[med$ID %in% prob_id,'Med_check'] <- 'Yes'


library(vegan)

med[is.na(med$Med_check),'Med_check'] <- 'No'

dist_prote <- dist(prot[, -c(1:3)])

# ==========================================
# TEST 1: medications Yes/No
# ==========================================

test_medcheck <- adonis2(
  dist_prote ~ Med_check,
  data = med,
  strata = med$ID,
  permutations = 999
)

print(test_medcheck)

# ==========================================
# TEST 2: drug category
# ==========================================

idx_med <- med$Med_check == 'Yes'

med_drug  <- med[idx_med, ]
prot_drug <- prot[idx_med, ]

dist_prote_drug <- dist(prot_drug[, -c(1:3)])

split_classes <- strsplit(
  as.character(med_drug$drug_category),
  ",\\s*"
)

all_classes <- trimws(unlist(split_classes))
all_classes <- all_classes[all_classes != ""]

freq <- sort(table(all_classes), decreasing = TRUE)

print(freq)

keep_drugs <- names(freq[freq >= 14])

print(keep_drugs)

drug_mat <- sapply(
  keep_drugs,
  function(cl)
    sapply(split_classes, function(x) cl %in% trimws(x))
)

drug_mat <- as.data.frame(drug_mat)
colnames(drug_mat) <- make.names(colnames(drug_mat))

med_drug <- cbind(med_drug, drug_mat)

form <- as.formula(
  paste(
    "dist_prote_drug ~",
    paste(colnames(drug_mat), collapse = " + ")
  )
)

test_drugclass <- adonis2(
  form,
  data = med_drug,
  strata = med_drug$ID,
  permutations = 999,
  by = "margin"
)

#drug category statistics per woman
n_per_farm <- rep(0, length(colnames(drug_mat)))
names(n_per_farm) <- colnames(drug_mat)

for(farm in unique(colnames(drug_mat))) {
  
  ids_farm <- c()
  
  for(id in unique(med_drug$ID)) {
    med_temp <- med_drug[med_drug$ID == id, ]
    
    if(any(med_temp[[farm]], na.rm = TRUE)) {
      ids_farm <- c(ids_farm, id)
      n_per_farm[[farm]] <- n_per_farm[[farm]] + 1
    }
  }
  
  if(length(ids_farm) > 0) {
    cat(farm, ":", paste(ids_farm, collapse = ", "), "\n")
  }
}

