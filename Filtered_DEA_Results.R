# Filtered DEA Results for IMP Exploration
# Author: Daniel Guevara

# 1. Call Libraries
library(openxlsx)
library(dplyr)

# 2. Read data
sampledata <- read.xlsx("~/Documents/SimonLab/Models/RU_Paired_Samples.xlsx", colNames = T)
normCounts <- readRDS("~/Documents/SimonLab/FLCdb/GSEA/NormCounts_Requena.rds")
dea_results <- read.xlsx("~/Documents/SimonLab/FLCdb/DEA/results/panel1/res_T_N_Requena_Panel_1_all.xlsx", colNames = T)

# 3. Pre processing
sampledata <- sampledata[, c("patient_internal", "id", "tumor", "library")]
sampledata <- sampledata[sampledata$library %in% c("RU-A", "RU-B", "RU-C"), ]

paired_samples <- unique(sampledata$patient_internal[duplicated(sampledata$patient_internal)])
sampledata <- sampledata[sampledata$patient_internal %in% paired_samples, ]

normal_samples <- unique(sampledata$patient_internal[sampledata$tumor == "N"])
sampledata <- sampledata[sampledata$patient_internal %in% normal_samples, ]

dea_results <- dea_results[, 1:8]

merged_data <- normCounts[, c(sampledata$id, "Mean.Normal", "Mean.Tumor", "Mean", "geneID", "symbol")]

## Merged data
merged_data <- merged_data %>% left_join(dea_results %>% select(-symbol), by = c("geneID" = "ensembl"))

# 4. Filtering
## FDR, log2FoldChange (from DEA results) and Mean.Tumor (from Normalized Counts)
merged_data <- merged_data %>% filter(padj < 0.05 & Mean.Tumor > 50)

## Fold Change per Patient (from Normalized Counts)
for (patient in sampledata$patient_internal) {
  colname <- paste0("l2FC_RU", patient)
  
  n_counts <- rowMeans(merged_data[sampledata$id[sampledata$patient_internal == patient & sampledata$tumor == "N"]])
  t_counts <- rowMeans(merged_data[sampledata$id[sampledata$patient_internal == patient & sampledata$tumor == "T"]])
  
  n_counts <- ifelse(n_counts == 0, log2(n_counts + 0.01), log2(n_counts))
  t_counts <- ifelse(t_counts == 0, log2(t_counts + 0.01), log2(t_counts))
  
  merged_data[, colname] <- t_counts - n_counts
}

merged_data <- merged_data %>% filter(if_all(starts_with("l2FC_RU"), ~ .x >= 2 & !is.na(.x)))

# 5. Save results
normCounts <- normCounts[normCounts$geneID %in% merged_data$geneID, ]
dea_results <- dea_results[dea_results$ensembl %in% merged_data$geneID, ]

saveRDS(normCounts, "Filtered_Normalized_Counts_IMP.rds")
saveRDS(dea_results, "Filtered_DEA_Results_IMP.rds")
