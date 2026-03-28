library(data.table)

cat("--- S3 Audit ---\n")
s3 <- fread("Dataset_S3.tsv")
cat("Total signatures:", nrow(s3), "\n")
cat("Signatures with Common_interaction not empty:\n")
print(table(s3$Common_interaction != "" & s3$Common_interaction != "None" & !is.na(s3$Common_interaction)))
cat("Signatures with Meaningful_interaction not empty:\n")
print(table(s3$Meaningful_interaction != "" & s3$Meaningful_interaction != "None" & !is.na(s3$Meaningful_interaction)))
cat("Multi-member breakdown:\n")
s3[, is_multi := grepl("\\+", Signatures)]
cat("Is multi:\n")
print(table(s3$is_multi))
cat("Multi with interactions:\n")
print(table(s3[is_multi == TRUE]$Common_interaction != "" & s3[is_multi == TRUE]$Common_interaction != "None"))
cat("Single with interactions:\n")
print(table(s3[is_multi == FALSE]$Common_interaction != "" & s3[is_multi == FALSE]$Common_interaction != "None"))

cat("\n--- S4 Audit ---\n")
s4 <- fread("Dataset_S4.tsv")
cat("Total interactions analyzed for convergence/divergence:", nrow(s4), "\n")
cat("Unique signatures in S4:", length(unique(s4$Signatures)), "\n")
cat("Unique Nomenclature in S4:", length(unique(s4$Nomenclature)), "\n")

cat("\n--- S5 Audit ---\n")
s5 <- fread("Dataset_S5.tsv")
cat("Total circuitries in S5:", nrow(s5), "\n")
if("Category" %in% names(s5)) {
    cat("Category breakdown:\n")
    print(table(s5$Category))
    print(prop.table(table(s5$Category)) * 100)
} else {
    cat("Signature_type / Interaction_type breakdown:\n")
    print(table(s5$Signature_type, s5$Interaction_type))
}

cat("\n--- S2 Audit ---\n")
s2 <- fread("Dataset_S2.tsv")
cat("Total associations in S2:", nrow(s2), "\n")
s2_prog <- s2[Cox_OS_type %in% c("Protective", "Risky") | 
              Cox_DSS_type %in% c("Protective", "Risky") | 
              Cox_DFI_type %in% c("Protective", "Risky") | 
              Cox_PFI_type %in% c("Protective", "Risky") | 
              OS_worst_prognosis_group %in% c("High", "Low", "Mutant worst", "Wild-type worst", "Deleted worst", "Duplicated/Normal worst", "Mixed/Compound") | 
              DSS_worst_prognosis_group %in% c("High", "Low", "Mutant worst", "Wild-type worst", "Deleted worst", "Duplicated/Normal worst", "Mixed/Compound") | 
              DFI_worst_prognosis_group %in% c("High", "Low", "Mutant worst", "Wild-type worst", "Deleted worst", "Duplicated/Normal worst", "Mixed/Compound") | 
              PFI_worst_prognosis_group %in% c("High", "Low", "Mutant worst", "Wild-type worst", "Deleted worst", "Duplicated/Normal worst", "Mixed/Compound")]

cat("Total prognostic associations in S2:", nrow(s2_prog), "\n")
cat("Immune classification breakdown (counts):\n")
print(table(s2_prog$Immune_classification))
cat("Immune classification breakdown (percentages):\n")
print(prop.table(table(s2_prog$Immune_classification)) * 100)
