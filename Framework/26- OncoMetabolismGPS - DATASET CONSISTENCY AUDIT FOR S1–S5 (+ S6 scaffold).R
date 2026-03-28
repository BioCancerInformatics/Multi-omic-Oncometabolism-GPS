###############################################################################
# ONCOMETABOLISMGPS — DATASET CONSISTENCY AUDIT (S1–S5, TSV ONLY)
# PIPELINE-ALIGNED REVISED VERSION
#
# PURPOSE
# - Audit internal consistency across Dataset S1–S5
# - Respect the declared pipeline semantics:
#     S1 = primary significant associations
#     S2 = expanded/annotated multi-mapped associations
#     S3 = reconstructed omic-specific metabolic signatures
#     S4 = regulator–signature interactions
#     S5 = final metabolic regulatory circuitries
#
# IMPORTANT DESIGN CHANGE
# - The following checks are REMOVED as hard failures because they are not
#   valid under the declared pipeline semantics:
#     * S1 row-wise projection represented in S2
#     * S5 sig nomenclatures present in S3
#     * S5 int nomenclatures present in S4
#
# - Instead, the script now performs:
#     * row-count validation
#     * S1 -> S2 domain coverage diagnostics
#     * S1 -> S2 mapping-rate diagnostics
#     * S3 distribution audits for Figure 2 / Results
#     * S4 unique-signature count audit
#     * optional S3 structure diagnostics (non-blocking)
#
# IMPORTANT MOLECULAR-CLASS HARMONIZATION
# - S1 and S2 use different naming granularities for Molecular_class.
# - To avoid false audit failures, Molecular_class is harmonized to a canonical
#   comparison space:
#       Coding gene                   -> Coding-related
#       Enzyme-coding gene            -> Coding-related
#       Coding Transcript Isoform     -> Coding-related
#       Long non-coding Transcript
#         Isoforms                    -> lncRNA
#       lncRNA                        -> lncRNA
#       miRNA                         -> miRNA
#       Protein                       -> Protein
#
# IMPORTANT S3 MEMBER-COUNT AMENDMENT
# - In Dataset S3, the variable 'Members' is the true signature cardinality.
# - Therefore, S3 member-count auditing must use 'Members' directly when it is
#   numeric or numeric-like, rather than tokenizing it as if it were a
#   delimited text field.
# - A robust parser is implemented:
#     * numeric / integer-like Members -> direct cardinality
#     * text-delimited Members         -> token count fallback
#     * invalid / non-positive values  -> hard stop
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

###############################################################################
# 0. WORKING DIRECTORY AND FILES
###############################################################################

rm(list = ls())

setwd("D:/OncoMetabolismGPShiny/Frontiers in Molecular Biosciences/ONCOMETABOLISMGPS — DATASET CONSISTENCY AUDIT FOR S1–S5")
base_dir <- getwd()

dataset_files <- list(
  S1 = file.path(base_dir, "Dataset_S1.tsv"),
  S2 = file.path(base_dir, "Dataset_S2.tsv"),
  S3 = file.path(base_dir, "Dataset_S3.tsv"),
  S4 = file.path(base_dir, "Dataset_S4.tsv"),
  S5 = file.path(base_dir, "Dataset_S5.tsv")
)

###############################################################################
# 1. EXPECTED COUNTS
###############################################################################

expected_counts <- list(
  S1 = 171782L,
  S2 = 463433L,
  S3 = 241415L,
  S4 = 204591L,
  S5 = 24796L,
  S4_unique_signatures = 84270L
)

###############################################################################
# 2. HELPER FUNCTIONS
###############################################################################

read_tsv_strict <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("File not found: ", path), call. = FALSE)
  }
  
  df <- readr::read_delim(
    file = path,
    delim = "\t",
    show_col_types = FALSE,
    progress = FALSE,
    trim_ws = TRUE
  )
  
  df <- as_tibble(df)
  
  names(df) <- names(df) |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_trim()
  
  df
}

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- stringr::str_replace_all(x, "[\u2010\u2011\u2012\u2013\u2014\u2212]", "-")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x <- trimws(x)
  x
}

normalize_numeric <- function(x, digits = 12) {
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- ""
  x_chr <- trimws(x_chr)
  
  num <- suppressWarnings(as.numeric(x_chr))
  
  out <- ifelse(
    x_chr == "",
    "",
    ifelse(
      is.na(num),
      x_chr,
      formatC(signif(num, digits = digits), digits = digits, format = "fg", flag = "#")
    )
  )
  
  out
}

canonical_nomenclature <- function(x) {
  x <- normalize_text(x)
  x <- stringr::str_replace_all(x, "\\s*([+.,;:/|()\\[\\]-])\\s*", "\\1")
  x <- toupper(x)
  x
}

canonical_molecular_class <- function(x) {
  x_raw <- normalize_text(x)
  x_low <- tolower(x_raw)
  
  dplyr::case_when(
    x_low %in% c("coding gene", "enzyme-coding gene", "coding transcript isoform") ~ "Coding-related",
    x_low %in% c("lncrna", "long non-coding transcript isoforms") ~ "lncRNA",
    x_low == "mirna" ~ "miRNA",
    x_low == "protein" ~ "Protein",
    x_low == "" ~ "",
    TRUE ~ x_raw
  )
}

assert_required_columns <- function(df, required, dataset_name) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      paste0(
        dataset_name, " is missing required columns:\n- ",
        paste(missing, collapse = "\n- ")
      ),
      call. = FALSE
    )
  }
}

count_members <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  
  sapply(x, function(xx) {
    if (xx == "") return(0L)
    
    parts <- stringr::str_split(
      xx,
      "\\s*\\+\\s*|\\s*,\\s*|\\s*;\\s*|\\s*\\|\\s*|\\s*/\\s*"
    )[[1]]
    
    parts <- trimws(parts)
    parts <- parts[parts != ""]
    length(parts)
  })
}

parse_member_count <- function(x) {
  # Robust parser for S3 'Members':
  # 1. If x is numeric/integer-like, use it directly as the true cardinality.
  # 2. Otherwise, fall back to token counting for delimited character strings.
  # 3. Reject missing, non-finite, or non-positive counts.
  
  x_chr <- as.character(x)
  x_chr[is.na(x_chr)] <- ""
  x_chr <- trimws(x_chr)
  
  x_num <- suppressWarnings(as.numeric(x_chr))
  
  is_integer_like_numeric <- !is.na(x_num) & is.finite(x_num) & (x_num == floor(x_num))
  is_blank <- x_chr == ""
  
  out <- integer(length(x_chr))
  
  # Direct numeric path
  out[is_integer_like_numeric] <- as.integer(x_num[is_integer_like_numeric])
  
  # Fallback textual token-count path
  fallback_idx <- !is_integer_like_numeric & !is_blank
  if (any(fallback_idx)) {
    out[fallback_idx] <- count_members(x_chr[fallback_idx])
  }
  
  # Preserve blanks / failed parses as NA for explicit validation
  out[is_blank] <- NA_integer_
  
  # Validation
  if (any(is.na(out))) {
    bad_vals <- unique(x_chr[is.na(out)])
    bad_vals <- bad_vals[bad_vals != ""]
    stop(
      paste0(
        "Invalid or missing values detected in S3 'Members'. ",
        "Examples: ", paste(utils::head(bad_vals, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (any(!is.finite(out))) {
    stop("Non-finite values detected in parsed S3 member counts.", call. = FALSE)
  }
  
  if (any(out <= 0L)) {
    bad_vals <- unique(x_chr[out <= 0L])
    stop(
      paste0(
        "Non-positive values detected in S3 'Members'. Examples: ",
        paste(utils::head(bad_vals, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  out
}

has_any_interaction_flag <- function(common_interaction, meaningful_interaction) {
  ci <- tolower(normalize_text(common_interaction))
  mi <- tolower(normalize_text(meaningful_interaction))
  
  false_tokens <- c("", "0", "no", "n", "false", "na", "nan", "null")
  
  ci_yes <- !(ci %in% false_tokens)
  mi_yes <- !(mi %in% false_tokens)
  
  ci_yes | mi_yes
}

safe_count_table <- function(df, column_name) {
  if (!column_name %in% names(df)) {
    return(tibble())
  }
  
  df %>%
    mutate(.tmp = normalize_text(.data[[column_name]])) %>%
    count(.tmp, sort = TRUE, name = "n") %>%
    rename(value = .tmp) %>%
    mutate(pct = 100 * n / sum(n))
}

build_s1s2_components <- function(df) {
  df %>%
    transmute(
      Target_key                    = normalize_text(Target),
      Molecular_class_key           = canonical_molecular_class(Molecular_class),
      Molecular_class_raw_key       = normalize_text(Molecular_class),
      Cancer_types_key              = normalize_text(Cancer_types),
      Omic_layer_key                = normalize_text(Omic_layer),
      Phenotypic_layer_key          = normalize_text(Phenotypic_layer),
      Correlation_rho_key           = normalize_numeric(Correlation_rho, digits = 12),
      Correlation_padj_key          = normalize_numeric(`Correlation_p.adj`, digits = 12)
    )
}

make_key <- function(df, cols) {
  apply(df[, cols, drop = FALSE], 1, paste, collapse = "||")
}

###############################################################################
# 3. FILE EXISTENCE CHECK
###############################################################################

file_check <- vapply(dataset_files, file.exists, logical(1))

cat("File existence check:\n")
print(file_check)
cat("\n")

if (!all(file_check)) {
  stop("One or more TSV files were not found. Check file names and paths.", call. = FALSE)
}

###############################################################################
# 4. LOAD DATASETS
###############################################################################

dataset_s1 <- read_tsv_strict(dataset_files$S1)
dataset_s2 <- read_tsv_strict(dataset_files$S2)
dataset_s3 <- read_tsv_strict(dataset_files$S3)
dataset_s4 <- read_tsv_strict(dataset_files$S4)
dataset_s5 <- read_tsv_strict(dataset_files$S5)

cat("Loaded datasets:\n")
cat("S1:", nrow(dataset_s1), "rows\n")
cat("S2:", nrow(dataset_s2), "rows\n")
cat("S3:", nrow(dataset_s3), "rows\n")
cat("S4:", nrow(dataset_s4), "rows\n")
cat("S5:", nrow(dataset_s5), "rows\n\n")

###############################################################################
# 5. REQUIRED COLUMNS
###############################################################################

required_s1 <- c(
  "Target", "Molecular_class", "Cancer_types", "Omic_layer",
  "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj"
)

required_s2 <- c(
  "Target", "Molecular_class", "Cancer_types", "Metabolism", "Pathways",
  "Metabolic_cell_death", "Omic_layer", "Phenotypic_layer",
  "Correlation_rho", "Correlation_p.adj"
)

required_s3 <- c(
  "Members", "Signatures", "Nomenclature", "Molecular_class",
  "Common_interaction", "Meaningful_interaction", "Metabolism",
  "Pathways", "Metabolic_cell_death", "CTAB", "Omic_layer",
  "Phenotypic_layer"
)

required_s4 <- c(
  "Nomenclature", "Signatures", "Meaningful_interaction", "CTAB",
  "Metabolism", "Pathways", "Metabolic_cell_death",
  "Omic_layer_signature", "Phenotypic_layer_signature",
  "Omic_layer_interaction", "Phenotypic_layer_interaction",
  "Final_concordance_summary"
)

required_s5 <- c(
  "Circuitries_id", "Nomenclature_sig", "Signatures", "Nomenclature_int",
  "Interaction", "CTAB", "Metabolism", "Pathways", "Metabolic_cell_death",
  "Omic_layer_sig", "Phenotypic_layer_sig", "Omic_layer_int",
  "Phenotypic_layer_int", "Final_concordance_summary"
)

assert_required_columns(dataset_s1, required_s1, "Dataset S1")
assert_required_columns(dataset_s2, required_s2, "Dataset S2")
assert_required_columns(dataset_s3, required_s3, "Dataset S3")
assert_required_columns(dataset_s4, required_s4, "Dataset S4")
assert_required_columns(dataset_s5, required_s5, "Dataset S5")

###############################################################################
# 6. BASIC ROW COUNT AUDIT
###############################################################################

audit_counts <- tibble(
  Dataset = c("S1", "S2", "S3", "S4", "S5"),
  Expected = c(
    expected_counts$S1,
    expected_counts$S2,
    expected_counts$S3,
    expected_counts$S4,
    expected_counts$S5
  ),
  Observed = c(
    nrow(dataset_s1),
    nrow(dataset_s2),
    nrow(dataset_s3),
    nrow(dataset_s4),
    nrow(dataset_s5)
  )
) %>%
  mutate(Pass = Expected == Observed)

cat("Basic row count audit:\n")
print(audit_counts)
cat("\n")

###############################################################################
# 7. S1 -> S2 DOMAIN COVERAGE AND MAPPING DIAGNOSTICS
#
# NOTE
# This section is diagnostic only. It is NOT a hard-fail row-wise projection
# audit because S2 is an expanded multi-mapped annotation layer, not a 1:1
# projection of S1.
###############################################################################

s1_comp <- build_s1s2_components(dataset_s1)
s2_comp <- build_s1s2_components(dataset_s2)

# Full key-based diagnostic mapping rate
# Uses harmonized Molecular_class_key to avoid false mismatch driven by label
# granularity differences between S1 and S2.
s1_full <- s1_comp %>%
  select(
    Target_key, Molecular_class_key, Cancer_types_key,
    Omic_layer_key, Phenotypic_layer_key,
    Correlation_rho_key, Correlation_padj_key
  ) %>%
  mutate(full_key = make_key(., names(.)))

s2_full <- s2_comp %>%
  select(
    Target_key, Molecular_class_key, Cancer_types_key,
    Omic_layer_key, Phenotypic_layer_key,
    Correlation_rho_key, Correlation_padj_key
  ) %>%
  mutate(full_key = make_key(., names(.)))

missing_in_s2 <- anti_join(
  s1_full,
  s2_full %>% select(full_key),
  by = "full_key"
)

s1_to_s2_mapping_rate <- 1 - (nrow(missing_in_s2) / nrow(dataset_s1))

cat("S1 -> S2 mapping diagnostics:\n")
cat("  S1 rows not found in S2 under full key: ", nrow(missing_in_s2), "\n", sep = "")
cat("  S1 -> S2 mapping rate under full key: ", round(100 * s1_to_s2_mapping_rate, 3), "%\n\n", sep = "")

write_tsv(missing_in_s2, file.path(base_dir, "audit_missing_S1_projection_in_S2.tsv"))

# Relaxed-key diagnostics
diagnostic_key_sets <- list(
  full_key = c(
    "Target_key", "Molecular_class_key", "Cancer_types_key",
    "Omic_layer_key", "Phenotypic_layer_key",
    "Correlation_rho_key", "Correlation_padj_key"
  ),
  no_numeric = c(
    "Target_key", "Molecular_class_key", "Cancer_types_key",
    "Omic_layer_key", "Phenotypic_layer_key"
  ),
  no_target = c(
    "Molecular_class_key", "Cancer_types_key", "Omic_layer_key",
    "Phenotypic_layer_key", "Correlation_rho_key", "Correlation_padj_key"
  ),
  target_class_cancer_only = c(
    "Target_key", "Molecular_class_key", "Cancer_types_key"
  )
)

s1_s2_relaxed_diag <- imap_dfr(diagnostic_key_sets, function(cols, nm) {
  s1_tmp <- s1_comp %>%
    mutate(tmp_key = make_key(., cols)) %>%
    distinct(tmp_key)
  
  s2_tmp <- s2_comp %>%
    mutate(tmp_key = make_key(., cols)) %>%
    distinct(tmp_key)
  
  tibble(
    Diagnostic = nm,
    Missing_rows = anti_join(s1_tmp, s2_tmp, by = "tmp_key") %>% nrow()
  )
})

cat("S1 -> S2 relaxed-key diagnostics:\n")
print(s1_s2_relaxed_diag)
cat("\n")

write_tsv(s1_s2_relaxed_diag, file.path(base_dir, "audit_S1_S2_relaxed_key_diagnostics.tsv"))

# Fieldwise domain coverage
field_coverage <- tibble(
  Field = c("Target", "Molecular_class", "Cancer_types", "Omic_layer", "Phenotypic_layer"),
  S1_unique = c(
    n_distinct(s1_comp$Target_key),
    n_distinct(s1_comp$Molecular_class_key),
    n_distinct(s1_comp$Cancer_types_key),
    n_distinct(s1_comp$Omic_layer_key),
    n_distinct(s1_comp$Phenotypic_layer_key)
  ),
  Values_in_S1_not_in_S2 = c(
    sum(!(unique(s1_comp$Target_key) %in% unique(s2_comp$Target_key))),
    sum(!(unique(s1_comp$Molecular_class_key) %in% unique(s2_comp$Molecular_class_key))),
    sum(!(unique(s1_comp$Cancer_types_key) %in% unique(s2_comp$Cancer_types_key))),
    sum(!(unique(s1_comp$Omic_layer_key) %in% unique(s2_comp$Omic_layer_key))),
    sum(!(unique(s1_comp$Phenotypic_layer_key) %in% unique(s2_comp$Phenotypic_layer_key)))
  )
)

cat("S1 -> S2 field value coverage diagnostics:\n")
print(field_coverage)
cat("\n")

write_tsv(field_coverage, file.path(base_dir, "audit_S1_S2_field_value_coverage.tsv"))

# Harmonized hard check retained only for domain consistency
molecular_class_missing <- setdiff(
  unique(canonical_molecular_class(dataset_s1$Molecular_class)),
  unique(canonical_molecular_class(dataset_s2$Molecular_class))
)

molecular_class_check <- tibble(
  Check = "S1 Molecular_class levels absent from S2 (harmonized)",
  Expected = 0L,
  Observed = length(molecular_class_missing)
) %>%
  mutate(Pass = Expected == Observed)

if (length(molecular_class_missing) > 0) {
  write_lines(
    molecular_class_missing,
    file.path(base_dir, "audit_missing_molecular_class_levels_S1_not_in_S2_harmonized.txt")
  )
}

# Optional transparency table: raw-to-canonical mapping observed in each dataset
molecular_class_harmonization_table <- bind_rows(
  tibble(
    Dataset = "S1",
    Molecular_class_raw = sort(unique(normalize_text(dataset_s1$Molecular_class)))
  ),
  tibble(
    Dataset = "S2",
    Molecular_class_raw = sort(unique(normalize_text(dataset_s2$Molecular_class)))
  )
) %>%
  mutate(Molecular_class_harmonized = canonical_molecular_class(Molecular_class_raw))

write_tsv(
  molecular_class_harmonization_table,
  file.path(base_dir, "audit_molecular_class_harmonization_table.tsv")
)

###############################################################################
# 8. S3 ANALYSIS
#
# NOTE
# The member/interactions logic is retained as a NON-BLOCKING diagnostic.
# However, the member-count logic is now corrected and validated:
# - If 'Members' is numeric or numeric-like, it is used directly.
# - If 'Members' is a delimited text field in a future dataset version,
#   token-count fallback is used.
###############################################################################

dataset_s3 <- dataset_s3 %>%
  mutate(
    member_count      = parse_member_count(Members),
    multi_member      = member_count > 1L,
    has_interaction   = has_any_interaction_flag(Common_interaction, Meaningful_interaction),
    Nomenclature_norm = canonical_nomenclature(Nomenclature),
    Omic_layer_norm   = normalize_text(Omic_layer),
    Metabolism_norm   = normalize_text(Metabolism),
    MCD_norm          = normalize_text(Metabolic_cell_death)
  )

# Explicit validation of S3 member counts after parsing
if (any(is.na(dataset_s3$member_count))) {
  stop("NA values detected in S3 member_count after parsing.", call. = FALSE)
}

if (any(dataset_s3$member_count <= 0L)) {
  stop("Non-positive values detected in S3 member_count after parsing.", call. = FALSE)
}

s3_total  <- nrow(dataset_s3)
s3_single <- sum(!dataset_s3$multi_member, na.rm = TRUE)
s3_multi  <- sum(dataset_s3$multi_member, na.rm = TRUE)

s3_with_interaction    <- sum(dataset_s3$has_interaction, na.rm = TRUE)
s3_without_interaction <- s3_total - s3_with_interaction

cat("S3 structural diagnostics (non-blocking):\n")
cat("  Total signatures: ", s3_total, "\n", sep = "")
cat("  Single-member signatures: ", s3_single, "\n", sep = "")
cat("  Multi-member signatures: ", s3_multi, "\n", sep = "")
cat("  With interaction: ", s3_with_interaction, "\n", sep = "")
cat("  Without interaction: ", s3_without_interaction, "\n", sep = "")
cat("  Member-count range: ", min(dataset_s3$member_count), " to ", max(dataset_s3$member_count), "\n\n", sep = "")

s3_structure_diagnostics <- tibble(
  Metric = c(
    "Total signatures",
    "Single-member signatures",
    "Multi-member signatures",
    "With interaction",
    "Without interaction",
    "Minimum member count",
    "Maximum member count"
  ),
  Value = c(
    s3_total,
    s3_single,
    s3_multi,
    s3_with_interaction,
    s3_without_interaction,
    min(dataset_s3$member_count),
    max(dataset_s3$member_count)
  )
)

write_tsv(
  s3_structure_diagnostics,
  file.path(base_dir, "audit_S3_structure_diagnostics.tsv")
)

# Optional detailed S3 member-count distribution
s3_member_count_distribution <- dataset_s3 %>%
  count(member_count, name = "n") %>%
  arrange(member_count) %>%
  mutate(pct = 100 * n / sum(n))

write_tsv(
  s3_member_count_distribution,
  file.path(base_dir, "audit_S3_member_count_distribution.tsv")
)

###############################################################################
# 9. S4 ANALYSIS
###############################################################################

dataset_s4 <- dataset_s4 %>%
  mutate(
    Nomenclature_norm           = canonical_nomenclature(Nomenclature),
    Signatures_norm             = normalize_text(Signatures),
    Meaningful_interaction_norm = canonical_nomenclature(Meaningful_interaction)
  )

s4_unique_signatures <- dataset_s4 %>%
  distinct(Nomenclature_norm) %>%
  filter(Nomenclature_norm != "") %>%
  nrow()

cat("S4 summary:\n")
cat("  Total interactions: ", nrow(dataset_s4), "\n", sep = "")
cat("  Unique signatures with interaction: ", s4_unique_signatures, "\n\n", sep = "")

###############################################################################
# 10. S5 DESCRIPTIVE DIAGNOSTICS ONLY
#
# NOTE
# No hard check is imposed on direct nomenclature containment because S5 is a
# circuitry atlas derived by mapping S4 composite interaction signatures onto S3.
###############################################################################

dataset_s5 <- dataset_s5 %>%
  mutate(
    Nomenclature_sig_norm = canonical_nomenclature(Nomenclature_sig),
    Nomenclature_int_norm = canonical_nomenclature(Nomenclature_int),
    Interaction_norm      = canonical_nomenclature(Interaction),
    Signatures_norm       = normalize_text(Signatures)
  )

blank_sig <- dataset_s5 %>%
  filter(Nomenclature_sig_norm == "") %>%
  nrow()

blank_int <- dataset_s5 %>%
  filter(Nomenclature_int_norm == "") %>%
  nrow()

cat("S5 descriptive diagnostics:\n")
cat("  Blank Nomenclature_sig: ", blank_sig, "\n", sep = "")
cat("  Blank Nomenclature_int: ", blank_int, "\n\n", sep = "")

s5_blank_diag <- tibble(
  Check = c(
    "S5 blank sig nomenclatures",
    "S5 blank int nomenclatures"
  ),
  Expected = c(0L, 0L),
  Observed = c(blank_sig, blank_int)
) %>%
  mutate(Pass = Expected == Observed)

###############################################################################
# 11. DISTRIBUTION TABLES FOR MANUSCRIPT / FIGURE 2 CHECKS
###############################################################################

dist_s3_omic <- safe_count_table(dataset_s3, "Omic_layer_norm")
dist_s3_path <- safe_count_table(dataset_s3, "Metabolism_norm")
dist_s3_mcd  <- safe_count_table(dataset_s3, "MCD_norm")

write_csv(dist_s3_omic, file.path(base_dir, "audit_S3_distribution_by_omic.csv"))
write_csv(dist_s3_path, file.path(base_dir, "audit_S3_distribution_by_pathway.csv"))
write_csv(dist_s3_mcd,  file.path(base_dir, "audit_S3_distribution_by_mcd.csv"))

cat("S3 distribution by omic layer:\n")
print(dist_s3_omic)
cat("\n")

cat("S3 distribution by metabolic category:\n")
print(dist_s3_path)
cat("\n")

cat("S3 distribution by metabolic cell death:\n")
print(dist_s3_mcd)
cat("\n")

###############################################################################
# 12. FINAL VALID CONSISTENCY CHECK TABLE
###############################################################################

audit_checks <- bind_rows(
  tibble(
    Check = c(
      "S1 row count",
      "S2 row count",
      "S3 row count",
      "S4 row count",
      "S5 row count",
      "S4 unique signatures count"
    ),
    Expected = c(
      expected_counts$S1,
      expected_counts$S2,
      expected_counts$S3,
      expected_counts$S4,
      expected_counts$S5,
      expected_counts$S4_unique_signatures
    ),
    Observed = c(
      nrow(dataset_s1),
      nrow(dataset_s2),
      nrow(dataset_s3),
      nrow(dataset_s4),
      nrow(dataset_s5),
      s4_unique_signatures
    )
  ) %>% mutate(Pass = Expected == Observed),
  molecular_class_check,
  s5_blank_diag
)

cat("Final valid audit check table:\n")
print(audit_checks)
cat("\n")

write_tsv(audit_checks, file.path(base_dir, "audit_report_S1_S5.tsv"))

###############################################################################
# 13. HARD STOP ONLY FOR VALID CHECKS
###############################################################################

if (any(!audit_checks$Pass)) {
  failed <- audit_checks %>% filter(!Pass)
  
  cat("FAILED VALID CHECKS:\n")
  print(failed)
  cat("\n")
  
  stop("One or more valid audit checks failed. Review the failed checks above.", call. = FALSE)
}

cat("ALL VALID AUDITS PASSED — DATASET CONSISTENCY VERIFIED\n")

###############################################################################
# 14. OPTIONAL INTERPRETIVE SUMMARY
###############################################################################

audit_summary <- tibble(
  Item = c(
    "S1_to_S2_full_key_missing_rows",
    "S1_to_S2_full_key_mapping_rate",
    "S3_total_signatures",
    "S3_single_member_signatures",
    "S3_multi_member_signatures",
    "S3_with_interaction",
    "S3_without_interaction",
    "S3_member_count_min",
    "S3_member_count_max",
    "S4_unique_signatures",
    "S5_blank_sig_nomenclatures",
    "S5_blank_int_nomenclatures"
  ),
  Value = c(
    nrow(missing_in_s2),
    round(s1_to_s2_mapping_rate, 6),
    s3_total,
    s3_single,
    s3_multi,
    s3_with_interaction,
    s3_without_interaction,
    min(dataset_s3$member_count),
    max(dataset_s3$member_count),
    s4_unique_signatures,
    blank_sig,
    blank_int
  )
)

write_tsv(audit_summary, file.path(base_dir, "audit_summary_interpretive.tsv"))

###############################################################################
# 15. OPTIONAL CONSOLE DIAGNOSTICS
###############################################################################

cat("Raw Molecular_class levels in S1:\n")
print(sort(unique(normalize_text(dataset_s1$Molecular_class))))
cat("\n")

cat("Raw Molecular_class levels in S2:\n")
print(sort(unique(normalize_text(dataset_s2$Molecular_class))))
cat("\n")

cat("Harmonized Molecular_class levels in S1:\n")
print(sort(unique(canonical_molecular_class(dataset_s1$Molecular_class))))
cat("\n")

cat("Harmonized Molecular_class levels in S2:\n")
print(sort(unique(canonical_molecular_class(dataset_s2$Molecular_class))))
cat("\n")

cat("S3 member-count summary:\n")
print(summary(dataset_s3$member_count))
cat("\n")

cat("S3 member-count range:\n")
print(range(dataset_s3$member_count))
cat("\n")

setdiff(
  unique(canonical_molecular_class(dataset_s1$Molecular_class)),
  unique(canonical_molecular_class(dataset_s2$Molecular_class))
)

###############################################################################
# 16. OPTIONAL RE-ESTIMATION OF MCD DISTRIBUTION EXCLUDING 'Unrelated'
###############################################################################

df_mcd <- read_csv(
  file.path(base_dir, "audit_S3_distribution_by_mcd.csv"),
  show_col_types = FALSE
)

df_mcd_excluding_unrelated <- df_mcd %>%
  filter(value != "Unrelated") %>%
  mutate(pct = (n / sum(n)) * 100)

write_csv(
  df_mcd_excluding_unrelated,
  file.path(base_dir, "audit_S3_distribution_by_mcd_excluding_unrelated.csv")
)

cat("ALL AUDITS PASSED — DATASET CONSISTENCY VERIFIED\n")

#### 
#### 
#### plot_audit_results.R
# plot_audit_results.R
# This script reads the dataset consistency audit results and generates comparative plots.

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# Define the base directory (where the script is located)
base_dir <- getwd()

# Create output directory for plots
output_dir <- file.path(base_dir, "audit_plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Helper function to plot percentage bar charts
plot_percentage_bar <- function(data_file, title, output_filename, color_fill) {
  if (file.exists(data_file)) {
    # Read data
    df <- read_csv(data_file, show_col_types = FALSE)
    
    # Ensure columns exist
    if(all(c("value", "n", "pct") %in% colnames(df))) {
      
      # Order factors by percentage descending
      df <- df %>% 
        arrange(desc(pct)) %>% 
        mutate(value = factor(value, levels = value)) %>%
        mutate(label_text = case_when(
          pct > 0 & pct < 0.1 ~ "<0.1%",
          TRUE ~ sprintf("%.1f%%", pct)
        ))
      
      # Create replacing underscores with spaces for readability in labels
      df$label_value <- str_replace_all(df$value, "_", " ")
      df$label_value <- factor(df$label_value, levels = df$label_value)
      
      # Plot
      p <- ggplot(df, aes(x = label_value, y = pct)) +
        geom_bar(stat = "identity", fill = color_fill) +
        geom_text(aes(label = label_text), vjust = -0.5, size = 3.5) +
        theme_minimal(base_size = 14) +
        labs(title = title,
             x = "",
             y = "Percentage (%)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.margin = margin(t = 20, r = 20, b = 20, l = 20))
      
      # Save plot
      output_path <- file.path(output_dir, output_filename)
      ggsave(output_path, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
      cat("Saved plot to:", output_path, "\n")
      
      return(p)
    } else {
      cat("Columns missing in", data_file, "- expected value, n, pct.\n")
      return(NULL)
    }
  } else {
    cat("File not found:", data_file, "\n")
    return(NULL)
  }
}

# Plot files
cat("Generating individual plots...\n")

# 1. MCD Distribution
file_mcd <- file.path(base_dir, "audit_S3_distribution_by_mcd.csv")
p1 <- plot_percentage_bar(
  file_mcd, 
  "Distribution by Molecular Cell Death (MCD)", 
  "distribution_by_mcd.png", 
  "#3498db"
)

# 2. MCD Distribution (Excluding Unrelated)
file_mcd_ex <- file.path(base_dir, "audit_S3_distribution_by_mcd_excluding_unrelated.csv")
p2 <- plot_percentage_bar(
  file_mcd_ex, 
  "Distribution by MCD (Excluding 'Unrelated')", 
  "distribution_by_mcd_excluding_unrelated.png", 
  "#e74c3c"
)

# 3. Omic Distribution
file_omic <- file.path(base_dir, "audit_S3_distribution_by_omic.csv")
p3 <- plot_percentage_bar(
  file_omic, 
  "Distribution by Omic Strategy", 
  "distribution_by_omic.png", 
  "#2ecc71"
)

# 4. Pathway Distribution
file_pathway <- file.path(base_dir, "audit_S3_distribution_by_pathway.csv")
p4 <- plot_percentage_bar(
  file_pathway, 
  "Distribution by Metabolic Pathway", 
  "distribution_by_pathway.png", 
  "#f39c12"
)

cat("Generating combined multi-panel plot...\n")
# Combine plots into a single figure with panels A, B, C, D
if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
  
  # Try to use patchwork (recommended for ggplot2 multi-panel plots)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    suppressPackageStartupMessages(library(patchwork))
    
    # Arrange 2x2 with tags A, B, C, D
    combined_plot <- (p1 | p2) / (p3 | p4) + 
      plot_annotation(tag_levels = 'A')
    
    combined_path <- file.path(output_dir, "combined_audit_plots.png")
    ggsave(combined_path, plot = combined_plot, width = 20, height = 14, dpi = 300, bg = "white")
    cat("Saved combined plot (using patchwork) to:", combined_path, "\n")
    
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    # Fallback to cowplot if patchwork is not installed
    suppressPackageStartupMessages(library(cowplot))
    
    combined_plot <- plot_grid(p1, p2, p3, p4, labels = "AUTO", ncol = 2)
    
    combined_path <- file.path(output_dir, "combined_audit_plots.png")
    ggsave(combined_path, plot = combined_plot, width = 20, height = 14, dpi = 300, bg = "white")
    cat("Saved combined plot (using cowplot) to:", combined_path, "\n")
  } else {
    cat("\n*** NOTE: To generate the combined 4-panel figure, please install the 'patchwork' package.\n")
    cat("You can install it by running: install.packages('patchwork') in your R console.\n")
  }
}

cat("Successfully finished script execution.\n")

