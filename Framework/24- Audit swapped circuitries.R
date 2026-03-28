###############################################################################
# AUDIT + NORMALIZE + DEDUPLICATE Circuitries_id
# Goal:
#  (1) Strictly audit whether Circuitries_id follows EXACT canonical format:
#      "CANCER-NNNN / CANCER-NNNN"  (ONE space, slash, ONE space)
#  (2) Detect whitespace/pathology variants (leading/trailing, multi-spaces, etc.)
#  (3) Detect redundant circuitries due to swapped order (A / B vs B / A)
#  (4) Provide an EXACT-MATCH filter that returns ONLY the 4 requested rows
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(rio)
})

circuitries <-  import("Dataset_S5.tsv") # original dataset

stopifnot(exists("circuitries"))
stopifnot(is.data.frame(circuitries))
stopifnot("Circuitries_id" %in% names(circuitries))

# -----------------------------------------------------------------------------
# 1) Define canonical strict pattern
#    EXACTLY: token-digits␠/␠token-digits
# -----------------------------------------------------------------------------
canon_pat <- "^[A-Za-z0-9]+-[0-9]+ / [A-Za-z0-9]+-[0-9]+$"

# -----------------------------------------------------------------------------
# 2) Core helpers
# -----------------------------------------------------------------------------
trim_ws <- function(x) str_trim(as.character(x), side = "both")

# Collapse any whitespace into single spaces, then enforce " / " around slash
# This is a *normalizer*; use only if you want to repair IDs, not just audit.
normalize_circuitry_id_strict <- function(x) {
  x <- trim_ws(x)
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_replace_all(x, "\\s*/\\s*", " / ")
  x
}

# Order-invariant key so "A / B" and "B / A" map to the same identity
# (Assumes normalized delimiter " / " is present)
make_order_invariant_key <- function(x) {
  x <- normalize_circuitry_id_strict(x)
  parts <- str_split_fixed(x, " / ", 2)
  # If split fails, return as-is (will be flagged elsewhere)
  if (ncol(parts) < 2) return(x)
  a <- parts[, 1]
  b <- parts[, 2]
  ifelse(a <= b, paste0(a, " / ", b), paste0(b, " / ", a))
}

# -----------------------------------------------------------------------------
# 3) AUDIT TABLE: identify format violations and whitespace anomalies
# -----------------------------------------------------------------------------
audit_tbl <- circuitries %>%
  transmute(
    Circuitries_id_raw = as.character(.data$Circuitries_id),
    Circuitries_id_trim = trim_ws(.data$Circuitries_id),
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    
    # STRICT format checks
    is_canonical_raw  = str_detect(.data$Circuitries_id_raw,  canon_pat),
    is_canonical_trim = str_detect(.data$Circuitries_id_trim, canon_pat),
    is_canonical_norm = str_detect(.data$Circuitries_id_norm, canon_pat),
    
    # Common pathology flags
    has_leading_or_trailing_ws = (.data$Circuitries_id_raw != .data$Circuitries_id_trim),
    has_any_tab_or_newline     = str_detect(.data$Circuitries_id_raw, "[\\t\\r\\n]"),
    has_multiple_spaces        = str_detect(.data$Circuitries_id_raw, " {2,}"),
    has_no_spaces_around_slash = str_detect(.data$Circuitries_id_raw, "/") &
      !str_detect(.data$Circuitries_id_raw, " / "),
    has_weird_slash_spacing    = str_detect(.data$Circuitries_id_raw, "\\s*/\\s*") &
      !str_detect(.data$Circuitries_id_raw, " / ")
  )

# Summary counts
audit_summary <- audit_tbl %>%
  summarise(
    n_total = n(),
    n_canonical_raw  = sum(is_canonical_raw,  na.rm = TRUE),
    n_canonical_trim = sum(is_canonical_trim, na.rm = TRUE),
    n_canonical_norm = sum(is_canonical_norm, na.rm = TRUE),
    
    n_leading_trailing_ws = sum(has_leading_or_trailing_ws, na.rm = TRUE),
    n_tabs_newlines       = sum(has_any_tab_or_newline,     na.rm = TRUE),
    n_multi_spaces        = sum(has_multiple_spaces,        na.rm = TRUE),
    n_no_spaces_around_sl = sum(has_no_spaces_around_slash, na.rm = TRUE),
    n_weird_slash_spacing = sum(has_weird_slash_spacing,    na.rm = TRUE)
  )

cat("\n================ Circuitries_id AUDIT SUMMARY ================\n")
print(audit_summary)

# Show a small set of non-canonical raw examples (if any)
noncanon_examples <- audit_tbl %>%
  filter(!is_canonical_raw) %>%
  distinct(Circuitries_id_raw) %>%
  slice_head(n = 20)

if (nrow(noncanon_examples) > 0) {
  cat("\n--- Examples of NON-canonical Circuitries_id (raw) ---\n")
  print(noncanon_examples, row.names = FALSE)
} else {
  cat("\nAll Circuitries_id values are canonical under the STRICT raw pattern.\n")
}

# -----------------------------------------------------------------------------
# 4) REDUNDANCY CHECK: swapped-order duplicates
# -----------------------------------------------------------------------------
dup_tbl <- circuitries %>%
  mutate(
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    circuitry_key = make_order_invariant_key(.data$Circuitries_id)
  ) %>%
  count(circuitry_key, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows))

cat("\n================ Swapped-order DUPLICATES (key-level) ================\n")
if (nrow(dup_tbl) == 0) {
  cat("No swapped-order duplicates detected (order-invariant keys are unique).\n")
} else {
  print(dup_tbl, row.names = FALSE)
}

# OPTIONAL: inspect the first duplicated key’s actual raw rows
if (nrow(dup_tbl) > 0) {
  first_key <- dup_tbl$circuitry_key[1]
  cat("\n--- Raw rows for first duplicated key ---\n")
  print(
    circuitries %>%
      mutate(
        Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
        circuitry_key = make_order_invariant_key(.data$Circuitries_id)
      ) %>%
      filter(circuitry_key == first_key) %>%
      select(Circuitries_id, Circuitries_id_norm, circuitry_key) %>%
      distinct(),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 5) EXACT-MATCH FILTER (STRICT): returns ONLY rows whose *raw* string matches
#    exactly one of your provided IDs (including spaces).
#    This will NOT match swapped order.
# -----------------------------------------------------------------------------
target_ids_exact <- c(
  "LGG-5422 / LGG-6519",
  "LGG-5904 / LGG-6708",
  "LUSC-6681 / LUSC-8604",
  "LUAD-10062 / LUAD-11695"
)

circuitries_exact4 <- circuitries %>%
  filter(.data$Circuitries_id %in% target_ids_exact)

cat("\n================ EXACT-MATCH FILTER RESULT ================\n")
cat("Rows returned: ", nrow(circuitries_exact4), "\n", sep = "")
print(circuitries_exact4 %>% select(Circuitries_id) %>% distinct(), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6) ORDER-INVARIANT FILTER (OPTIONAL): matches either orientation by key.
#    Use this if you want to retrieve rows regardless of "A / B" vs "B / A".
# -----------------------------------------------------------------------------
target_keys <- make_order_invariant_key(target_ids_exact)

circuitries_by_key <- circuitries %>%
  mutate(circuitry_key = make_order_invariant_key(.data$Circuitries_id)) %>%
  filter(circuitry_key %in% target_keys)

cat("\n================ ORDER-INVARIANT FILTER RESULT (optional) ================\n")
cat("Rows returned: ", nrow(circuitries_by_key), "\n", sep = "")
print(
  circuitries_by_key %>%
    select(Circuitries_id, circuitry_key) %>%
    distinct() %>%
    arrange(circuitry_key, Circuitries_id),
  row.names = FALSE
)
#####
##### Circuitries_id redundancy due to swapped component ordering: 
##### bidirectional assignment during interaction–signature pairing generates duplicated circuitries.
##### 

# --- helpers assumed to exist (as in your previous code) ----------------------
# normalize_circuitry_id_strict(x)  -> enforces strict "A / B" formatting
# make_order_invariant_key(x)       -> returns canonical key "min(A,B) / max(A,B)"

# -----------------------------------------------------------------------------
# 4) REDUNDANCY CHECK: swapped-order duplicates circuitries (KEY-LEVEL)
#     - dup_tbl now INCLUDES the swapped/raw Circuitries_id variants per key
# -----------------------------------------------------------------------------
library(rio)
library(dplyr)
circuitries_original <- import("Dataset_S5.tsv")

dup_tbl <- circuitries_original %>%
  dplyr::mutate(
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    circuitry_key       = make_order_invariant_key(.data$Circuitries_id)
  ) %>%
  dplyr::group_by(.data$circuitry_key) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    
    # how many DISTINCT raw orientations exist for this key
    n_distinct_raw  = dplyr::n_distinct(.data$Circuitries_id),
    n_distinct_norm = dplyr::n_distinct(.data$Circuitries_id_norm),
    
    # attach the actual raw and normalized variants (this is what you asked for)
    raw_ids  = list(sort(unique(.data$Circuitries_id))),
    norm_ids = list(sort(unique(.data$Circuitries_id_norm))),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(.data$n_rows > 1) %>%
  dplyr::arrange(dplyr::desc(.data$n_rows), dplyr::desc(.data$n_distinct_raw))

cat("\n================ Swapped-order DUPLICATES (key-level, with IDs) ================\n")
if (nrow(dup_tbl) == 0) {
  cat("No swapped-order duplicates detected (order-invariant keys are unique).\n")
} else {
  print(dup_tbl, row.names = FALSE)
}

# -----------------------------------------------------------------------------
# OPTIONAL BUT RECOMMENDED: expanded raw rows for ALL duplicated keys
#     - This prevents the “it became smaller” confusion: this table keeps ALL rows
# -----------------------------------------------------------------------------
dup_rows <- circuitries_original %>%
  dplyr::mutate(
    Circuitries_id_norm = normalize_circuitry_id_strict(.data$Circuitries_id),
    circuitry_key       = make_order_invariant_key(.data$Circuitries_id)
  ) %>%
  dplyr::semi_join(dup_tbl %>% dplyr::select(.data$circuitry_key), by = "circuitry_key") %>%
  dplyr::arrange(.data$circuitry_key, .data$Circuitries_id_norm)

cat("\n--- Raw rows belonging to duplicated keys (dup_rows) ---\n")
cat("nrow(dup_rows) = ", nrow(dup_rows), "\n", sep = "")

# ------------------------------------------------------------------------------
# Snippet to generate and annotate swapped regulatory circuitries
# ------------------------------------------------------------------------------

library(dplyr)

circuitries_annotated <- circuitries %>%
  mutate(
    # (1) Observed circuitry identifier, strictly normalized
    #     This preserves the original ordering as reported in the dataset
    circuitry_id_normalized = normalize_circuitry_id_strict(.data$Circuitries_id),
    
    # (2) Canonical biological identity (order-invariant)
    #     This collapses equivalent circuitries regardless of signature–interaction order
    circuitry_id_canonical  = make_order_invariant_key(.data$Circuitries_id),
    
    # (3) Observed orientation relative to the canonical representation
    #     Indicates whether the circuitry appears in canonical or swapped order
    circuitry_orientation = if_else(
      circuitry_id_normalized == circuitry_id_canonical,
      "canonical_order",
      "swapped_order"
    )
  ) %>%
  group_by(.data$circuitry_id_canonical) %>%
  mutate(
    # (4) Final flag reported in the manuscript
    #     Marks circuitries that appear more than once with opposite orientations
    swapped_circuitry = if_else(n() > 1, "Yes", "No")
  ) %>%
  ungroup()

# ------------------------------------------------------------------------------
# Export annotated circuitry tables
# ------------------------------------------------------------------------------

export(circuitries_annotated, "Dataset_S5_02.tsv")
export(circuitries_annotated, "Dataset_S5_amended.tsv")
export(circuitries_annotated, "Dataset_S5_amended.rds")
