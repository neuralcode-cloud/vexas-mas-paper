# ==============================================================================
# PCA of Olink Reveal proteomic data -- Figure 2E
#
# Manuscript: "VEXAS and Macrophage Activation Syndrome: clinical and multiomic
#              approach of a unique hyperinflammatory clinical dyad"
# Original script as received from the study team: original/PCA_Olink_data.R
#
# This version starts from the de-identified, already-QC-filtered NPX table that
# ships with the repository (data/proteomics/olink_long.csv), so it runs without
# the raw plate export and without the sample-identity metadata. See
# README.md for the five changes relative to the original.
#
# Usage:  Rscript analysis/proteomics/run_pca.R
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", ".."),
         mustWork = FALSE)
if (is.na(ROOT) || !dir.exists(ROOT)) ROOT <- normalizePath(".")

DATA <- file.path(ROOT, "data", "proteomics")
OUT  <- file.path(ROOT, "outputs", "proteomics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LOD_THRESHOLD <- 95   # discard an assay if > this % of samples fall below LOD

# ------------------------------------------------------------------ 1. import
npx  <- read.csv(file.path(DATA, "olink_long.csv"), check.names = FALSE)
meta <- read.csv(file.path(DATA, "sample_metadata.csv"), check.names = FALSE)

# The exported CSV carries an unnamed row-index column; drop it (dplyr refuses
# to operate on a frame with blank column names).
npx <- npx[, names(npx) != "" & !is.na(names(npx)), drop = FALSE]

stopifnot(all(c("SampleID", "Assay", "NPX", "LOD") %in% names(npx)))
stopifnot(setequal(unique(npx$SampleID), meta$sample_id))

cat("samples:", length(unique(npx$SampleID)),
    " assays:", length(unique(npx$Assay)), "\n")

# ------------------------------------------------- 2. QC (already applied)
# The shipped table contains only AssayQC == "PASS" / SampleQC == "PASS" rows
# and no Olink internal controls, so the original script's QC step is a no-op
# here. Asserted rather than assumed.
if ("AssayQC" %in% names(npx))
  stopifnot(all(npx$AssayQC == "PASS"))
if ("SampleType" %in% names(npx))
  stopifnot(all(npx$SampleType == "SAMPLE"))

# --------------------------------------------------------- 3. LOD filtering
lod_status <- function(data, threshold) {
  total_samples <- length(unique(data$SampleID))
  below <- data %>% filter(NPX < LOD)

  per_assay <- below %>%
    count(Assay, name = "n_below") %>%
    mutate(pct_below = round(n_below / total_samples * 100, 2))

  discard <- per_assay$Assay[per_assay$pct_below > threshold]

  cat("\n  values below LOD\n")
  cat("    samples with >=1 below LOD:", length(unique(below$SampleID)),
      "/", total_samples, "\n")
  cat("    assays  with >=1 below LOD:", length(unique(below$Assay)),
      "/", length(unique(data$Assay)), "\n")
  cat("    assays >", threshold, "% below LOD:", length(discard), "\n")

  list(per_assay = per_assay, discard = discard, total_samples = total_samples)
}

lod <- lod_status(npx, LOD_THRESHOLD)

# With 24 samples one sample is 4.17%, so a ">95%" rule is satisfied by 23/24
# (95.83%) as well as by 24/24. The count reported in the manuscript (15
# discarded, 1018 retained) corresponds to the 24/24 rule. Both are printed.
n_all_below <- sum(lod$per_assay$n_below == lod$total_samples)
cat("    of which 100% below LOD  :", n_all_below, "\n")
cat("    -> retained at >", LOD_THRESHOLD, "%: ",
    length(unique(npx$Assay)) - length(lod$discard), "\n", sep = "")
cat("    -> retained at 100% rule : ",
    length(unique(npx$Assay)) - n_all_below, "\n", sep = "")

npx_final <- npx %>% filter(!Assay %in% lod$discard)

write.csv(lod$per_assay, file.path(OUT, "lod_per_assay.csv"), row.names = FALSE)
writeLines(as.character(lod$discard), file.path(OUT, "assays_discarded.txt"))

# ------------------------------------------------------------------- 4. PCA
mat <- npx_final %>%
  select(SampleID, Assay, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# align metadata to the matrix row order (never assume the orders agree)
meta_ord <- meta[match(rownames(mat), meta$sample_id), ]
stopifnot(identical(meta_ord$sample_id, rownames(mat)))
stopifnot(!any(is.na(mat)))

pca <- prcomp(mat, scale. = TRUE)
ve <- round(summary(pca)$importance[2, 1:2] * 100, 1)

scores <- as.data.frame(pca$x[, 1:4]) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta_ord, by = "sample_id")

write.csv(scores, file.path(OUT, "pca_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = seq_len(min(10, ncol(pca$x))),
                     pct_variance = round(summary(pca)$importance[2, 1:min(10, ncol(pca$x))] * 100, 3)),
          file.path(OUT, "pca_variance.csv"), row.names = FALSE)

# ---------------------------------------------------------------- 5. figure
p <- ggplot(scores, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3.4, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = subject), size = 2.6,
                           show.legend = FALSE, max.overlaps = 20) +
  labs(x = paste0("PC1 (", ve[1], "%)"),
       y = paste0("PC2 (", ve[2], "%)"),
       colour = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.72)

ggsave(file.path(OUT, "figure_2e_pca.png"), p, width = 6.6, height = 5.2, dpi = 300)

cat("\nPC1:", ve[1], "%  PC2:", ve[2], "%\n")
cat("wrote:", file.path(OUT, "figure_2e_pca.png"), "\n")

# ------------------------------------------------------- 6. reproducibility
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
