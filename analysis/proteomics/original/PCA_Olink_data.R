# ==============================================================================
# Script Name: PCA for Olink Data
# Manuscript: "VEXAS and Macrophage Activation Syndrome: clinical and multiomic approach of a unique hyperinflammatory clinical dyad"
# Provenance: script as received from the study team. The only change is the
# removal of personal credit lines and internal laboratory identifiers; no
# code or logic was altered. See ../README.md.

# 1. LIBRARIES and FUNCTIONS ---------------------------------------------------------
if (!requireNamespace("OlinkAnalyze", quietly = TRUE))
  install.packages("OlinkAnalyze")
library(OlinkAnalyze)

library(readxl)
library(ggplot2)
library(tibble)
library(tidyverse)
library(ggfortify)

#### Function for getting info regarding LOD status ####
analyze_lod_status <- function(data, LOD_thesh) {
  # Filter only values below LOD
  df_below <- data %>% filter(NPX < LOD)
  
  # Get statistical info
  n_samples_lod <- length(unique(df_below$SampleID))
  n_assays_lod  <- length(unique(df_below$Assay))
  total_samples <- length(unique(data$SampleID))
  total_assays  <- length(unique(data$Assay))
  
  cat("\t Distribution of values < LOD \n")
  cat("Samples with at least 1 value < LOD: ", n_samples_lod, "/", total_samples, "\n")
  cat("Assays with at least 1 value < LOD:  ", n_assays_lod, "/", total_assays, "\n")
  
  # Get percentage of values < LOD per Assay
  df_assay_stats <- as.data.frame(table(df_below$Assay))
  colnames(df_assay_stats) <- c("Assay", "Freq_BelowLOD")
  df_assay_stats$Perc_BelowLOD <- round((df_assay_stats$Freq_BelowLOD / total_samples) * 100, 2)
  
  # Get list of discardable proteins 
  discardable <- df_assay_stats$Assay[df_assay_stats$Perc_BelowLOD > LOD_thesh]
  cat("Proteins with >", LOD_thesh,"% values below LOD:", length(discardable), "\n")
  
  # Statistics on percentage of values < LOD per SampleID
  df_sample_stats <- as.data.frame(table(df_below$SampleID))
  colnames(df_sample_stats) <- c("SampleID", "Freq_BelowLOD")
  original_order <- unique(data$SampleID) 
  df_sample_stats$SampleID <- factor(df_sample_stats$SampleID, levels = original_order)
  
  df_sample_stats$Perc_BelowLOD <- round((df_sample_stats$Freq_BelowLOD / total_assays) * 100, 2)
  
  return(list(
    summary_assay = df_assay_stats,
    summary_sample = df_sample_stats,
    discard_list = discardable,
  ))
}

#### Function for PCA ####
PCA_plot <- function(db_pca_npx, db_meta, type, var_for_color){
  pca <- prcomp(db_pca_npx, scale. = T)
  
  pca_plot_ann <- autoplot(pca, data = db_meta, label = F, colour = var_for_color,
                           label.colour = "black", size = 5, shape = 16, alpha = 0.8) +
    geom_text(aes(label = PatientID), size = 2) + theme_minimal() + coord_cartesian(clip = "off") + 
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16), title = element_text(size = 16),
          legend.text = element_text(size = 8), legend.title = element_blank(), 
          legend.box.just = "left", legend.justification = c(0.9,1), 
          legend.position = "top"
    )
  
  print(pca_plot_ann)
  ggsave(paste0("./results/pca_", type, "_", var_for_color, ".png"), pca_plot_ann, width = 6, height = 4, dpi = 1200)
  
  return(pca = pca)
}

# 2. DATA IMPORT ---------------------------------------------------------------

# Upload npx matrix of Reveal plate
import <- OlinkAnalyze::read_NPX("path/to/olinkplate/export")
# Upload annotation file for timepoint
metadata <- readxl::read_xlsx("./data/ListOfSamples 18-03-2026.xlsx", sheet = 2, range = "A1:E87")
# Rename columns
colnames(metadata) <- c("PatientID", "ExpSampleID", "Age", "Sex", "Annotation")
# Keep the Annotation column order unchanged
metadata$Annotation <- factor(metadata$Annotation, levels = unique(metadata$Annotation))
# Generate variable for timepoint
metadata <- metadata %>%
  group_by(PatientID) %>%
  mutate(Timepoint = ifelse(
    grepl("post", Annotation),
    paste0("T", cumsum(grepl("post", Annotation))),
    "T0")
  ) %>% 
  ungroup()
# Simplify Annotation name of UBA1- group
metadata$Annotation <- gsub("UBA1-negative with no hematological disease", "UBA1-", metadata$Annotation)
# Generate unique variable for sample
metadata$SampleID <- paste0(metadata$PatientID, "_", metadata$Timepoint)

# 3. QUALITY CONTROL and DATA FILTERING ----------------------------------------

# Select internal Olink controls 
ctrl_samples <- import %>% 
  dplyr::filter(SampleType %in% c("SAMPLE_CONTROL", "PLATE_CONTROL", "NEGATIVE_CONTROL")) %>% 
  pull(SampleID) %>% 
  unique()
# Collect the Assays for which the Quality Control produce a WARN
assayQC <- which(import$AssayQC == "WARN")
# Remove proteins for which Warning is observed
import <- import %>%
  filter(Assay != unique(import$Assay[assayQC]))

# LOD calculation and filtering
# Upload the file from Olink website, collecting fixedLODs for Olink Reveal
fixedLOD_path <- "data/Reveal Fixed LOD - csv file.csv"
# Calculation of LOD using the FixedLOD
npx <- olink_lod(import, lod_file_path = fixedLOD_path, lod_method = "FixedLOD")
# Get list of Assays for which NA value is observed in LOD
na_assays <- unique(npx$Assay[is.na(npx$LOD)])

# Remove plate sample controls and plate assay controls from the dataset
npx_reveal_filt <-  npx %>%
  filter(!SampleID %in% ctrl_samples) %>%
  filter(!Assay %in% na_assays) %>%
  arrange(match(SampleID, metadata$ExpSampleID))

# Convert npx matrix SampleID to match with metadata
npx_reveal_filt$SampleID <- metadata$SampleID[match(npx_reveal_filt$SampleID, metadata$ExpSampleID)]
npx_reveal_filt <- npx_reveal_filt %>%
  arrange(match(SampleID, metadata$SampleID))

# Get info regarding LOD status
lod_results <- analyze_lod_status(npx_reveal_filt, 95) # there would be 15 discardable Assays

# Select proteins having more than 95% of LOD values
discardable_prot <- lod_results$discard_list

# Filter out proteins having more than 95% of LOD values
npx_final <- npx_reveal_filt %>% 
  filter(!Assay %in% discardable_prot)


# 4. PRINCIPAL COMPONENT ANALYSIS (PCA) ----------------------------------------

# Select metadata of samples of interest
meta_pca_long <- metadata %>%
  group_by(PatientID) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(match(SampleID, unique(npx_final$SampleID))) %>%
  mutate(Annotation = as.factor(Annotation))

rownames(meta_pca_long) <- meta_pca_long$SampleID

# Select raw data of samples of interest
pca_npx_long <- npx_final %>%
  dplyr::filter(SampleID %in% meta_pca_long$SampleID) %>%
  dplyr::select(SampleID, Assay, NPX) %>%
  pivot_wider(names_from = Assay, values_from = NPX) %>%
  column_to_rownames("SampleID") %>%
  as.matrix()

# Get PCA
pca_long <- PCA_plot(pca_npx_long, meta_pca_long, "long", "Annotation") 
pca_long <- PCA_plot(pca_npx_long, meta_pca_long, "long", "Timepoint")


# 5. REPRODUCIBILITY INFO ------------------------------------------------------
sessionInfo()