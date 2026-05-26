# title: Create MEDICC2 input files by combining cluster-level TSVs into region-level TSVs
# - takes cluster-level TSV files as input (from 10X DLP data)
# - outputs region-level TSV files for MEDICC2 input

# Load libraries
library(data.table)

# Set folder structure
minseg     <- "2.5"
input_dir  <- paste0("/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/Haixi/scDNA/MEDICC2_input/per_cluster_10X_DLP/minseg_", minseg)
output_dir <- paste0("/mnt/iribhm/homes/aste0033/projects/CNA-PhyloAnalysis/Haixi/scDNA/MEDICC2_input/per_region_10X_DLP/minseg_", minseg)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define regions and read cluster-level TSV files, combine them into region-level TSV files
regions <- c("P", "R1", "R2", "R3", "R4", "R5")
tsv_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Loop through regions, find corresponding cluster-level TSV files, combine them, and write region-level TSV files
for (region in regions) {
  region_files <- tsv_files[grepl(paste0("_", region, "_"), basename(tsv_files))]

  if (length(region_files) == 0) {
    message("No files found for region: ", region)
    next
  }

  combined <- rbindlist(lapply(region_files, fread, sep = "\t"))

  out_file <- file.path(output_dir, paste0(region, ".tsv"))
  fwrite(combined, out_file, sep = "\t")
  message("Written: ", out_file, " (", nrow(combined), " rows from ", length(region_files), " clusters)")
}
