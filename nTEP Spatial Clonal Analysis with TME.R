#############################################
# nTEP BARCODE + TME MARKER SPATIAL ANALYSIS
#
# PURPOSE:
#   - Load nuclei measurements (nTEP barcodes) from CellProfiler.
#   - Normalize epitope intensities by mCherry.
#   - Apply QC filters and optimize threshold for doublets (2-barcode cells).
#   - Cluster nuclei by normalized epitope profile and assign nTEP codes.
#   - Export QC tables and cluster/nTEP calls.
#   - Load TME marker datasets (PDPN, aSMA, CD31, IBA1) from CellProfiler.
#   - Combine tumour clones + TME markers into a single object.
#
# INPUT:
#   - Nuclei CSV exported from CellProfiler
#   - Cellular marker CSVs (PDPN, aSMA, CD31, IBA1) exported from CellProfiler
#
# OUTPUT:
#   - CSV: nTEP codes + QC-passed nuclei
#   - CSV: QC summary table for all nuclei
#   - Spatial plots for tumour clones and TME markers
#############################################


#############################################
# 1. SETUP
#############################################

# ---- Working directory (change to your own) ----
setwd("/Users/silva01/Desktop/")

# ---- PACKAGE INSTALLATION (RUN ONCE, COMMENT OUT AFTER) ----
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("Rtsne")
# install.packages("RColorBrewer")
# install.packages("viridisLite")
# install.packages("ComplexHeatmap")
# install.packages("psych")
# install.packages("ggplot2")
# install.packages("cowplot")
# install.packages("preprocessCore")
# install.packages("tidyverse")
# install.packages("fastcluster")
#
# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("preprocessCore")

# ---- Load libraries ----
library(data.table)
library(dplyr)
library(Rtsne)
library(RColorBrewer)
library(viridisLite)
library(ComplexHeatmap)
library(psych)
library(ggplot2)
library(cowplot)
library(preprocessCore)
library(tidyverse)
library(fastcluster)

# ---- Load custom helper functions (if needed) ----
source("/Users/silva01/Desktop/Scripts and Dataset for Spatial Analysis pt.1/core.R")


#############################################
# 2. USER INPUTS / PARAMETERS
#############################################

# ---- Core measurement settings ----
nuclei_channel   <- "mCherry"                  # nuclei stain channel
measure          <- "Intensity_MeanIntensity"  # base measurement prefix

# ---- Epitope channels (names as in CellProfiler) ----
epitopes <- c("VSVG", "STAG", "V5", "HA", "FLAG", "MYC")

channels     <- epitopes
data_columns <- paste0(measure, "_", channels)

# ---- QC + thresholding parameters ----
intensity_thresh <- 0.145      # intensity threshold for nuclei + epitope
pts_alpha        <- 0.05       # alpha for scatter plot points
min_clust_size   <- 200        # minimum cluster size (not heavily used here)
n_clust_cut      <- 16         # number of clusters in hierarchical clustering

# ---- File paths (update these for new datasets) ----

# Nuclei measurements (CellProfiler output)
nuclei_file <- "/Users/silva01/Desktop/Scripts and Dataset for Spatial Analysis pt.1/BL195112Nuclei.csv"

# Output directory / prefix for this dataset
output_prefix <- "10 Tumour Analysis/Comprehensive Dataset S28 BL195111/BL195112_S28"

# TME marker files (CellProfiler exports)

asma_file  <- "Scripts and Dataset for Spatial Analysis pt.2/MyExpt_Test_aSMA_aSMA_Memb.csv"
cd31_file  <- "Scripts and Dataset for Spatial Analysis pt.3/MyExpt_Test_CD31_Memb.csv"
iba1_file  <- "Scripts and Dataset for Spatial Analysis pt.3/MyExpt_Test_IBA1_IBA1_Memb.csv"
pdpn_file  <- "Scripts and Dataset for Spatial Analysis pt.3/MyExpt_Test_PDPN_PDPN_Memb.csv"

#############################################
# 3. LOAD NUCLEI DATA
#############################################

# Read nuclei measurements from CellProfiler
d_list <- list()
d_list[[1]] <- fread(nuclei_file)

# Assign ImageNumber to each nucleus (field index)
posVect <- unlist(lapply(1:length(d_list), function(x) {
  rep(x, nrow(d_list[[x]]))
}))

# Combine all positions into one data frame
d <- bind_rows(d_list)
d$ImageNumber <- posVect

class(d)


#############################################
# 4. BUILD EPITOPE INTENSITY MATRIX & BASIC QC PLOTS
#############################################

# Extract matrix of epitope intensities
mat <- d %>%
  select(all_of(data_columns)) %>%
  data.matrix()

# Plot epitope vs nuclei intensity and correlation per epitope
par(mfrow = c(3, 3))

for (i in seq_along(channels)) {
  
  x <- data.frame(d)[, data_columns[i]]
  y <- d[[paste0("Intensity_MeanIntensity_", nuclei_channel)]]
  
  # Nuclei that pass both epitope and nuclei threshold
  idx <- x > intensity_thresh & y > intensity_thresh
  
  cor_test <- cor.test(x[idx], y[idx])
  print(cor_test)
  
  fit <- lm(y[idx] ~ x[idx])
  
  pts_col <- rep(rgb(0, 0, 0, pts_alpha), length(x))
  pts_col[!idx] <- rgb(55, 126, 184, pts_alpha * 255, maxColorValue = 255)  # blue
  
  plot(
    x, y,
    xlab = paste0(epitopes[i], " (mean nuclei intensity)"),
    ylab = paste0(nuclei_channel, " (mean nuclei intensity)"),
    xlim = c(0, 1), ylim = c(0, 1),
    col  = pts_col,
    cex  = 0.5
  )
  
  abline(fit, col = brewer.pal(9, "Set1")[1])
  legend("topright",
         legend = paste0("r=", format(cor_test$estimate, digits = 3)),
         bty    = "n"
  )
}

# Barplot: number of epitope channels above threshold per nucleus
barplot(
  table(apply(mat > intensity_thresh, 1, sum)),
  ylab = "Nuclei",
  xlab = paste0("nTEP channels > ", intensity_thresh, " intensity")
)


#############################################
# 5. HELPER: CHANNEL HISTOGRAMS
#############################################

plotHistChannels <- function(mat, names, ...) {
  # Plot histograms of per-channel intensities for QC
  dark_gray  <- rgb(120, 120, 120, maxColorValue = 255)
  
  par(mfrow = c(3, 3))
  par(lwd = 0.2)
  
  for (i in 1:ncol(mat)) {
    hist(mat[, i],
         ylab   = "Nuclei frequency",
         xlab   = names[i],
         ylim   = c(0, 1000),
         breaks = 200,
         col    = dark_gray,
         border = dark_gray,
         ...
    )
  }
}

# Inspect raw mean intensities per epitope
plotHistChannels(mat, names = epitopes, main = "Mean nuclei intensity")


#############################################
# 6. NORMALIZATION BY NUCLEI STAIN (mCherry)
#############################################

# Normalize epitope intensities by nuclei (mCherry)
mat_norm <- sweep(mat, 1, d[[paste0("Intensity_MeanIntensity_", nuclei_channel)]], "/")

plotHistChannels(mat_norm, names = epitopes,
                 main = paste0(nuclei_channel, " normalized"))

# Clamp low-intensity values to NA to avoid over-normalizing dim cells
mat_norm[mat < intensity_thresh] <- NA

plotHistChannels(mat_norm, names = epitopes,
                 main = paste0("Clamped ", nuclei_channel, " normalized"))

# Quantile normalize across epitope channels
mat_norm <- normalize.quantiles(mat_norm)
colnames(mat_norm) <- paste0(
  sapply(strsplit(colnames(mat), "_"), function(x) x[3]),
  "_norm"
)

plotHistChannels(mat_norm, names = epitopes, main = "Quantile normalized")

# Replace NA by zero post-normalization
mat_norm[is.na(mat_norm)] <- 0

# Append normalized epitope intensities to nuclei data
d <- cbind(d, mat_norm)


#############################################
# 7. QC FILTERS FOR NUCLEI
#############################################

# Start QC index: bright enough nuclei
idx <- d[[paste0("Intensity_MeanIntensity_", nuclei_channel)]] > intensity_thresh

message("QC excluding dim ", nuclei_channel, ". Nuclei remaining: ", sum(idx), " / ", length(idx))

# Store QC flag for low nuclei intensity
d[[paste0("QC_low_", nuclei_channel)]] <-
  d[[paste0("Intensity_MeanIntensity_", nuclei_channel)]] <= intensity_thresh

# Require at least one epitope channel above threshold
idx <- idx & apply(mat > intensity_thresh, 1, sum) > 0
message("QC excluding nuclei without nTEP signal. Nuclei remaining: ", sum(idx), " / ", length(idx))

# Store QC flag: no epitope channel above threshold
d$QC_low_epitope <- apply(mat > intensity_thresh, 1, sum) <= 0


#############################################
# 8. DISTRIBUTION OF NORMALIZED SIGNALS AFTER QC
#############################################

summary(mat_norm)

par(mfrow = c(2, 3))
dark_gray  <- rgb(120, 120, 120, maxColorValue = 255)

# Histogram of max normalized signal per nucleus (QC-passed only)
par(lwd = 0.2)
hist(apply(mat_norm[idx, ], 1, max, na.rm = TRUE),
     ylab   = "Nuclei frequency",
     xlab   = "Max normalized channel",
     main   = "",
     breaks = 200,
     col    = dark_gray,
     border = dark_gray
)
par(lwd = 1)


#############################################
# 9. OPTIMIZE THRESHOLD FOR DOUBLETS (2-BARCODE CELLS)
#############################################

thresh_range <- seq(0, 1.5, length.out = 100)

# Cell code multiplicity across threshold range
cell_multiplicity <- sapply(thresh_range, function(thresh) {
  table(
    factor(
      apply(mat_norm > thresh, 1, sum, na.rm = TRUE),
      levels = 0:ncol(mat_norm[idx, ])
    )
  )
})

# Threshold that maximizes number of doublets
opt_thresh <- thresh_range[which.max(cell_multiplicity[3, ])]

# Number of channels above optimal threshold per nucleus
n_channels <- apply(mat_norm > opt_thresh, 1, sum, na.rm = TRUE)
d$QC_n_channels <- n_channels

# Additional QC: keep only doublets (2 channels positive)
idx <- idx & n_channels == 2
message("QC excluding non-binary nuclei. Nuclei remaining: ", sum(idx), " / ", length(idx))

# Optimization curves (singles/doubles/none vs threshold)
colors     <- brewer.pal(9, "Set1")[-6]
line_width <- 1.5

plot(
  thresh_range, cell_multiplicity[2, ],
  ylim = range(cell_multiplicity[1:3, ]),
  col  = colors[1],
  type = "l",
  lwd  = line_width,
  xlab = paste0("Channel detection threshold (x / ", nuclei_channel, ")"),
  ylab = "Nuclei",
  main = ""
)

lines(thresh_range, cell_multiplicity[3, ], col = colors[2], lwd = line_width, type = "l")
lines(thresh_range, cell_multiplicity[1, ], col = colors[3], lwd = line_width, type = "l")

abline(v = opt_thresh, lty = 2, col = dark_gray)

legend("topleft",
       legend = c("Singles", "Doubles", "None"),
       col    = colors[1:3],
       pch    = 15
)

# Barplot of multiplicity at optimal threshold
light_gray <- rgb(240, 240, 240, maxColorValue = 255)
bar_colors <- rep(light_gray, length(table(n_channels)))
bar_colors[3] <- brewer.pal(9, "Set1")[2]   # highlight doublets

barplot(table(n_channels),
        space = 0,
        col   = bar_colors,
        ylab  = "Nuclei frequency",
        xlab  = paste0("nTEP code multiplicity at x / ", nuclei_channel,
                       " > ", format(opt_thresh, digits = 3))
)

# Barplot nuclei per image region: QC+ vs QC-
d_sub <- d[idx, ]

image_counts <- rbind(
  table(d_sub$ImageNumber),
  table(d$ImageNumber) - table(d_sub$ImageNumber)
)

barplot(image_counts,
        legend = c("QC +", "QC -"),
        ylab   = "Nuclei",
        xlab   = "Image region (field)"
)

# Extract normalized matrix for QC-passed nuclei
mat_sub <- mat_norm[idx, ]

# Boxplot of normalized signals for QC-passed nuclei
boxplot(mat_sub,
        names = colnames(mat_sub),
        las   = 2,
        ylab  = paste0("x / ", nuclei_channel, " (QC passed)")
)


#############################################
# 10. FINAL NORMALIZATION FOR CLUSTERING
#############################################

# Store final QC inclusion flag
d$QC_include <- idx

# Recompute QC-passed objects for clarity
d_sub  <- d[idx, ]
mat_sub <- mat_norm[idx, ]

# Pseudo log2 transform
mat_sub <- log2(mat_sub + 1)


#############################################
# 11. tSNE VISUALIZATION OF BARCODE SPACE
#############################################

set.seed(42)
tsne <- Rtsne(mat_sub, check_duplicates = FALSE)

d_sub$tSNE_X <- tsne$Y[, 1]
d_sub$tSNE_Y <- tsne$Y[, 2]

# tSNE coloured by epitope intensity
plts <- lapply(epitopes, function(ch) {
  ggplot(d_sub, aes_string(
    x   = "tSNE_X",
    y   = "tSNE_Y",
    col = paste0(ch, "_norm")
  )) +
    scale_y_reverse() +
    scale_colour_gradientn(colors = rev(magma(20))) +
    geom_point(size = 0.5) +
    theme_classic()
})

plot_grid(plotlist = plts)

# tSNE coloured by ImageNumber (batch effect check)
ggplot(d_sub, aes(
  x   = tSNE_X,
  y   = tSNE_Y,
  col = ImageNumber
)) +
  scale_y_reverse() +
  geom_point() +
  theme_classic()


#############################################
# 12. HIERARCHICAL CLUSTERING & CLUSTER CENTERS
#############################################

min_sig_noise_ratio <- 2   # (defined but not used here)

# Run hierarchical clustering on QC-passed nuclei
clust <- cutree(
  fastcluster::hclust(
    dist(mat_sub, method = "manhattan"),
    method = "average"
  ),
  k = n_clust_cut
)

table(clust)

# Cluster centers: mean log-normalized signal per epitope
centers <- aggregate(mat_sub, list(clust), mean)[, -1]
centers <- data.matrix(centers)

rownames(centers) <- 1:nrow(centers)
colnames(centers) <- sapply(strsplit(colnames(centers), "_"), function(x) x[1])

# Save cluster ID in d_sub
d_sub$clust <- factor(clust)

# Signal-to-noise ratio per cluster (top 2 vs rest)
sig_noise_ratio <- apply(centers, 1, function(values) {
  sorted_values <- sort(values, decreasing = TRUE)
  sum(sorted_values[1:2]) / sum(sorted_values[-1:-2])
})


#############################################
# 13. ASSIGN nTEP CODES TO CLUSTERS
#############################################

# Rank epitope names per cluster by mean signal
code_rank <- apply(centers, 1, function(values) {
  names(sort(values, decreasing = TRUE))
})
code_rank <- t(code_rank)

code_basis <- 2   # encode cluster by top 2 epitopes

ntep_codes <- data.frame(
  clust     = factor(1:nrow(centers)),
  ntep_code = apply(code_rank[, 1:code_basis], 1, function(elem) {
    paste0(sort(elem), collapse = "_")   # sort so A_B == B_A
  })
)

# Add cluster size suffix (e.g. CODE_123)
num_clusts <- as.integer(table(clust))
ntep_codes <- ntep_codes %>%
  mutate(ntep_code_n = paste0(ntep_code, "_", num_clusts))

# Merge back into nuclei data
d_sub <- merge(d_sub, ntep_codes, by = "clust")


#############################################
# 14. FILTER CLUSTERS / CODES & HEATMAP OF CENTERS
#############################################

# Keep clusters with at least a minimum number of nuclei
include_clusts <- table(d_sub$clust) > 150
include_clusts
sum(include_clusts)

# Mark excluded clusters as NA for downstream plots
d_sub$clust_min <- d_sub$clust
d_sub$clust_min[d_sub$clust %in% which(!include_clusts)] <- NA
d_sub$clust_min <- factor(d_sub$clust_min)

d_sub$ntep_code_n[d_sub$clust %in% which(!include_clusts)] <- NA

# Determine row order for heatmap using base-10 expansion heuristic
row_order <- order(
  round(centers[include_clusts, ]) %*% 10^(ncol(centers):1),
  decreasing = TRUE
)

# Annotate heatmap with cluster sizes
ha <- rowAnnotation(
  Nuclei = anno_barplot(
    as.vector(table(d_sub$clust)[include_clusts][row_order])
  )
)

# Heatmap of cluster mean intensities
Heatmap(
  centers[include_clusts, ][row_order, ],
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  name            = "Relative mean intensity (a.u.)",
  col             = viridis(20),
  right_annotation = ha
)


#############################################
# 15. tSNE & SPATIAL PLOTS COLOURED BY nTEP CODE
#############################################

# tSNE coloured by nTEP code
ggplot(d_sub, aes(
  x   = tSNE_X,
  y   = tSNE_Y,
  col = ntep_code_n
)) +
  scale_y_reverse() +
  geom_point() +
  theme_classic()

# Spatial distribution of nTEP codes by image field
plts <- lapply(unique(d$ImageNumber), function(k) {
  d_field <- d_sub %>%
    filter(ImageNumber == k)
  
  ggplot(d_field, aes(
    x   = Location_Center_X,
    y   = Location_Center_Y,
    col = ntep_code_n,
    alpha = 1
  )) +
    scale_y_reverse() +
    scale_color_discrete(drop = FALSE) +  # consistent colors across fields
    geom_point(size = 0.01) +
    theme_minimal()
})

plot_grid(plotlist = plts)


#############################################
# 16. EXPORT nTEP CODES & QC TABLES
#############################################

# Export codes + basic morphology for QC-passed nuclei
export_columns_codes <- c(
  "ImageNumber", "ObjectNumber",
  "clust_min", "ntep_code",
  "Location_Center_X", "Location_Center_Y",
  "AreaShape_Area", "AreaShape_Perimeter",
  paste0(epitopes, "_norm")
)

write.csv(
  data.frame(d_sub)[, export_columns_codes],
  paste0(output_prefix, "_nTEP_codes.csv"),
  row.names = FALSE
)

# Export QC summary table for all nuclei
export_columns_qc <- c(
  "ImageNumber", "ObjectNumber",
  "Location_Center_X", "Location_Center_Y",
  "QC_include", "AreaShape_Area", "AreaShape_Perimeter",
  paste0("QC_low_", nuclei_channel),
  "QC_low_epitope", "QC_n_channels",
  "Intensity_MeanIntensity_mCherry",
  paste0(epitopes, "_norm")
)

write.csv(
  data.frame(d)[, export_columns_qc],
  paste0(output_prefix, "_nTEP_QC.csv"),
  row.names = FALSE
)


#############################################
# 17. LOAD TME MARKERS & BASIC PLOTS
#############################################

# Load PDPN+ fibroblast markers
R3 <- read.csv(file = pdpn_file)

ggplot(R3, aes(
  x = Location_Center_X,
  y = Location_Center_Y
)) +
  scale_y_reverse() +
  geom_point(size = 0.05, alpha = 0.5, shape = 3) +
  theme_minimal()

# Load aSMA fibroblast markers
R4 <- read.csv(file = asma_file)

ggplot(R4, aes(
  x = Location_Center_X,
  y = Location_Center_Y
)) +
  scale_y_reverse() +
  geom_point(size = 0.005, alpha = 0.5, shape = 6, colour = "red") +
  theme_minimal()

# Tag each TME dataset with a "stain" label
R3.1 <- mutate(R3, stain = "PDPN")
R4.1 <- mutate(R4, stain = "aSMA")

R5 <- bind_rows(R3.1, R4.1)
R5$stain <- as.factor(R5$stain)

ggplot(R5, aes(
  x = Location_Center_X,
  y = Location_Center_Y,
  group = stain,
  shape = stain,
  colour = stain
)) +
  scale_y_reverse() +
  geom_point(size = 0.05, alpha = 1) +
  theme_minimal()

# Load CD31 (endothelial) markers
R6 <- read.csv(file = cd31_file)

ggplot(R6, aes(
  x = Location_Center_X,
  y = Location_Center_Y
)) +
  scale_y_reverse() +
  geom_point(size = 0.05, alpha = 0.5, shape = 7) +
  theme_minimal()

R6.1 <- mutate(R6, stain = "CD31")

R7 <- bind_rows(R5, R6.1)

ggplot(R7, aes(
  x = Location_Center_X,
  y = Location_Center_Y,
  group = stain,
  shape = stain,
  colour = stain
)) +
  scale_colour_manual(values = c("#1AE9F8", "#3FE533", "#AE56CE")) +
  scale_y_reverse() +
  geom_point(size = 0.05, alpha = 1) +
  theme_minimal()

# Load IBA1 (macrophage) markers
R8 <- read.csv(file = iba1_file)

ggplot(R8, aes(
  x = Location_Center_X,
  y = Location_Center_Y
)) +
  scale_y_reverse() +
  geom_point(size = 0.01, alpha = 0.5, shape = 7) +
  theme_minimal()

R8.1 <- mutate(R8, stain = "IBA1")

R9 <- bind_rows(R7, R8.1)

ggplot(R9, aes(
  x = Location_Center_X,
  y = Location_Center_Y,
  group = stain,
  shape = stain,
  colour = stain
)) +
  scale_colour_manual(values = c("#1AE9F8", "#3FE533", "#FBDD21", "#AE56CE")) +
  scale_shape_manual(values = c(0, 2, 5, 4)) +
  scale_y_reverse() +
  geom_point(size = 0.01, alpha = 1) +
  theme_minimal()


#############################################
# 18. MERGE TME MARKERS WITH BARCODED CANCER CELLS
#############################################

# Remove a specific false-positive nTEP code if needed
d_sub_1 <- d_sub[!d_sub$ntep_code_n == "FLAG_VSVG_5248", ]
table(d_sub$ntep_code_n)
table(d_sub_1$ntep_code_n)

# Quick spatial plot of all tumour clones
ggplot(d_sub_1, aes(
  x   = Location_Center_X,
  y   = Location_Center_Y,
  col = ntep_code_n
)) +
  scale_y_reverse() +
  scale_color_discrete(drop = FALSE) +
  geom_point(size = 0.01, alpha = 1) +
  theme_minimal()

# Tag cancer cells with "Cancer"
d_sub_1 <- mutate(d_sub_1, stain = "Cancer")

# Combine cancer cells and TME markers
R10 <- bind_rows(R9, d_sub_1)
R10$stain <- as.factor(R10$stain)

# Combine stain + nTEP code into a single label
R10 <- R10 %>%
  mutate(stain = paste0(stain, "_", ntep_code_n))

table(R10$stain)

# Combined spatial overview of clones + TME markers
ggplot(R10, aes(
  x     = Location_Center_X,
  y     = Location_Center_Y,
  shape = stain,
  col   = stain
)) +
  scale_shape_manual(values = c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 5, 4)) +
  scale_colour_manual(values = c("#18F8DF", "#EB8262", "#D47525", "#5C8611",
                                 "#2E9687", "#2777A4", "#33A2E0", "#6161D2",
                                 "#C564D6", "#DB2CC7", "#DE3B26",
                                 "#3FE533", "#FBDD21", "#9D33DD")) +
  scale_y_reverse() +
  geom_point(size = 0.01, alpha = 0.5) +
  theme_minimal()

