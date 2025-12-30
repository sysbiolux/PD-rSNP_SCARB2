suppressPackageStartupMessages({
  library(plotgardener)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
  library(org.Hs.eg.db)
  library(png)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggmanh)
  library(RColorBrewer)
})

# 2025-12-10 Aurélien Ginolhac based on Déborah Gérard initial work
# Rendered using a container: `docker run -u $(id -u):$(id -g) -ti -v ${HOME}/Work/250603/Gerard_et_al_2026:/mnt 250603:24.04`
# Be sure that renv is not invoked if you use it: `mv .Rprofile Rprofile`
# Once in the container: cd /mnt/FIGURE2;  Rscript --vanilla Figure2.R

########## FIGURE 1 ##########

pdf(
  "FIGURE2.pdf",
  width = 8.27,
  height = 11.67
)

# Create a A4 blank page
pageCreate(
  width = 8.27,
  height = 11.67,
  default.units = "inches",
  showGuides = FALSE
)

plotText(
  label = "Figure 2",
  fontsize = 14,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 0.25,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL A -                                                              ####
# text A
plotText(
  label = "A",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 0.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL B -                                                           #######
# text B
plotText(
  label = "B",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 3.125,
  y = 0.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

#### PANEL C -                ##################################################
# text C
plotText(
  label = "C",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 5.5,
  y = 0.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL D             ######################################################

plotText(
  label = "D",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 2.0,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

#### PANEL E -                                       ###########################
# text E
plotText(
  label = "E",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 4.0,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

# Load the picture prepared in Biorender
Fig2E_pip <- readPNG("Figure2_PanelE.png")

# Plot Figure 2A - Pipeline
plotRaster(
  image = Fig2E_pip,
  x = 0.5,
  y = 4.1,
  width = 3,
  height = 3,
  just = c("left", "top"),
  interpolate = TRUE
)


#### PANEL F ###################################################################
plotText(
  label = "F",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 3.75,
  y = 4,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


# Define parameters for the regions
region <- pgParams(
  chrom = "chr4",
  chromstart = 89600354,
  chromend = 90020000,
  assembly = "hg38"
)

interact_bed <- read_tsv("../SNPs/arches_in_snps.tsv", show_col_types = FALSE)

snps_254 <- read_tsv(
  "../SNPs/254rSNPs_GWAS_pval_rsID.bed",
  col_names = c("chr", "start", "end", "name", "rsID"),
  show_col_types = FALSE
)

interact_bed |>
  filter(sourceChrom == region$chrom) |>
  select(
    chrom1 = sourceChrom,
    start1 = sourceStart,
    end1 = sourceEnd,
    chrom2 = targetChrom,
    start2 = targetStart,
    end2 = targetEnd,
    value,
    sample
  ) -> STARE_intervals
STARE_intervals |>
  filter(sample == "STARE_ABCpp_mDAN_day30") |>
  ## Translate lengths into heights
  mutate(
    lengths = abs(start2 - start1),
    heights = (lengths / max(lengths)) + 100
  ) -> interval_mDAN
STARE_intervals |>
  filter(sample == "STARE_ABCpp_smNPC") |>
  ## Translate lengths into heights
  mutate(
    lengths = abs(start2 - start1),
    heights = (lengths / max(lengths)) + 100
  ) -> interval_smNPC


bw_path <- "../BIGWIG/"

bw_samples <- c(
  "smNPC_I.bw",
  "D30_POS_I.bw"
)

bw_region <- lapply(bw_samples, \(x) {
  readBigwig(
    file = file.path(bw_path, x),
    chrom = region$chrom,
    chromstart = region$chromstart,
    chromend = region$chromend
  )
})
scale_y_max <- 100

results <- vector(mode = "list", length = length(bw_region))

BW_COLORS <- c(
  "#ff85ff",
  "#8efa00"
)

Y_BW <- c(5.05, 6.1)

SAMPLES <- c(
  "smNPC",
  "mDAN D30"
)


for (i in seq_along(bw_region)) {
  results[[i]] <- plotSignal(
    data = bw_region[[i]],
    params = region,
    fill = BW_COLORS[i],
    linecolor = NA,
    x = 4.15,
    y = Y_BW[i],
    width = 3.75,
    height = 0.25,
    range = c(0, scale_y_max),
    just = c("left", "top"),
    default.units = "inches"
  )

  plotText(
    label = SAMPLES[i],
    fontsize = 6,
    fontcolor = BW_COLORS[i],
    fontfamily = "Helvetica",
    x = 3.8,
    y = Y_BW[i] + 0.1,
    just = c("left", "top"),
    default.units = "inches",
    fontface = "bold"
  )
  annoYaxis(
    plot = results[[i]],
    params = region,
    at = c(0, scale_y_max),
    fontsize = 4
  )
}

plotGenes(
  params = region,
  chrom = region$chrom,
  assembly = assembly(
    Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
    OrgDb = "org.Hs.eg.db"
  ),
  fill = "black",
  fontcolor = "black",
  x = 4.15,
  y = 6.5,
  width = 3.75,
  height = 0.5,
  just = c("left", "top"),
  default.units = "inches"
)

archPlot <- plotPairsArches(
  data = interval_smNPC,
  params = region,
  fill = colorby(
    "value",
    range = c(0.02, 1),
    palette = colorRampPalette(c("grey90", "#ff85ff"))
  ),
  linecolor = "fill",
  archHeight = "heights",
  alpha = 1,
  x = 4.15,
  y = 4.3,
  height = 0.7,
  width = 3.75,
  just = c("left", "top"),
  default.units = "inches"
)

## Annotate genome label
annoGenomeLabel(plot = archPlot, x = 4.15, y = 4.1, scale = "Mb")

## Annotate heatmap legend
annoHeatmapLegend(
  plot = archPlot,
  fontcolor = "black",
  x = 3.9,
  y = 4.4,
  width = 0.08,
  height = 0.5,
  fontsize = 5
)


## arches mDAN
archPlot2 <- plotPairsArches(
  data = interval_mDAN,
  params = region,
  fill = colorby(
    "value",
    # force original range
    range = c(0.02, 1),
    palette = colorRampPalette(c("grey90", "#8efa00"))
  ),
  linecolor = "fill",
  archHeight = "heights",
  alpha = 1,
  x = 4.15,
  y = 5.35,
  height = 0.7,
  width = 3.75,
  just = c("left", "top"),
  default.units = "inches"
)

## Annotate heatmap legend
annoHeatmapLegend(
  plot = archPlot2,
  fontcolor = "black",
  x = 3.9,
  y = 5.4,
  width = 0.08,
  height = 0.5,
  fontsize = 5
)


my_snps <- filter(
  snps_254,
  chr == region$chrom,
  start >= region$chromstart,
  end <= region$chromend
)
for (i in seq_len(nrow(my_snps))) {
  # Create the region using pgParams() for the current SNP
  SNP <- pgParams(
    chrom = my_snps[i, "chr"],
    chromstart = my_snps[i, "start"] - 200,
    chromend = my_snps[i, "end"] + 200,
    label = my_snps[i, "rsID"]
  )

  annoHighlight(
    plot = archPlot2,
    params = SNP,
    fill = "orange2",
    y = 6.4,
    height = 0.1,
    alpha = 1,
    just = c("left", "top"),
    default.units = "inches"
  )
}

plotText(
  params = SNP,
  label = "PD SNPs",
  x = 3.8,
  y = 6.45,
  fontcolor = "orange2",
  just = c("left", "center"),
  fontsize = 6
)

#### PANEL G ###################################################################

plotText(
  label = "G",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 7.75,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

#### PANEL H ###################################################################

plotText(
  label = "H",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 2,
  y = 7.75,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)
PDrSNPs_54_het.ASE.chisq <- read_tsv(
  "PDrSNPs_54_het.ASE.chisq.tsv",
  show_col_types = FALSE
)

PDrSNPs_54_het.ASE.chisq |>
  pivot_longer(cols = -Sample, names_to = "snp", values_to = "p.value") |>
  slice_min(n = 1, order_by = p.value, by = snp) |>
  arrange(p.value) |>
  mutate(cum_dist = cume_dist(p.value), p.value = p.adjust(p.value)) |>
  ggplot(aes(x = p.value, y = cum_dist, colour = p.value < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  geom_point() +
  annotate(geom = "text", label = "74%", x = 0.03, y = 0.95, colour = "red") +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "red") +
  guides(colour = guide_legend(position = "inside")) +
  theme_classic(7) +
  theme(legend.position.inside = c(0.5, 0.4)) +
  labs(
    x = "Chi-square adjusted p-values",
    y = "Cumulative distribution",
    colour = "padj < 0.05",
    title = "Allelic imbalance\nof the 54 PD-SNPs"
  ) -> snp_imbalance

plotGG(
  plot = snp_imbalance,
  x = 1.95,
  y = 8,
  width = 2.3,
  height = 2.4,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL I ###################################################################

plotText(
  label = "I",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 4.2,
  y = 7.75,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


# Manhattan as PNG
plotRaster(
  image = readPNG("manhattan_plot.png"),
  x = 4.25,
  y = 8.5,
  width = 3.75,
  height = 2,
  just = c("left", "top"),
  interpolate = TRUE
)


#pageGuideHide()
dev.off()
