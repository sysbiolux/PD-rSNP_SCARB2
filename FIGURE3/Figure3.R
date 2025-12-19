suppressPackageStartupMessages({
  library(plotgardener)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
  library(org.Hs.eg.db)
  library(png)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  library(ggseqlogo)
  library(TFBSTools)
  library(JASPAR2020)
})

# 2025-12-16 Aurélien Ginolhac based on Déborah Gérard initial work
# Rendered using a container: `docker run -u $(id -u):$(id -g) -ti -v ${HOME}/Work/250603/Gerard_et_al_2026:/mnt 250603:24.04`
# Be sure that renv is not invoked if you use it: `mv .Rprofile Rprofile`
# Once in the container: cd /mnt/FIGURE3;  Rscript --vanilla Figure3.R

########## FIGURE 1 ##########

pdf(
  "FIGURE3.pdf",
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
  label = "Figure 3",
  fontsize = 14,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 0.25,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL A - ####

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

# Add the SNP position that zoom in the TFBS
SNP_param_SCARB2 <- pgParams(
  chrom = "chr4",
  chromstart = 76213717,
  chromend = 76213717,
  assembly = "hg38"
)

# ATACseq signal SCARB2
# path where bigwig files are
bw_path <- "../BIGWIG"

bw_samples <- c(
  "smNPC_I.bw",
  "smNPC_II.bw",
  "smNPC_III.bw",
  "D15_POS_I.bw",
  "D15_POS_II.bw",
  "D15_POS_III.bw",
  "D30_POS_I.bw",
  "D30_POS_II.bw",
  "D30_POS_III.bw",
  "HFFTHmCherry_D50_possort_S1.bw",
  "HFFTHmCherry_D50_possort_S3.bw",
  "HFFTHmCherry_D50_possort_S4.bw"
)

bw_SCARB2 <- lapply(bw_samples, \(x) {
  readBigwig(
    file = file.path(bw_path, x),
    chrom = "chr4",
    chromstart = 76173617,
    chromend = 76253817
  )
})


# Add a scale next to the bigwig files - check the maximum to choose
scale_SCARB2_max <- vapply(
  bw_SCARB2,
  \(x) max(x$score),
  FUN.VALUE = numeric(1)
) |>
  max() |>
  ceiling()

# Define parameters for the regions
region.p.SCARB2 <- pgParams(
  chrom = "chr4",
  chromstart = 76173717,
  chromend = 76253717,
  assembly = "hg38",
  range = c(0, scale_SCARB2_max)
)

#Add the genomic label
plotGenomeLabel(
  chrom = "chr4",
  chromstart = 76173717,
  chromend = 76253717,
  assembly = "hg38",
  x = 0.5,
  y = 0.6,
  length = 4.0,
  default.units = "inches",
  scale = "Mb"
)


NR2C2_SCARB2 <- vector(mode = "list", length = length(bw_SCARB2))

BW_COLORS <- c(
  rep("#A9A9A9", 3),
  rep("#B22222", 3),
  rep("#DC143C", 3),
  rep("#E9967A", 3)
)

Y_SCARB2 <- c(
  rep(1, 3),
  rep(1.3, 3),
  rep(1.6, 3),
  rep(1.9, 3)
)

SAMPLES <- c(
  "smNPC",
  "mDAN D15",
  "mDAN D30",
  "mDAN D50"
)


for (i in seq_along(bw_SCARB2)) {
  NR2C2_SCARB2[[i]] <- plotSignal(
    data = bw_SCARB2[[i]],
    params = region.p.SCARB2,
    fill = BW_COLORS[i],
    linecolor = NA,
    x = 0.5,
    y = Y_SCARB2[i],
    width = 4,
    height = 0.25,
    gapdistance = 0.1,
    just = c("left", "top"),
    default.units = "inches"
  )
  if (i %% 3 == 0) {
    plotText(
      label = SAMPLES[i / 3],
      fontsize = 6,
      fontcolor = BW_COLORS[i],
      fontfamily = "Helvetica",
      x = 0.25,
      y = Y_SCARB2[i] + 0.05,
      just = c("left", "top"),
      default.units = "inches",
      fontface = "bold"
    )
    annoYaxis(
      plot = NR2C2_SCARB2[[i]],
      at = c(0, scale_SCARB2_max),
      fontsize = 6
    )
  }
}


annoHighlight(
  plot = NR2C2_SCARB2[[1]],
  params = SNP_param_SCARB2,
  fill = "#404788FF",
  y = 1.1,
  height = 1.75,
  just = c("left", "top"),
  default.units = "inches"
)

plotGenes(
  chrom = "chr4",
  chromstart = 76173717,
  chromend = 76253717,
  assembly = assembly(
    Genome = "hg38refGene",
    TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene",
    OrgDb = "org.Hs.eg.db"
  ),
  fontcolor = "black",
  fill = "black",
  x = 0.5,
  y = 2.2,
  width = 4.0,
  height = 0.375,
  just = c("left", "top"),
  default.units = "inches"
)

# Add the transcription factor binding motif (TFBS) for NR2C2
pfm.NR2C2 <- getMatrixByID(JASPAR2020, ID = "MA1536.1")

# Reverse to position weight matrix as the SNP is on the reverse strand
pfm.NR2C2 <- reverseComplement(pfm.NR2C2)

NR2C2 <- new.env()

NR2C2$NR2C2 <- pfm.NR2C2@profileMatrix

NR2C2 <- as.list(NR2C2)

# Remove nucleotide
NR2C2_TFBS <- ggseqlogo(NR2C2) +
  annotate(
    "rect",
    xmin = 3.5,
    xmax = 4.5,
    ymin = -0.05,
    ymax = 2.0,
    alpha = .1,
    col = "black",
    fill = "#404788FF"
  ) +
  theme(axis.text.x = element_blank())

# Place the plot
plotGG(
  plot = NR2C2_TFBS,
  x = 1.5975,
  y = 2.6,
  width = 1.5,
  height = 1.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL B ###################################################################
# text B
plotText(
  label = "B",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 4.75,
  y = 0.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

dat.RPKM.filt <- read_rds("mat_RPKM.rds")

SCARB2.NR2C2.exp <- ggboxplot(
  dat.RPKM.filt |>
    dplyr::filter(gene_name %in% c("SCARB2", "NR2C2", "FAM47E")) |>
    mutate(
      gene_name = factor(gene_name, levels = c("SCARB2", "NR2C2", "FAM47E"))
    ),
  x = "Cond",
  y = "RPKM",
  add = "jitter",
  fill = "Cond",
  palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"),
  facet.by = "gene_name",
  xlab = "",
  ylab = "Reads Per Kilobase Million (RPKM)"
) +
  scale_x_discrete(labels = c("smNPC", "mDAN D15", "mDAN D30", "mDAN D50")) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


# Place the plot
plotGG(
  plot = SCARB2.NR2C2.exp,
  x = 5.0,
  y = 0.5,
  width = 3,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL C - #################################################################
# text C
plotText(
  label = "C",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.5,
  y = 4,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


gluc_seap <- tribble(
  ~`SCARB2\nother allele` , ~`SCARB2\neffect allele` ,
  1.11786288              , 0.82518605               ,
  1.10069041              , 0.83822192               ,
  1.11533243              , 0.99908577               ,
  1.31302042              , 1.12219885
) |>
  pivot_longer(cols = everything(), names_to = "sample")

# Plot
gluc_seap |>
  ggbarplot(
    x = "sample",
    y = "value",
    width = 0.3,
    fill = "sample",
    palette = c("grey40", "grey60"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = str_wrap("Gluc/SEAP normalized to positive control (mini CMV)", 30),
    title = "TH-REP1 neurons D11"
  ) +
  geom_hline(yintercept = 1, lty = "dashed", color = "black") +
  theme(
    axis.title.y = element_text(face = "bold", size = 9),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 8,
      angle = 45,
      vjust = 1,
      hjust = 1,
      family = "Helvetica"
    ),
    legend.position = "none",
    axis.text.y = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
  stat_compare_means(
    aes(label = paste0("p = ", after_stat(p.format))),
    method = "t.test",
    paired = TRUE,
    method.args = list(var.equal = TRUE),
    comparisons = list(c(1, 2))
  ) -> gluc_plot

plotGG(
  plot = gluc_plot,
  x = 0.5,
  y = 4.1,
  width = 1.7,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL D  ##################################################################

plotText(
  label = "D",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 4.2,
  y = 4,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

# Knockdown of NR2C2 and quantification of NR2C2 and SCARB2 expression 3 days post transduction
KD.SCARB2.NR2C2.3days <- tibble(
  Sple = rep("shNR2C2", 6),
  Gene = c(rep("NR2C2", 3), rep("SCARB2", 3)),
  Exp = c(
    0.345210918,
    0.512864636,
    0.548665969,
    1.239922544,
    1.130295504,
    1.170885498
  )
)

# Stats
KD.SCARB2.NR2C2.3days.stats <- tibble(
  Sple = c(rep("shNR2C2", 6), rep("shCTRL", 6)),
  Gene = rep(c(rep("NR2C2", 3), rep("SCARB2", 3)), 2),
  Exp = c(
    9.43575,
    9.3043,
    9.263,
    5.6841,
    5.92295,
    5.9289,
    7.9013,
    8.34095,
    8.397,
    5.99435,
    6.09965,
    6.1565
  )
)

stat_NR2C2.KD <- compare_means(
  Exp ~ Sple,
  data = KD.SCARB2.NR2C2.3days.stats,
  group.by = "Gene",
  method = "t.test",
  paired = TRUE
)

# Plot
pl.KD_SCARB2_NR2C2 <- KD.SCARB2.NR2C2.3days |>
  mutate(Gene = factor(Gene, levels = c("NR2C2", "SCARB2"))) |>
  ggbarplot(
    x = "Gene",
    y = "Exp",
    fill = c("Gene"),
    width = 0.4,
    palette = c("#440154FF", "#287C8EFF"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = "Relative expression to shSCRAMBLE",
    title = "shNR2C2 - 3 days"
  ) +
  geom_hline(yintercept = 1.0, lty = "dashed", color = "black") +
  theme(
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 10,
      angle = 90,
      family = "Helvetica"
    ),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
  stat_pvalue_manual(
    stat_NR2C2.KD,
    label = "p = {p.format}",
    x = "Gene",
    y.position = c(0.65, 1.5),
    size = 3.0,
    hide.ns = FALSE
  )

plotGG(
  plot = pl.KD_SCARB2_NR2C2,
  x = 4.3,
  y = 4.1,
  width = 1.7,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)

# Knockdown of NR2C2 and quantification of NR2C2 and SCARB2 expression 6 days post transduction
KD.SCARB2.NR2C2.6days <- tibble(
  Sple = rep("shNR2C2", 6),
  Gene = c(rep("NR2C2", 3), rep("SCARB2", 3)),
  Exp = c(
    0.662710186,
    0.456029798,
    0.787662775,
    1.410640142,
    1.425976296,
    1.332143962
  )
)

# Stats
KD.SCARB2.NR2C2.6days.stats = tibble(
  Sple = c(rep("shNR2C2", 6), rep("shCTRL", 6)),
  Gene = rep(c(rep("NR2C2", 3), rep("SCARB2", 3)), 2),
  Exp = c(
    8.34355,
    9.0917,
    8.32445,
    5.3595,
    5.19535,
    6.5266,
    7.75,
    7.9589,
    7.9801,
    5.85585,
    5.7073,
    6.94035
  )
)

stat_NR2C2.KD.6days <- compare_means(
  Exp ~ Sple,
  data = KD.SCARB2.NR2C2.6days.stats,
  group.by = "Gene",
  method = "t.test",
  paired = TRUE
)

# Plot
pl.KD_SCARB2_NR2C2_6days <- KD.SCARB2.NR2C2.6days |>
  mutate(Gene = factor(Gene, levels = c("NR2C2", "SCARB2"))) |>
  ggbarplot(
    x = "Gene",
    y = "Exp",
    fill = c("Gene"),
    width = 0.4,
    palette = c("#440154FF", "#287C8EFF"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = "Relative expression to shSCRAMBLE",
    title = "shNR2C2 - 6 days"
  ) +
  geom_hline(yintercept = 1.0, lty = "dashed", color = "black") +
  theme(
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 10,
      angle = 90,
      family = "Helvetica"
    ),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
  stat_pvalue_manual(
    stat_NR2C2.KD.6days,
    label = "p = {p.format}",
    x = "Gene",
    y.position = c(0.85, 1.5),
    size = 3.0,
    hide.ns = FALSE
  )

plotGG(
  plot = pl.KD_SCARB2_NR2C2_6days,
  x = 6.2,
  y = 4.1,
  width = 1.7,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL E - ##################################################

plotText(
  label = "E",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.5,
  y = 7.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

scarb2_degrouped <- read_tsv("scarb2_degrouped.tsv", show_col_types = FALSE)
snp_glmm <- read_tsv("snp_glmm.tsv", show_col_types = FALSE)

scarb2_degrouped |>
  ggplot(aes(x = sample, y = counts)) +
  geom_jitter(
    position = position_dodge(width = 0.45),
    alpha = 0.6,
    size = 3,
    aes(color = nucl, shape = "replicates")
  ) +
  stat_summary(
    fun.data = "mean_se",
    aes(color = nucl, shape = "mean + std error"),
    alpha = 0.6,
    position = position_dodge(width = 0.45)
  ) +
  geom_text(
    data = snp_glmm,
    aes(sample, x, label = scales::label_pvalue(add_p = TRUE)(padj)),
    nudge_y = 20,
    nudge_x = 0.4,
    size = 2.5
  ) +
  scale_colour_manual(values = c("red2", "grey60")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_shape_manual(values = c(19, 4)) +
  guides(
    color = guide_legend(
      override.aes = list(shape = c(16)),
      position = "inside"
    ),
    shape = guide_legend(position = "inside")
  ) +
  theme_classic(10) +
  theme(legend.position.inside = c(0.92, 0.87)) +
  labs(
    x = NULL,
    colour = "Allele",
    shape = NULL,
    caption = "adjusted pvalues (GLMM binomial)"
  ) -> scarb3_plot

plotGG(
  plot = scarb3_plot,
  x = 0.1,
  y = 7.65,
  width = 4.7,
  height = 3.9,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL F  ######################################################

plotText(
  label = "F",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 5.2,
  y = 7.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

GTEX_SNPs_interest <- read_csv(
  "20251219_GTEX_data.csv",
  show_col_types = FALSE
) |>
  mutate(
    GTEX_sig = if_else(
      P.Value < 0.05,
      "*",
      ""
    ),
    Tissue = str_remove(Tissue, "Brain_") |>
      str_replace_all('_', ' ') |>
      str_wrap(20)
  )


GTEX_eqtl <- ggplot(
  data = GTEX_SNPs_interest,
  mapping = aes(
    x = factor(
      SNP.Id,
      levels = c(
        "rs6825004",
        "rs7697073",
        "rs11547135",
        "rs6812193",
        "rs1465922"
      )
    ),
    y = Tissue,
    fill = NES
  )
) +
  geom_tile() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 7),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "eQTL \neffect size"
  ) +
  geom_text(aes(label = GTEX_sig), size = 6) +
  guides(
    fill = guide_colourbar(
      title.hjust = .5,
      barwidth = unit(.7, "lines"),
      barheight = unit(12, "lines")
    )
  )


plotGG(
  plot = GTEX_eqtl,
  x = 5.3,
  y = 7.6,
  width = 3,
  height = 4,
  just = c("left", "top"),
  default.units = "inches"
)

dev.off()
