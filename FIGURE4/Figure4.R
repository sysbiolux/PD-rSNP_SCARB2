suppressPackageStartupMessages({
  library(plotgardener)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
  library(org.Hs.eg.db)
  library(png)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
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
# Once in the container: cd /mnt/FIGURE4;  Rscript --vanilla Figure4.R

########## FIGURE 1 ##########

pdf(
  "FIGURE4.pdf",
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
  label = "Figure 4",
  fontsize = 14,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 0.25,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL A - BAG3 locus ######################################################

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


# ATACseq signal BAG3
# path where bigwig files are
bw_path <- "../BIGWIG"


bw_samples <- c(
  "smNPC_I.bw",
  "smNPC_II.bw",
  "smNPC_III.bw",
  "D15_POS_II.bw",
  "D15_POS_III.bw",
  "D30_POS_I.bw",
  "D30_POS_II.bw",
  "D30_POS_III.bw",
  "HFFTHmCherry_D50_possort_S1.bw",
  "HFFTHmCherry_D50_possort_S3.bw",
  "HFFTHmCherry_D50_possort_S4.bw"
)

bw_BAG3 <- lapply(bw_samples, \(x) {
  readBigwig(
    file = file.path(bw_path, x),
    chrom = "chr10",
    chromstart = 119611305,
    chromend = 119691505
  )
})

# Add a scale next to the bigwig files - check the maximum to choose
scale_BAG3_max <- vapply(
  bw_BAG3,
  \(x) max(x$score),
  FUN.VALUE = numeric(1)
) |>
  max() |>
  ceiling()

# Define parameters for the regions
region.p.BAG3 = pgParams(
  chrom = "chr10",
  chromstart = 119611405,
  chromend = 119691405,
  assembly = "hg38",
  range = c(0, scale_BAG3_max)
)

# Add the genomic label
plotGenomeLabel(
  chrom = "chr10",
  chromstart = 119611405,
  chromend = 119691405,
  assembly = "hg38",
  x = 0.5,
  y = 0.6,
  length = 4.0,
  default.units = "inches",
  scale = "Mb"
)

BAG3 <- vector(mode = "list", length = length(bw_BAG3))
BW_COLORS <- c(
  rep("#A9A9A9", 3),
  rep("#B22222", 3),
  rep("#DC143C", 3),
  rep("#E9967A", 3)
)
Y_BAG3 <- c(
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

for (i in seq_along(bw_BAG3)) {
  BAG3[[i]] <- plotSignal(
    data = bw_BAG3[[i]],
    params = region.p.BAG3,
    fill = BW_COLORS[i],
    linecolor = NA,
    x = 0.5,
    y = Y_BAG3[i],
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
      y = Y_BAG3[i] + 0.05,
      just = c("left", "top"),
      default.units = "inches",
      fontface = "bold"
    )
    annoYaxis(
      plot = BAG3[[i]],
      at = c(0, scale_BAG3_max),
      fontsize = 6
    )
  }
}


# Gene track
plotGenes(
  chrom = "chr10",
  chromstart = 119611405,
  chromend = 119691405,
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

# Add the SNP position that zoom in the TFBS
SNP_param <- pgParams(
  chrom = "chr10",
  chromstart = 119651405,
  chromend = 119651405,
  assembly = "hg38"
)

annoHighlight(
  plot = BAG3[[1]],
  params = SNP_param,
  fill = "#404788FF",
  y = 1.1,
  height = 1.5,
  just = c("left", "top"),
  default.units = "inches"
)

# Add the transcription factor binding motif (TFBS) for LHX1
pfm.LHX1 <- getMatrixByID(JASPAR2020, ID = "MA1518.1")
LHX1 <- new.env()

LHX1$LHX1 = pfm.LHX1@profileMatrix

LHX1 <- as.list(LHX1)

# Remove the nucleotide position
LHX1_TFBS = ggseqlogo(LHX1) +
  annotate(
    "rect",
    xmin = 4.5,
    xmax = 5.5,
    ymin = -0.05,
    ymax = 2.5,
    alpha = .1,
    col = "black",
    fill = "#404788FF"
  ) +
  theme(axis.text.x = element_blank())

# Place the plot
plotGG(
  plot = LHX1_TFBS,
  x = 1.5975,
  y = 2.6,
  width = 1.5,
  height = 1.5,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL B - BAG3 and LHX1 expression (RNAseq)     ###########################

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

dat.RPKM.filt <- read_rds("../FIGURE3/mat_RPKM.rds")

# Expression for BAG3 and LHX1
BAG3.LHX1.exp <- ggboxplot(
  dat.RPKM.filt |>
    dplyr::filter(gene_name %in% c("BAG3", "LHX1")),
  x = "Cond",
  y = "RPKM",
  add = "jitter",
  fill = "Cond",
  palette = c("#A9A9A9", "#B22222", "#DC143C", "#E9967A"),
  facet.by = "gene_name",
  xlab = "",
  ylab = "Reads Per Kilobase Million (RPKM)"
) +
  scale_x_discrete(
    labels = c(
      "smNPC",
      expression("mDAN D15"),
      expression("mDAN D30"),
      expression("mDAN D50")
    )
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )


# Place the plot
plotGG(
  plot = BAG3.LHX1.exp,
  x = 5.0,
  y = 0.5,
  width = 3,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL C - ##################################################
# text C
plotText(
  label = "C",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 4.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


gluc_seap <- tribble(
  ~`BAG3\nother allele` , ~`BAG3\neffect allele` ,
  1.330145828           , 0.941542666            ,
  1.24940239            , 0.89488201             ,
  1.52514654            , 0.98599237
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
    title = ""
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
    size = 3,
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


#### PANEL D  ######################################################

plotText(
  label = "D",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 3,
  y = 4.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


# Knockdown of LHX1 and quantification of LHX1 and BAG3 expression 3 days post transduction
KD.BAG3.LHX1.3days <- tibble(
  Sple = rep("shLHX1", 6),
  Gene = c(rep("LHX1", 3), rep("BAG3", 3)),
  Exp = c(
    0.411253539,
    0.359459217,
    0.367190337,
    0.972688658,
    0.965902853,
    0.88813445
  )
)

# Stats
KD.BAG3.LHX1.3days.stats <- tibble(
  Sple = c(rep("shLHX1", 6), rep("shCTRL", 6)),
  Gene = rep(c(rep("LHX1", 3), rep("BAG3", 3)), 2),
  Exp = c(
    9.5167,
    8.39675,
    8.1855,
    9.13605,
    8.33795,
    8.76415,
    8.2348,
    6.92065,
    6.7401,
    9.0961,
    8.2879,
    8.593
  )
)

stat_LHX1.KD <- compare_means(
  Exp ~ Sple,
  data = KD.BAG3.LHX1.3days.stats,
  group.by = "Gene",
  method = "t.test",
  paired = TRUE
)
# Plot
pl.KD_BAG3_LHX1 <- KD.BAG3.LHX1.3days |>
  mutate(Gene = factor(Gene, levels = c("LHX1", "BAG3"))) %>%
  ggbarplot(
    x = "Gene",
    y = "Exp",
    fill = c("Gene"),
    position = position_dodge(0.9),
    width = 0.4,
    palette = c("#440154FF", "#287C8EFF"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = "Relative expression to shSCRAMBLE",
    title = "shLHX1 - 3 days"
  ) +
  geom_hline(yintercept = 1.0, lty = "dashed", color = "black") +
  theme(
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 10,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      family = "Helvetica"
    ),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
  stat_pvalue_manual(
    stat_LHX1.KD,
    label = "p = {p.format}",
    x = "Gene",
    y.position = c(0.5, 1.2),
    size = 3.0,
    hide.ns = FALSE
  )

plotGG(
  plot = pl.KD_BAG3_LHX1,
  x = 3.4,
  y = 4.2,
  width = 1.8,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL E  ######################################################

plotText(
  label = "E",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 5.25,
  y = 4.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

tribble(
  ~`ALLELE C_BAG3` , ~`C/C`   , ~`C/T`  , ~NTC     ,
  2.46657          , 0.621758 , NA      , NA       ,
  2.48847          , 0.545051 , NA      , NA       ,
  2.34544          , 0.500021 , NA      , NA       ,
  2.32725          , 0.499881 , NA      , NA       ,
  2.29662          , 0.505719 , NA      , NA       ,
  2.28812          , 0.493492 , NA      , NA       ,
  2.36095          , 0.508781 , NA      , NA       ,
  2.37598          , 0.515607 , NA      , NA       ,
  2.29781          , 0.501133 , NA      , NA       ,
  2.42384          , 0.572841 , NA      , NA       ,
  2.31259          , 0.510102 , NA      , NA       ,
  2.20307          , 0.566714 , NA      , NA       ,
  2.38593          , 0.511483 , NA      , NA       ,
  2.39012          , 0.502775 , NA      , NA       ,
  2.2386           , 0.503122 , NA      , NA       ,
  2.29505          , 0.48322  , NA      , NA       ,
  2.10849          , 0.448639 , NA      , NA       ,
  2.06545          , 0.466778 , NA      , NA       ,
  2.15762          , 0.465845 , NA      , NA       ,
  1.95428          , NA       , 1.21359 , NA       ,
  2.36055          , 0.503087 , NA      , NA       ,
  2.47763          , 0.528391 , NA      , NA       ,
  2.28943          , 0.490323 , NA      , NA       ,
  2.43726          , 0.524907 , NA      , NA       ,
  2.47414          , 0.529655 , NA      , NA       ,
  2.46027          , 0.510073 , NA      , NA       ,
  2.11818          , 0.478299 , NA      , NA       ,
  2.38019          , 0.471852 , NA      , NA       ,
  2.30632          , 0.482912 , NA      , NA       ,
  2.11668          , 0.443099 , NA      , NA       ,
  2.36308          , 0.470144 , NA      , NA       ,
  2.21907          , 0.471723 , NA      , NA       ,
  2.38006          , 0.529287 , NA      , NA       ,
  2.40473          , 0.514012 , NA      , NA       ,
  0.477947         , NA       , NA      , 0.392688 ,
  0.478545         , NA       , NA      , 0.383172 ,
  0.473244         , NA       , NA      , 0.360477
) |>
  pivot_longer(
    cols = -`ALLELE C_BAG3`,
    names_to = "name",
    values_to = "ALLELE T_BAG3"
  ) |>
  drop_na(`ALLELE T_BAG3`) |>
  ggplot(aes(x = `ALLELE C_BAG3`, y = `ALLELE T_BAG3`, colour = name)) +
  geom_point(alpha = 0.6) +
  theme_classic(8) +
  guides(colour = guide_legend(position = "inside")) +
  theme(
    legend.position.inside = c(0.3, 0.7),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(title = "Allelic discremination", colour = NULL) -> allele_BAG3


plotGG(
  plot = allele_BAG3,
  x = 5.8,
  y = 4.2,
  width = 1.8,
  height = 1.8,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL F  ######################################################

plotText(
  label = "F",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 5.25,
  y = 6,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


#### PANEL G  ######################################################

plotText(
  label = "G",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 0.25,
  y = 7.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)

# Knockdown of LHX1 and quantification of LHX1 and BAG3 expression 2 days post transduction
KD.BAG3.LHX1.2days <- tibble(
  Sple = "shLHX1",
  Gene = c(rep("LHX1", 4), rep("BAG3", 4)),
  Exp = c(
    0.327601,
    0.294289,
    0.693776,
    0.57871,
    1.217409,
    0.784028,
    1.374824,
    0.814344
  )
)

# One sample t-test
stat.1.test <- KD.BAG3.LHX1.2days |>
  group_by(Gene) |>
  rstatix::t_test(Exp ~ 1, mu = 1) |>
  rstatix::adjust_pvalue()

# Plot
pl.KD_BAG3_LHX1_2 <- KD.BAG3.LHX1.2days |>
  mutate(Gene = factor(Gene, levels = c("LHX1", "BAG3"))) %>%
  ggbarplot(
    x = "Gene",
    y = "Exp",
    fill = c("Gene"),
    width = 0.4,
    palette = c("#440154FF", "#287C8EFF"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = "Relative expression\nto shSCRAMBLE",
    title = "shLHX1 - 2 days"
  ) +
  geom_hline(yintercept = 1, lty = "dashed", color = "black") +
  theme(
    axis.title.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(
      face = "bold.italic",
      size = 10,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      family = "Helvetica"
    ),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.8)) +
  stat_pvalue_manual(
    stat.1.test,
    label = "p = {p.adj}",
    size = 3,
    y.position = c(1.5, 0.8),
    xmin = "Gene",
    xmax = NULL
  )

plotGG(
  plot = pl.KD_BAG3_LHX1_2,
  x = 0.2,
  y = 7.9,
  width = 2,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#### PANEL H  ######################################################

plotText(
  label = "H",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 2.75,
  y = 7.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


gluc_seap_smNPC <- tribble(
  ~`BAG3\nother allele` , ~`BAG3\neffect allele` ,
  0.832108              , 0.738428               ,
  1.060865              , 0.49314                ,
  0.896677              , 0.797812
) |>
  pivot_longer(cols = everything(), names_to = "sample")


# Plot
gluc_seap_smNPC |>
  ggbarplot(
    x = "sample",
    y = "value",
    width = 0.3,
    fill = "sample",
    palette = c("grey40", "grey60"),
    add = c("mean_se", "jitter"),
    xlab = "",
    ylab = str_wrap("Gluc/SEAP normalized to positive control (mini CMV)", 30),
    title = "BAG3 Luciferase in smNPC"
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
    size = 3,
    method.args = list(var.equal = TRUE),
    comparisons = list(c(1, 2))
  ) -> gluc_smNPC_plot

plotGG(
  plot = gluc_smNPC_plot,
  x = 2.75,
  y = 7.9,
  width = 1.7,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)

#### PANEL I  ######################################################

plotText(
  label = "I",
  fontsize = 12,
  fontfamily = "Helvetica",
  x = 5.25,
  y = 7.5,
  just = "left",
  default.units = "inches",
  fontface = "bold"
)


BAG3_wider <- read_tsv("BAG3_wider.tsv", show_col_types = FALSE)
glmm_bag3 <- read_tsv("glmm_bag3.tsv", show_col_types = FALSE)

BAG3_wider |>
  pivot_longer(cols = C:T, names_to = "nucl", values_to = "counts") |>
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
    data = glmm_bag3,
    aes(sample, x, label = scales::label_pvalue(add_p = TRUE)(p.value)),
    nudge_y = 0,
    nudge_x = 0.1,
    size = 3
  ) +
  scale_colour_manual(values = c("red2", "grey60")) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(),
    breaks = c(0, 5, 10, 25, 50, 100, 200)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  scale_shape_manual(values = c(19, 4)) +
  guides(
    color = guide_legend(
      override.aes = list(shape = c(16)),
      position = "inside"
    ),
    shape = "none"
  ) +
  theme_classic(10) +
  theme(
    legend.position.inside = c(0.5, 0.8),
    plot.title = element_text(face = "bold", size = 9)
  ) +
  labs(
    x = NULL,
    colour = "Allele",
    y = "counts (pseudo log)",
    title = "TH-REP1 with SNP-BAG3",
    shape = NULL,
    caption = "adj pvalue GLMM binomial"
  ) -> bag3_plot

plotGG(
  plot = bag3_plot,
  x = 5.2,
  y = 7.9,
  width = 2.9,
  height = 3.5,
  just = c("left", "top"),
  default.units = "inches"
)


#pageGuideHide()
dev.off()
