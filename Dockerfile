FROM rocker/r2u:24.04

# no need to install system libraries
RUN install.r  plotgardener BiocManager tidyverse \
  ggpubr ggseqlogo TFBSTools
RUN installBioc.r TxDb.Hsapiens.UCSC.hg38.knownGene \
  TxDb.Hsapiens.UCSC.hg38.refGene org.Hs.eg.db JASPAR2020 \
  AllelicImbalance
RUN installBioc.r ggmanh
# docker buildx build -t 250603:24.04 .
