

## Manhattan plot

Produce as a PNG to be used in the figure panel.

``` r
# Associated then the significant GWAS pvalue
snp_2_high <- read_tsv("sneep_254_gwas_p_rsid.tsv", show_col_types = FALSE) |> 
  filter(!is.na(name),
                gwas_pvalue <= 5e-08) |> 
  select(chr,start, snp = name, gwas_pvalue)
  


snp_to_pl <- read_tsv("nalls_allSNPs_hg38.tsv.gz", show_col_types = FALSE) |>
  select(chr = seqnames, start, gwas_pvalue = p) |> 
  bind_rows(snp_2_high) |>
  mutate(p = -log10(gwas_pvalue),
         chr = factor(chr,
                        levels = paste0("chr",
                                        1:22)),
         color = if_else(is.na(snp), "no_high", "high"),
         label = if_else(snp %in% c("rs1465922", "rs144814361"), snp, NA_character_))

highlight_colormap <- c("no_high" = adjustcolor( "grey", 
                                                 alpha.f = 0.2), 
                        "high" = "#6600FF")

no_high <- manhattan_data_preprocess(snp_to_pl,
                                     pval.colname = "gwas_pvalue",
                                     chr.colname = "chr",
                                     pos.colname = "start",
                                     highlight.colname = "color",
                                     highlight.col = highlight_colormap,
                                     signif = 5e-08)
# Plot
snp.pl <- manhattan_plot(x = no_high,
                         color.by.highlight = TRUE,
                         rescale = TRUE,
                         label.font.size = 2.5,
                         label.colname = "label") +
  theme(axis.text.x = element_text(size = 7, hjust = 1, vjust = 0.5,
                                   angle = 90),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9))

ggsave("manhattan_plot.png", plot = snp.pl, dpi = "print", height = 2, width = 3.75)
```