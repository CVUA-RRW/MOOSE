#!/usr/bin/env Rscript


# logging
log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

# Imports
library(tidyverse)
library(RColorBrewer)
library(treemapify)
# library(plotly)
# library(htmlwidgets)

# Get snakemake parameters
input <- snakemake@input[['rpkm']]
out_pdf <- snakemake@output[['pdf']]
out_png <- snakemake@output[['png']]
out_html <- snakemake@output[['html']]
sample_name <- snakemake@params[['sample']]

# Load and plot data
rpkm <- read_tsv(input)
rpkm_p <- 
    ggplot(
        rpkm,
        aes(
            area=RPKM, 
            fill=organism, 
            label=organism
            )
        ) +
    geom_treemap() +
    geom_treemap_text(
        fontface='italic', 
        colour="black", 
        place="center", 
        grow=TRUE
        ) +
    theme_void() +
    scale_color_brewer(palette='Set1') +
    theme(legend.position='none') +
    labs(
        title=paste0("RPKM repartition over host genomes for ", sample_name),
        caption="Box areas are proportional to the RPKM (reads per kilobase per million) value for each host genome."
        )

ggsave(out_pdf)
ggsave(out_png)

# Interactive html chart

# rpkm_ply <-
    # ggplotly(
        # p=rpkm_p,
        # width=600,
        # height=600,
        # tooltip=c('organism', 'RPKM')
        # )

# saveWidget(rpkm_ply, file=out_html)
