#!/usr/bin/env Rscript

# logging
log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')

# Imports
library(tidyverse)
library(RColorBrewer)
library(cowplot)
# library(plotly)
# library(htmlwidgets)

# Get snakemake parameters
coverage <- snakemake@input[['coverage']]
elements <- snakemake@input[['elements']]
events <- snakemake@input[['events']]
out_folder <- snakemake@output[['plot_dir']]
# out_html <- snakemake@output[['html']]
sample_name <- snakemake@params[['sample']]

# Aestethics
id_breaks <- c(0,60,70,80,90,100)
id_labels <- c("<60%", ">60%", ">70%", ">80%", ">90%")
color_scale <- c(
    "<60%"='#d7191c',
    ">60%"='#fdae61',
    ">70%"='#ffffbf',
    ">80%"='#a6d96a', 
    ">90%"='#1a9641'
    )
    

dir.create(out_folder)

# Load data and format
covdepth <- 
    read_tsv(
        coverage, 
        col_names=c('contig_id', 'pos', 'depth'),
        col_types="cii",
        skip=1
    )

belements <- 
    read_tsv(elements) %>%
    select(
        query_id, 
        subject_id, 
        query_start,
        query_end,
        subjectlength,
        subject_start,
        subject_end,
        strand,
        identity
    ) %>%
    rename(
        contig_id=query_id,
        seq_id=subject_id, 
        contig_start=query_start,
        contig_end=query_end,
        seq_length=subjectlength,
        seq_start=subject_start,
        seq_end=subject_end)
        
bevents <- 
    read_tsv(events) %>%
    select(
        subject_id,
        query_id,
        subject_start,
        subject_end,
        query_length,
        query_start,
        query_end,
        strand,
        identity
    ) %>%
    rename(
        contig_id=subject_id,
        seq_id=query_id,
        contig_start=subject_start,
        contig_end=subject_end,
        seq_length=query_length,
        seq_start=query_start,
        seq_end=query_end
    )

contig_lengths <-
    covdepth %>% 
    group_by(contig_id) %>% 
    slice_max(pos) %>% 
    select(-depth) %>% 
    rename(contig_length=pos)

blast <- bind_rows(belements, bevents)

contigs <- 
    covdepth %>%
    pull('contig_id') %>%
    unique()

# map over contigs
walk(contigs, ~{
    total_length <- 
        filter(contig_lengths, contig_id==.x) %>%
        pull(contig_length)
    
    # Prepping blast table for tiles
    blast_sorted <-
        blast %>%
        filter(contig_id==.x) %>%
        arrange(desc(identity), seq_id) %>%
        rowwise() %>%
        # get mapped parts width and center for geom_tiles()
        mutate(
            width=abs(seq_end-seq_start),
            center=min(contig_start, contig_end)+width/2
        ) %>%
        # get the row number for arranging on y axis
        rowid_to_column() %>%
        # discretize identity level
        mutate(identity_interval=cut_interval(
                                    identity, 
                                    n=5, 
                                    breaks=id_breaks,
                                    labels=id_labels)
        ) %>%
        # get the start and stop of the query sequence relative to the contig
        # and add stuffs for labels
        mutate(
            seq_start_contig=if_else(strand=="plus", contig_start-seq_start, contig_start-(seq_length-seq_end)),
            seq_end_contig=if_else(strand=="plus", contig_end+(seq_length-seq_end), contig_end+seq_start),
            direction=if_else(strand=="plus", '(+)', '(-)')
        )
    
    # Contig map
    contig_map <-
        ggplot(blast_sorted) +
        geom_segment(
            aes(x=1, xend=total_length, y=0, yend=0),
            # arrow=arrow(type='closed'),
            lineend='butt',
            # linejoin='mitre',
            show.legend=FALSE,
            inherit.aes=FALSE,
            linewidth=2
        ) +
        geom_tile(
            aes(x=center, y=0, width=width, height=0.6),
            fill='#5e3c99',
            linewidth=0,
            show.legend=FALSE,
            linewidth=1) +
        theme_void() +
        theme(legend.position="none")
    
    # Coverage density
    depth_p <-
        covdepth %>%
        filter(contig_id==.x) %>%
        arrange(pos) %>%
        ggplot(aes(x=pos, y=depth)) +
        geom_bar(
            show.legend=FALSE, 
            width=1, 
            stat='identity',
            fill='#2b83ba',
            color='#2b83ba') +
        theme_classic() +
        ylab("Depth (number of reads)") +
        theme(
            axis.line.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            legend.position="none"
        )
    
    # Blast 
    blast_map <-
        ggplot(blast_sorted) +
        geom_segment(
            aes(
                x=seq_start_contig,
                xend=seq_end_contig,
                y=rowid,
                yend=rowid
            ), 
            lineend='butt',
            show.legend=FALSE,
            linewidth=2
        ) +
        geom_tile(
            aes(
                x=center,
                y=rowid,
                width=width,
                height=0.6,
                fill=identity_interval
            ),
            linewidth=0,
            linewidth=1
        ) +
        geom_text(
            aes(x=center, y=rowid, label=paste(direction, seq_id)),
            hjust='inward',
            nudge_y=0.4
        ) +
        labs(fill="Identity [%]") +
        ylab("Local sequence alignement") +
        xlab("Position") +
        scale_fill_manual(values=color_scale, drop=FALSE) +
        theme_classic() +
        theme(
            legend.position="none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.y=element_blank()
            ) +
        xlim(0,total_length)

    # merge graphs
    n_blast <- 
        blast_sorted %>% 
        pull(rowid) %>%
        max()

    main <- 
        plot_grid(
            depth_p, contig_map, blast_map, 
            rel_heights=c(2, 1, n_blast),
            labels="",
            ncol=1,
            align="v")
    
    # Add title
    title <- 
        ggdraw() + 
        draw_label(
            paste0(sample_name, "\n", .x),
            fontface = 'bold',
            x = 0,
            hjust = 0
        ) +
        theme(
            # add margin on the left of the drawing canvas,
            # so title is aligned with left edge of first plot
            plot.margin = margin(0, 0, 0, 7)
        )
    
    # Add legend for blast ID
    legend_b <- 
        get_legend(
            blast_map +
            guides(color=guide_legend(nrow=1)) +
            theme(
                legend.position='bottom',
                legend.key.height= unit(0.5, 'cm'),
                legend.key.width= unit(1, 'cm'))
        )
    
    # Make figure
    plot_grid(
        title, main, legend_b,
        ncol = 1,
        rel_heights = c(0.1, 1, 0.1)
    )
    
    ggsave(
        file.path(out_folder, paste0(.x, ".pdf")),
        width =18,
        units='cm'
    )
    ggsave(
        file.path(out_folder, paste0(.x, ".png")),
        width= 18,
        units='cm'
    )
})

