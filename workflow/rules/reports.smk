shell.executable("bash")


# Plottings


rule host_prop:
    input:
        rpkm="{sample}/host_mapping/{sample}_host_rpkm.tsv",
    output:
        pdf="{sample}/reports/{sample}_hosts_detection.pdf",
        png="{sample}/reports/{sample}_hosts_detection.png",
        # html="{sample}/reports/{sample}_hosts_detection.html",
    params:
        sample=lambda wildcards: wildcards.sample
    message:
        "[{wildcards.sample}] plotting host genome repartition"
    log:
        "logs/{sample}/host_genome_plot.log"
    conda:
        "../envs/tidyr.yaml"
    script:
        "../scripts/host_treemap.R"


rule blast_graph:
    input:
        coverage="{sample}/assembly/{sample}_coverage_depth.tsv",
        elements="{sample}/blast/{sample}_blast.tsv",
        events="{sample}/blast/{sample}_events.tsv",
    output:
        plot_dir=directory("{sample}/reports/blast_graph"),
        placeholder=touch("{sample}/reports/blast_graph_artifact"),
        # html="{sample}/reports/{sample}_blast_graph.html",
    params:
        sample=lambda wildcards: wildcards.sample
    message:
        "[{wildcards.sample}] plotting BLAST results"
    log:
        "logs/{sample}/blast_graph.log"
    conda:
        "../envs/tidyr.yaml"
    script:
        "../scripts/blast_graph.R"