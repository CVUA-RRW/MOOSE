shell.executable("bash")


# reports


rule assembly_report:
    input:
        all_contigs="{sample}/assembly/{sample}_contigs.fasta",
        clus_contigs="{sample}/assembly/{sample}_collapsed_contigs.fasta",
        fil_contigs="{sample}/reports/{sample}_coverage_filtered.tsv",
        flagstats_before="{sample}/assembly/{sample}_flagstats_unfiltered.tsv",
        flagstats_after="{sample}/assembly/{sample}_flagstats_filtered.tsv",
    output:
        report="{sample}/reports/{sample}_assembly_report.tsv",
    params:
        sample=lambda w: w.sample,
    message:
        "[{wildcards.sample}] reporting on assembly"
    log:
        "logs/{sample}/assembly_report.tsv",
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        exec 2> {log}
        echo "sample\ttotal reads\ttotal contigs\tcollapsed contigs\tfiltered contigs\tmean contig length\tmean coverage depth\tmean reads in contigs\tprimary mapped all\t primarymapped all [%]\t primary mapped filtered\tprimary mapped filteed [%]" \
            > {output.report}
        treads=$(grep "total" {input.flagstats_before} | cut -d"\t" -f1)
        total=$(grep -c ">" {input.all_contigs})
        clustered=$(grep -c ">" {input.clus_contigs})
        filtered=$(tail -n +2 {input.fil_contigs} | wc -l | cut -d" " -f1)
        length=$(awk '{{ sum += $4 }} END {{ if (NR > 0) print sum / NR }}' {input.fil_contigs} || 0)
        depth=$(awk '{{ sum += $8 }} END {{ if (NR > 0) print sum / NR }}' {input.fil_contigs} || 0)
        nreads=$(awk '{{ sum += $5 }} END {{ if (NR > 0) print sum / NR }}' {input.fil_contigs} || 0)
        map_before=$(grep "primary mapped$" {input.flagstats_before} | cut -d"\t" -f1)
        map_before_p=$(grep "primary mapped %$" {input.flagstats_before} | cut -d"\t" -f1 | tr -d "%")
        map_after=$(grep "primary mapped$" {input.flagstats_after} | cut -d"\t" -f1)
        map_after_p=$(echo "scale=2; 100*$map_after/$treads" | bc)
        echo "{params.sample}\t$treads\t$total\t$clustered\t$filtered\t$length\t$depth\t$nreads\t$map_before\t$map_before_p\t$map_after\t$map_after_p" \
            >> {output.report}
        """

 
rule gather_assembly_reports:
    input:
        report=expand("{sample}/reports/{sample}_assembly_report.tsv", sample=samples.index),
    output:
        agg="reports/assembly_report.tsv"
    message:
        "[all] Aggregating assembly stats"
    log:
        "logs/report_assembly.log"
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """       


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
        # events="{sample}/blast/{sample}_events.tsv",
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
