shell.executable("bash")


rule map_reads_host:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        aln="{sample}/host_mapping/{sample}_aln.sam",
    params:
        refs=config["indexed_ref"],
    threads:
        config['threads_sample']
    message:
        "Mapping {wildcards.sample} to host genomes"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/mapping_host.log"
    shell:
        """
        exec 2> {log}
        bwa mem -t {threads} \
            {params.refs} {input.r1} {input.r2} > {output.aln}
        """


rule mapping_stats_host:
    input:
        sam="{sample}/host_mapping/{sample}_aln.sam",
    output:
        flagstats="{sample}/host_mapping/{sample}_flagstats.tsv",
        mapping_stats="{sample}/host_mapping/{sample}_mapstats.tsv",
        unmapped_sam="{sample}/host_mapping/{sample}_unmapped.sam",
        unmapped_stats="{sample}/host_mapping/{sample}_unmapstats.tsv",
    message:
        "Retrieving host_mapping stats for {wildcards.sample}"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{sample}/host_mapping_stats.log"
    shell:
        """
        exec 2> {log}
        samtools flagstat -O tsv {input.sam} > {output.flagstats}
        samtools stats {input.sam} > {output.mapping_stats}
        # -F filtering UNMAP and MATE UNMAP reads (aka not mapped in proper pair)
        samtools view -h -F 2 {input.sam} > {output.unmapped_sam}
        samtools stats {output.unmapped_sam} > {output.unmapped_stats}
        """


rule map_seq2org:
    output:
        map="common/seqid_mapping.txt",
    params:
        refs=config["indexed_ref"],
    message:
        "Mapping seqIDs to organisms"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/common/seqid_mapping.log"
    shell:
        """
        exec 2> {log}
        # make table seqid org
        grep ">" {params.refs} | tr -d ">" | cut -f1,2,3 -d" " | sort -k1 > {output.map}
        """


rule mapped_count_host:
    input:
        sam="{sample}/host_mapping/{sample}_aln.sam",
        map="common/seqid_mapping.txt",
    output:
        perseq="{sample}/host_mapping/{sample}_read_per_seq.txt",
        perorg="{sample}/host_mapping/{sample}_read_per_org_long.txt",
    message:
        "Counting reads per organism for {wildcards.sample}"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{sample}/count_reads.log"
    shell:
        """
        # exec 2> {log}
        # get count per seq
        samtools view -f 2 -F 104 {input.sam} | cut -f3 | sort | uniq -c | sort -k2 > {output.perseq}
        # Count per org
        join -1 1 -2 2 -o 1.2,1.3,2.1 {input.map} {output.perseq} > {output.perorg}
        """


rule collapse_counts_host:
    input:
        counts="{sample}/host_mapping/{sample}_read_per_org_long.txt",
    output:
        counts="{sample}/host_mapping/{sample}_read_per_org.tsv",
    message:
        "Collapsing read count per organism for {wildcards.sample}"
    log:
        "logs/{sample}/collapse_counts.log"
    run:
        import sys
        sys.stderr = open(log[0], "w")
        
        import pandas as pd
        
        df = pd.read_csv(input.counts, sep=" ", names=['genus', 'species', 'count'])
        df=df.groupby(['genus', 'species']).sum().reset_index().sort_values('count', ascending=False)
        df['organism'] = df['genus'].map(str) + ' ' + df['species'].map(str)
        df[['organism', 'count']].to_csv(output.counts, sep="\t", index=False)


    