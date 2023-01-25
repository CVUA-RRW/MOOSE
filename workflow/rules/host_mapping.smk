shell.executable("bash")


rule map_seq2org:
    output:
        map="common/seqid_mapping.txt",
    params:
        refs=config["indexed_ref"],
    message:
        "[all] mapping seqIDs to organisms"
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


rule map_reads_host:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        aln="{sample}/host_mapping/{sample}_host_aln.bam",
    params:
        refs=config["indexed_ref"],
    threads:
        config['threads_sample']
    message:
        "[{wildcards.sample}] mapping reads to host genomes with BWA"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/mapping_host.log"
    shell:
        """
        exec 2> {log}
        bwa mem -t {threads} \
            {params.refs} {input.r1} {input.r2} \
            | samtools sort -@{threads} -o {output.aln} 
        """


rule mapping_stats_host:
    input:
        aln="{sample}/host_mapping/{sample}_host_aln.bam",
    output:
        flagstats="{sample}/host_mapping/{sample}_host_flagstats.tsv",
        mapping_stats="{sample}/host_mapping/{sample}_host_mapstats.tsv",
    message:
        "[{wildcards.sample}] retrieving host mapping statistics"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/host_mapping_stats.log"
    shell:
        """
        exec 2> {log}
        samtools flagstat -O tsv {input.aln} > {output.flagstats}
        samtools stats {input.aln} > {output.mapping_stats}
        """


rule coverage_host:
    input:
        aln="{sample}/host_mapping/{sample}_host_aln.bam",
    output:
        cov="{sample}/host_mapping/{sample}_host_coverage.tsv",
    message:
        "[{wildcards.sample}] calculating host coverage"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/host_coverage.log"
    shell:
        """
        exec 2> {log}
        # Coverage and calculate RPK and add header
        samtools coverage {input.aln} > {output.cov}
        """


rule host_rpkm:
    input:
        cov="{sample}/host_mapping/{sample}_host_coverage.tsv",
        mapping_stats="{sample}/host_mapping/{sample}_host_mapstats.tsv",
        map="common/seqid_mapping.txt",
    output:
        rpkm="{sample}/host_mapping/{sample}_host_rpkm.tsv",
    message:
        "[{wildcards.sample}] calculating RPKM per host"
    log:
        "logs/{sample}/rpkm.log"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/host_reads_to_rpkm.py"