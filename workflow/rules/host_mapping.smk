shell.executable("bash")


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


rule map_reads_host:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        aln="{sample}/host_mapping/{sample}_aln.bam",
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
            {params.refs} {input.r1} {input.r2} \
            | samtools sort -@{threads} -o {output.aln} 
        """


rule mapping_stats_host:
    input:
        aln="{sample}/host_mapping/{sample}_aln.bam",
    output:
        flagstats="{sample}/host_mapping/{sample}_flagstats.tsv",
        mapping_stats="{sample}/host_mapping/{sample}_mapstats.tsv",
    message:
        "Retrieving host_mapping stats for {wildcards.sample}"
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
        aln="{sample}/host_mapping/{sample}_aln.bam",
    output:
        cov="{sample}/host_mapping/{sample}_hosts_coverage.tsv",
    message:
        "Calculating host coverage for {wildcards.sample}"
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


rule normalized_read_count:
    input:
        cov="{sample}/host_mapping/{sample}_hosts_coverage.tsv",
        mapping_stats="{sample}/host_mapping/{sample}_mapstats.tsv",
        map="common/seqid_mapping.txt",
    output:
        rpkm="{sample}/host_mapping/{sample}_rpkm.tsv",
    message:
        "Calculating RPKM per host for {wildcards.sample}"
    log:
        "logs/{sample}/rpkm.log"
    run:
        import sys
        sys.stderr = open(log[0], "w")
        
        import pandas as pd
        
        # get total reads
        with open(input.mapping_stats, 'r') as fi:
            for line in fi:
                l = line.split("\t")
                if l[0] == "SN" and l[1] == "reads mapped:":
                    miReads = int(l[2])/1000000
        
        # Adding organism info to coverage table, then collapse and finally calculate rpkm
        map = pd.read_csv(input.map, sep=" ", names=['id', 'genus', 'species'])
        cov = pd.read_csv(input.cov, sep="\t"
            ).merge(
                map, left_on="#rname", right_on='id', how='left'
            ).groupby(
                ['genus', 'species']
            ).sum().reset_index()
        cov['RPM'] = cov.apply(lambda x: x['numreads']/miReads, axis = 1)
        cov['RPKM'] = cov.apply(lambda x: x['RPM']/(x['endpos']/1000), axis = 1)
        cov = cov.sort_values('RPKM', ascending=False)
        cov['organism'] = cov['genus'].map(str) + ' ' + cov['species'].map(str)
        cov[['organism', 'RPKM']].to_csv(output.rpkm, sep="\t", index=False)

