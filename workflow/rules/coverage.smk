shell.executable("bash")


rule index_contigs:
    input:
        contigs="{sample}/assembly/contigs.fasta",
    output:
        multiext("{sample}/assembly/contigs.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    message:
        "Indexing contigs for {wildcards.sample}"
    log:
        "logs/{sample}/index_contigs.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input.contigs} &> {log}
        """


rule remap_to_contigs:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
        contigs="{sample}/assembly/contigs.fasta",
        index=lambda wildcards: expand("{sample}/assembly/contigs.fasta{ext}", sample=wildcards.sample, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
    output:
        aln="{sample}/assembly/{sample}_aln.bam",
    threads:
        config['threads_sample']
    message:
        "Mapping {wildcards.sample} to host genomes"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/mapping_to_contigs.log"
    shell:
        """
        exec 2> {log}
        bwa mem -t {threads} \
            {input.contigs} {input.r1} {input.r2} \
            | samtools sort -@{threads} -o {output.aln} 
        """


rule coverage_contigs:
    input:
        aln="{sample}/assembly/{sample}_aln.bam",
    output:
        cov="{sample}/assembly/{sample}_coverage.tsv",
    message:
        "Calculating contigs coverage for {wildcards.sample}"
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


rule contigs_rpkm:
    input:
        mapping_stats="{sample}/host_mapping/{sample}_mapstats.tsv",
        cov="{sample}/assembly/{sample}_coverage.tsv",
    output:
        rpkm="{sample}/assembly/{sample}_rpkm.tsv",
    message:
        "Calculation RPKM per contigs for {wildcards.sample}"
    log:
        "logs/{sample}/rpkm_contigs.log"
    run:
        import sys
        sys.stderr = open(log[0], "w")
        
        import pandas as pd
        
        # get total reads
        with open(input.mapping_stats, 'r') as fi:
            for line in fi:
                l = line.split("\t")
                if l[0] == "SN" and l[1] == "raw total sequences:":
                    miReads = int(l[2])/1000000
        
        # Calculate RPKM for each contigs
        cov = pd.read_csv(input.cov, sep="\t")
        cov['RPM'] = cov.apply(lambda x: x['numreads']/miReads, axis = 1)
        cov['RPKM'] = cov.apply(lambda x: x['RPM']/(x['endpos']/1000), axis = 1)
        cov = cov.sort_values('RPKM', ascending=False)
        cov[['#rname', 'RPKM']].to_csv(output.rpkm, sep="\t", index=False)