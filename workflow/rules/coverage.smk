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
