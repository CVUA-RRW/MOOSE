shell.executable("bash")


rule target_enriched_assembly:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        gfa="{sample}/assembly/{sample}_assembly_graph.gfa",
        contigs="{sample}/assembly/{sample}_contigs.fasta",
        all_contigs="{sample}/assembly/contigs_allvariants.fasta",
    message:
        "[{wildcards.sample}] read assembly with SAUTE"
    params:
        panel=config['panel'],
    threads:
        config['threads_sample']
    conda:
        "../envs/saute.yaml"
    log:
        "logs/{sample}/assembly_saute.log"
    shell:
        """
        saute --cores {threads} \
            --gfa {output.gfa} --reads {input.r1},{input.r2} \
            --targets {params.panel} --all_variants {output.all_contigs} --selected_variants {output.contigs} \
            --vector_percent 1 --extend_ends &> {log}
        """


rule index_contigs:
    input:
        contigs="{sample}/assembly/{sample}_contigs.fasta",
    output:
        multiext("{sample}/assembly/{sample}_contigs.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    message:
        "[{wildcards.sample}] indexing contigs"
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
        contigs="{sample}/assembly/{sample}_contigs.fasta",
        index=lambda wildcards: expand("{sample}/assembly/{sample}_contigs.fasta{ext}", sample=wildcards.sample, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
    output:
        aln="{sample}/assembly/{sample}_aln.bam",
    threads:
        config['threads_sample']
    message:
        "[{wildcards.sample}] mapping reads to contigs with BWA"
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
        depth="{sample}/assembly/{sample}_coverage_depth.tsv",
    message:
        "[{wildcards.sample}] calculating contigs coverage"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/host_coverage.log"
    shell:
        """
        exec 2> {log}
        # Coverage and calculate RPK and add header
        samtools coverage {input.aln} > {output.cov}
        samtools depth -aa -H {input.aln} > {output.depth}
        """


rule contigs_rpkm:
    input:
        mapping_stats="{sample}/host_mapping/{sample}_mapstats.tsv",
        cov="{sample}/assembly/{sample}_coverage.tsv",
    output:
        rpkm="{sample}/assembly/{sample}_rpkm.tsv",
    message:
        "[{wildcards.sample}] Calculation RPKM per contigs"
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


# rule connect_assembly:
    # input:
        # r1="{sample}/trimmed/{sample}_R1.fastq",
        # r2="{sample}/trimmed/{sample}_R2.fastq",
        # contigs="{sample}/assembly/contigs.fasta",
    # output:
        # gfa="{sample}/assembly/{sample}_assembly_graph_connected.gfa",
        # csv="{sample}/assembly/{sample}_gfa_connection.csv",
    # conda:
        # "../envs/saute.yaml"
    # message:
        # "Attempting to connect contigs for {wildcards.sample}"
    # threads:
        # config['threads_sample']
    # log:
        # "logs/{sample}/gfa_connector.log"
    # shell:
        # """
        # gfa_connector --cores {threads} \
            # --reads {input.r1},{input.r2} --contigs {input.contigs}\
            # --gfa {output.gfa} --csv {output.csv} \
            # --vector_percent 1 &> {log}
        # """

# rule contigs_qc:
    # input:
        # contigs="{sample}/assembly/{sample}_contigs.fasta",
    # output:
        # quast_dir=directory("{sample}/assembly/quast/"),
        # report="{sample}/assembly/quast/report.tsv",
    # message:
        # "Quality control of {wildcards.sample} assembly"
    # conda:
        # "../envs/quast.yaml"
    # log:
        # "logs/{sample}/assembly_qc.log"
    # shell:
        # """
        # quast.py --output-dir {output.quast_dir} {input.contigs} &> {log}
        # """

