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
        assembly_kmer_min_count=config['assembly_kmer_min_count'],
        assembly_noise_to_signal=config['assembly_noise_to_signal'],
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
            --min_count {params.assembly_kmer_min_count} --fraction {params.assembly_noise_to_signal} --protect_reference_ends \
            --vector_percent 1 --extend_ends &> {log}
        """


rule collapse_identical_contigs:
    input:
        contigs="{sample}/assembly/{sample}_contigs.fasta",
    output:
        contigs="{sample}/assembly/{sample}_collapsed_contigs.fasta",
        table="{sample}/assembly/{sample}_clustering.tsv",
    message:
        "[{wildcards.sample}] collapsing highly similar contigs with VSEARCH"
    params:
        id=config["assembly_identity_collapse"],
    conda:
        "../envs/vsearch.yaml"
    log:
        "logs/{sample}/collapse_vsearch.log"
    shell:
        """
        exec 2> {log}
        vsearch --cluster_fast {input.contigs} --consout {output.contigs} --id {params.id} \
            --strand both -uc {output.table} --fasta_width 0 --xsize
        #Fix Headers
        sed -i -E  -e 's/^>centroid=/>/' -e 's/;seqs=[0-9]+$//' {output.contigs}
        """


rule index_contigs:
    input:
        contigs="{sample}/assembly/{sample}_collapsed_contigs.fasta",
    output:
        multiext("{sample}/assembly/{sample}_collapsed_contigs.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
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
        contigs="{sample}/assembly/{sample}_collapsed_contigs.fasta",
        index=lambda wildcards: expand("{sample}/assembly/{sample}_collapsed_contigs.fasta{ext}", sample=wildcards.sample, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
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

