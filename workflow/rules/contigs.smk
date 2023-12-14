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
        assembly_template=config['assembly_template'],
        assembly_kmer_min_count=config['assembly_kmer_min_count'],
        assembly_noise_to_signal=config['assembly_noise_to_signal'],
        assembly_sec_kmer_threshold=config['assembly_sec_kmer_threshold'],
        assembly_target_coverage=config['assembly_target_coverage'],
    threads:
        config['threads_sample']
    conda:
        "../envs/saute.yaml"
    log:
        "logs/{sample}/assembly_saute.log"
    shell:
        """
        saute --cores {threads} \
            --gfa {output.gfa} \
            --reads {input.r1},{input.r2} \
            --targets {params.assembly_template} \
            --all_variants {output.all_contigs} \
            --selected_variants {output.contigs} \
            --min_count {params.assembly_kmer_min_count} \
            --fraction {params.assembly_noise_to_signal} \
            --secondary_kmer_threshold {params.assembly_sec_kmer_threshold} \
            --target_coverage {params.assembly_target_coverage} \
            --protect_reference_ends \
            --use_ambiguous_na \
            --extend_ends &> {log}
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
        # Fix Headers
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
        aln=temp("{sample}/assembly/{sample}_aln.bam"),
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
        # Warning, with -a will output all mapping locations
        if [ -s {input.contigs} ]; then
            bwa mem -t {threads} -a \
                {input.contigs} {input.r1} {input.r2} \
                | samtools sort -@{threads} -o {output.aln}
        else
            touch {output.aln}
        fi
        """


rule coverage_contigs:
    input:
        aln="{sample}/assembly/{sample}_aln.bam",
    output:
        cov="{sample}/assembly/{sample}_coverage.tsv",
        depth="{sample}/assembly/{sample}_coverage_depth.tsv",
    params:
        sample=lambda w: w.sample,
        mean_mapQ=config["mean_mapQ"],
    message:
        "[{wildcards.sample}] calculating contigs coverage"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/host_coverage.log"
    shell:
        """
        exec 2> {log}
        if [ -s {input.aln} ]; then
            # Coverage and calculate and add header
            # coverage specify ecl falgs so that SECONDARY are included
            samtools coverage {input.aln} \
                --min-MQ {params.mean_mapQ} \
                --excl-flags UNMAP,QCFAIL,DUP \
                | awk -v OFS='\t' 'NR==1 {{ print "sample", $0}} NR>1 {{ print "{params.sample}", $0}}' \
                > {output.cov}
            # with depth can add secondary aln 
            samtools depth -a \
                --min-MQ {params.mean_mapQ} \
                -g SECONDARY \
                {input.aln} \
                | awk -v OFS='\t' '{{ print "{params.sample}", $0}}' \
                > {output.depth}
        else
            echo "sample\t#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq" \
                > {output.cov}
            touch {output.depth}
        fi
        sed -i '1isample\trname\tpos\tnreads' {output.depth}
        """


rule filter_contigs:
    input:
        contigs="{sample}/assembly/{sample}_collapsed_contigs.fasta",
        cov="{sample}/assembly/{sample}_coverage.tsv",
    output:
        filter="{sample}/reports/{sample}_coverage_filtered.tsv",
        fasta="{sample}/assembly/{sample}_selected_contigs.fasta",
        seqids="{sample}/assembly/{sample}_selected_seqids.txt",
    params:
        minlength=config["contig_min_length"],
        mindepth=config["contig_min_depth"],
    message:
        "[{wildcards.sample}] filtering contigs"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/filter_contigs.log"
    shell:
        """
        exec 2> {log}
        cat {input.cov} | head -n 1 > {output.filter}
        tail -n +2 {input.cov} |
            awk '$4>{params.minlength} && $8>{params.mindepth} {{ print $0 }}' \
            >> {output.filter}
        # Now filter fasta
        tail -n +2 {output.filter} \
            | cut -d"\t" -f2 \
            > {output.seqids}
        seqkit grep --by-name -f {output.seqids} {input.contigs} > {output.fasta}
        """


rule flagstats:
    input:
        aln="{sample}/assembly/{sample}_aln.bam",
        selected="{sample}/reports/{sample}_coverage_filtered.tsv",
    output:
        aln="{sample}/assembly/{sample}_aln_filtered.bam",
        bed="{sample}/assembly/{sample}_filtered.bed",
        flags_before="{sample}/assembly/{sample}_flagstats_unfiltered.tsv",
        flags_after="{sample}/assembly/{sample}_flagstats_filtered.tsv",
    message:
        "[{wildcards.sample}] Getting alignement flag statitics"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/flagstats.log"
    shell:
        """
        exec 2> {log}
        # make bed
        tail -n +2 {input.selected} |
            cut -d"\t" -f 2-4 \
            > {output.bed}
        # filter bam
        samtools view --bam \
            --target-file {output.bed} \
            --output {output.aln} {input.aln}
        # flagstats on both mappings
        samtools flagstat -O tsv {input.aln} > {output.flags_before}
        samtools flagstat -O tsv {output.aln} > {output.flags_after}
        """
 

rule report_contigs:
    input:
        report=expand("{sample}/reports/{sample}_coverage_filtered.tsv", sample=samples.index),
    output:
        agg="reports/contigs.tsv"
    message:
        "[all] Aggregating contig stats"
    log:
        "logs/report_contigs.log"
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """
