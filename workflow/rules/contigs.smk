shell.executable("bash")


rule target_enriched_assembly:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        gfa="{sample}/assembly/assembly_graph.gfa",
        contigs="{sample}/assembly/contigs.fasta",
        all_contigs="{sample}/assembly/contigs_allvariants.fasta",
    message:
        "Assembling all reads for {wildcards.sample} with SAUTE"
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


rule connect_assembly:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
        contigs="{sample}/assembly/contigs.fasta",
    output:
        gfa="{sample}/assembly/assembly_graph_connected.gfa",
        csv="{sample}/assembly/gfa_connection.csv",
    conda:
        "../envs/saute.yaml"
    message:
        "Attempting to connect contigs for {wildcards.sample}"
    threads:
        config['threads_sample']
    log:
        "logs/{sample}/gfa_connector.log"
    shell:
        """
        gfa_connector --cores {threads} \
            --reads {input.r1},{input.r2} --contigs {input.contigs}\
            --gfa {output.gfa} --csv {output.csv} \
            --vector_percent 1 &> {log}
        """

rule contigs_qc:
    input:
        contigs="{sample}/assembly/contigs.fasta",
    output:
        quast_dir=directory("{sample}/assembly/quast/"),
        report="{sample}/assembly/quast/report.tsv",
    message:
        "Quality control of {wildcards.sample} assembly"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/{sample}/assembly_qc.log"
    shell:
        """
        quast.py --output-dir {output.quast_dir} {input.contigs} &> {log}
        """

