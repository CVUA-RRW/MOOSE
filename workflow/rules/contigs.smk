shell.executable("bash")


rule get_fastq:
    input:
        host_unmapped="{sample}/host_mapping/{sample}_unmapped.sam",
        panel_mapped="{sample}/panel_mapping/{sample}_mapped.sam"
    output:
        merged="{sample}/assembly/{sample}_merged_reads.sam",
        paired1="{sample}/assembly/{sample}_pair1.fastq",
        paired2="{sample}/assembly/{sample}_pair2.fastq",
        unpaired="{sample}/assembly/{sample}_unpaired.fastq",
        mergedstats="{sample}/assembly/{sample}_merged_reads_stats.tsv",
    message:
        "Collecting reads and creating fastq for {wildcards.sample}"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{sample}/create_fastq.log"
    shell:
        """
        exec 2> {log}
        
        samtools merge -o {output.merged} {input.host_unmapped} {input.panel_mapped}
        samtools stats {output.merged} > {output.mergedstats}
        samtools fastq -1 {output.paired1} -2 {output.paired2} -0 /dev/null -s {output.unpaired} -n {output.merged}
        """


rule create_contigs:
    input:
        paired1="{sample}/assembly/{sample}_pair1.fastq",
        paired2="{sample}/assembly/{sample}_pair2.fastq",
        unpaired="{sample}/assembly/{sample}_unpaired.fastq",
    output:
        spades_dir=directory("{sample}/assembly/denovo/"),
        contigs="{sample}/assembly/denovo/contigs.fasta",
        scaffold="{sample}/assembly/denovo/scaffolds.fasta",
    message:
        "Assembling filtered reads for {wildcards.sample} with skesa"
    params:
        memory=config['memory'],
        temp=config['temp'],
    conda:
        "../envs/spades.yaml"
    threads:
        config['threads_sample']
    log:
        "logs/{sample}/assembly_denovo.log"
    shell:
        """
        spades.py -o {output.spades_dir} \
            -k 21,33,55,77 \
            -1 {input.paired1} -2 {input.paired2} -s {input.unpaired} \
            -t {threads} -m {params.memory} --tmp-dir {params.temp} &> {log}
        """


rule contigs_qc:
    input:
        contigs="{sample}/assembly/denovo/contigs.fasta",
    output:
        quast_dir=directory("{sample}/assembly/denovo/quast/"),
        report="{sample}/assembly/denovo/quast/report.tsv",
    message:
        "Quality control of {wildcards.sample} assembly"
    conda:
        "../envs/quast.yaml"
    log:
        "logs/{sample}/assembly_denovo_qc.log"
    shell:
        """
        quast.py --output-dir {output.quast_dir} {input.contigs} &> {log}
        """


rule create_contigs_saute:
    input:
        paired1="{sample}/assembly/{sample}_pair1.fastq",
        paired2="{sample}/assembly/{sample}_pair2.fastq",
    output:
        gfa="{sample}/assembly/saute/assembly_graph.gfa",
        contigs="{sample}/assembly/saute/contigs.fasta",
    message:
        "Assembling filtered reads for {wildcards.sample} with saute"
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
            --gfa {output.gfa} --reads {input.paired1},{input.paired2} \
            --targets {params.panel} --all_variants {output.contigs} \
            --vector_percent 1 --extend_ends &> {log}
        """


rule create_contigs_saute_all_reads:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        gfa="{sample}/assembly/saute_all/assembly_graph.gfa",
        contigs="{sample}/assembly/saute_all/contigs.fasta",
    message:
        "Assembling all reads for {wildcards.sample} with saute"
    params:
        panel=config['panel'],
    threads:
        config['threads_sample']
    conda:
        "../envs/saute.yaml"
    log:
        "logs/{sample}/assembly_saute_all.log"
    shell:
        """
        saute --cores {threads} \
            --gfa {output.gfa} --reads {input.r1},{input.r2} \
            --targets {params.panel} --all_variants {output.contigs} \
            --vector_percent 1 --extend_ends &> {log}
        """
