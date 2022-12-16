shell.executable("bash")


rule map_reads_panel:
    input:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
    output:
        aln="{sample}/panel_mapping/{sample}_aln.sam",
    params:
        refs=config["panel"],
        diag=config['bwa_d'],
        reseed=config['bwa_r'],
        A=config['bwa_A'],
        B=config['bwa_B'],
        O=config['bwa_O'],
        E=config['bwa_E'],
        L=config['bwa_L'],
        U=config['bwa_U'],
    threads:
        config['threads_sample']
    message:
        "Mapping {wildcards.sample} to panel sequences"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/mapping_panel.log"
    shell:
        """
        exec 2> {log}
        bwa mem -t {threads} \
            -d {params.diag} -r {params.reseed} \
            -A {params.A} -B {params.B} -O {params.O} \
            -E {params.E} -L {params.L} -U {params.U} \
            {params.refs} {input.r1} {input.r2} > {output.aln}
        """


rule get_mapped_panel:
    input:
        sam="{sample}/panel_mapping/{sample}_aln.sam",
    output:
        mapped="{sample}/panel_mapping/{sample}_mapped.sam",
        munmap="{sample}/panel_mapping/{sample}_munmapped.sam",
        unmap="{sample}/panel_mapping/{sample}_unmapped.sam",
        proper="{sample}/panel_mapping/{sample}_propermapped.sam",
    message:
        "Extracting panel mapped reads for {wildcards.sample}"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{sample}/extract_panel_mapped.log"
    shell:
        """
        exec 2> {log}
        # -collect all reads where at least one mate maps
        samtools view -h -F 4 -f 8 -o {output.munmap} {input.sam}
        samtools view -h -F 8 -f 4 -o {output.unmap} {input.sam}
        samtools view -h -f 2 -o {output.proper} {input.sam}
        
        samtools merge {output.mapped} {output.munmap} {output.unmap} {output.proper}
        """



rule mapping_stats_panel:
    input:
        mapped="{sample}/panel_mapping/{sample}_mapped.sam",
    output:
        flagstats="{sample}/panel_mapping/{sample}_flagstats.tsv",
        mapping_stats="{sample}/panel_mapping/{sample}_mapstats.tsv",
        index_stats="{sample}/panel_mapping/{sample}_index_stats.tsv",
    message:
        "Retrieving panel mapping stats for {wildcards.sample}"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/{sample}/panel_mapping_stats.log"
    shell:
        """
        exec 2> {log}
        samtools flagstat -O tsv {input.mapped} > {output.flagstats}
        samtools stats {input.mapped} > {output.mapping_stats}
        samtools view -f 2 -q 60 {input.mapped} | cut -f3 | sort | uniq -c | sort -k2 > {output.index_stats}
        """
