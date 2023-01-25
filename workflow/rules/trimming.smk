shell.executable("bash")


# Rules quality trimming ------------------------------------------------------


rule run_fastp:
    input:
        r1=lambda wildcards: get_fastq(wildcards, "fq1"),
        r2=lambda wildcards: get_fastq(wildcards, "fq2"),
    output:
        r1="{sample}/trimmed/{sample}_R1.fastq",
        r2="{sample}/trimmed/{sample}_R2.fastq",
        json="{sample}/trimmed/{sample}.json",
        html="{sample}/trimmed/{sample}.html",
    params:
        length_required=config["read_length_required"],
        qualified_quality_phred=config["qualified_quality_phred"],
        window_size=config["qctrim_window_size"],
        mean_qual=config["qctrim_mean_quality"],
        umi="--umi" if config['umi'] else "",
        umi_loc=f"--umi_loc {config['umi_loc']}" if config['umi'] else "",
        umi_len=f"--umi_len {config['umi_len']}"if config['umi'] else "",
    threads: config["threads_sample"]
    message:
        "[{wildcards.sample}] quality trimming with FASTP"
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/{sample}/fastp.log",
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --cut_by_quality3 \
            --cut_window_size {params.window_size} \
            --cut_mean_quality {params.mean_qual} \
            --detect_adapter_for_pe \
            --thread {threads} \
            --report_title 'Sample {wildcards.sample}' \
            {params.umi} {params.umi_loc} {params.umi_len} \
        > {log} 2>&1
        """