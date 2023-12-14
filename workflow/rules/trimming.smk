shell.executable("bash")


# Rules quality trimming ------------------------------------------------------


rule run_fastp:
    input:
        r1=lambda wildcards: get_fastq(wildcards, "fq1"),
        r2=lambda wildcards: get_fastq(wildcards, "fq2"),
    output:
        r1=temp("{sample}/trimmed/{sample}_R1.fastq"),
        r2=temp("{sample}/trimmed/{sample}_R2.fastq"),
        json="{sample}/trimmed/{sample}.json",
        html="{sample}/trimmed/{sample}.html",
    params:
        length_required=config["read_length_required"],
        qualified_quality_phred=config["qualified_quality_phred"],
        window_size=config["qctrim_window_size"],
        mean_qual=config["qctrim_mean_quality"],
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
            --thread {threads} \
            --report_title 'Sample {wildcards.sample}' \
        > {log} 2>&1
        """


rule parse_fastp:
    input:
        json="{sample}/trimmed/{sample}.json",
        html="{sample}/trimmed/{sample}.html",
    output:
        tsv=temp("{sample}/reports/{sample}_fastp.tsv"),
    params:
        sample=lambda w: w.sample,
    message:
        "[{wildcards.sample}] parsing FASTP json report"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/{sample}/parse_fatsp.log",
    script:
        "../scripts/parse_fastp.py"


rule collect_trimming_stats:
    input:
        report=expand("{sample}/reports/{sample}_fastp.tsv", sample=samples.index),
    output:
        agg="reports/fastp_stats.tsv",
    message:
        "[All] aggregating fastp stats"
    log:
        "logs/all/trimming_stats.log",
    conda:
        "../envs/pandas.yaml"
    shell:
        """
        exec 2> {log}
        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """

