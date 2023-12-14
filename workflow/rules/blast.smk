shell.executable("bash")


rule blast_contigs:
    input:
        contigs="{sample}/assembly/{sample}_selected_contigs.fasta",
    output:
        blast="{sample}/blast/{sample}_blast.asn",
        tabular="{sample}/blast/{sample}_blast.tsv",
        sam="{sample}/blast/{sample}_blast.sam",
        pairs="{sample}/blast/{sample}_blast.txt",
    params:
        elements_detection=config['elements_detection'],
        blast_evalue=config['blast_evalue'],
        blast_identity=config['blast_identity'],
        blast_coverage=config['blast_coverage'],
    threads:
        config['threads_sample']
    message:
        "[{wildcards.sample}] aligning panel elements to contigs with BLAST"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}/blast.log"
    shell:
        """
        exec 2> {log}
        if [ -s {input.contigs} ]; then
            blastn -query  {params.elements_detection} -out {output.blast} -subject {input.contigs} \
                -task megablast -strand both \
                -max_hsps 150 \
                -evalue {params.blast_evalue} \
                -perc_identity {params.blast_identity} \
                -qcov_hsp_perc {params.blast_coverage} \
                -outfmt 11 -num_threads {threads}
            
            blast_formatter -archive {output.blast} -out {output.tabular} \
                -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send bitscore sstrand length pident mismatch"
            echo "\n" >> {output.tabular}
            
            blast_formatter -archive {output.blast} -out {output.sam} -outfmt 17
            
            blast_formatter -archive {output.blast} -out {output.pairs} -outfmt 0
        else
            touch {output.blast}
            touch {output.tabular}
            touch {output.sam}
            touch {output.pairs}
        fi
        sed -i '1 i\query_id\tsubject_id\tquery_length\tsubjectlength\tquery_start\tquery_end\tsubject_start\tsubject_end\tbitscore\tstrand\tlength\tidentity\tmismatch' {output.tabular}
        """


# rule find_events:
    # input:
        # contigs="{sample}/assembly/{sample}_contigs.fasta",
    # output:
        # blast="{sample}/blast/{sample}_events.asn",
        # tabular="{sample}/blast/{sample}_events.tsv",
        # sam="{sample}/blast/{sample}_events.sam",
        # pairs="{sample}/blast/{sample}_events.txt",
    # params:
        # events=config['event_detection']
    # message:
        # "[{wildcards.sample}] detecting event-specific sequences with BLAST"
    # conda:
        # "../envs/blast.yaml"
    # log:
        # "logs/{sample}/event_blast.log"
    # shell:
        # """
        # exec 2> {log}
        
        # blastn -query {params.events} -out {output.blast} -subject {input.contigs} \
            # -task blastn -outfmt 11 -subject_besthit -perc_identity 100
        
        # blast_formatter -archive {output.blast} -out {output.tabular} \
            # -outfmt "6 qseqid sseqid qlen slen qstart qend sstart send bitscore sstrand length pident mismatch"
        # echo "\n" >> {output.tabular}
        # sed -i '1 i\query_id\tsubject_id\tquery_length\tsubjectlength\tquery_start\tquery_end\tsubject_start\tsubject_end\tbitscore\tstrand\tlength\tidentity\tmismatch' {output.tabular}
        
        # blast_formatter -archive {output.blast} -out {output.sam} -outfmt 17
        
        # blast_formatter -archive {output.blast} -out {output.pairs} -outfmt 0
        # """


rule blast_coverage:
    input:
        blast="{sample}/blast/{sample}_blast.tsv",
        bam="{sample}/assembly/{sample}_aln.bam",
    output:
        bed="{sample}/blast/{sample}_blast.bed",
        cov="{sample}/reports/{sample}_blastcoverage.tsv",
    params:
        sample=lambda w: w.sample,
        mean_mapQ=config["mean_mapQ"],
    message:
        "[{wildcards.sample}] recovering coverage of blast results"
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{sample}/blast_coverage.log"
    shell:
        """
        # make blast as bed file https://www.biostars.org/p/317112/
        # BED output fields are: contig start_pos end_pos blast_query
        # Note in BED start must be lower than end
        cat {input.blast} \
            | tail -n +2 \
            | sed '/^$/d' \
            | awk -v OFS='\t' '{{print $2, ($7<$8)?$7-1:$8, ($7>$8)?$7-1:$8, $1}}' \
            > {output.bed}
        # coverage with bedtools - or empty file
        if [ -s {output.bed} ]; then
            samtools index {input.bam}
            samtools bedcov {output.bed} {input.bam} \
                -g SECONDARY \
                | awk -v OFS='\t' '{{print "{params.sample}", $1, $2, $3, $4, $5/($3-$2)}}' \
                > {output.cov}
        else
            touch {output.cov}
        fi
        # add header
        sed -i '1isample\trname\tstart\tend\tquery\tmeandepth' {output.cov}
        """


rule report_blast:
    input:
        report=expand("{sample}/reports/{sample}_blastcoverage.tsv", sample=samples.index),
    output:
        agg="reports/blast.tsv"
    message:
        "[all] Aggregating blast stats"
    log:
        "logs/report_blast.log"
    shell:
        """
        exec 2> {log}

        cat {input.report[0]} | head -n 1 > {output.agg}
        for i in {input.report}; do 
            cat ${{i}} | tail -n +2 >> {output.agg}
        done
        """