shell.executable("bash")


rule blast_contigs:
    input:
        contigs="{sample}/assembly/{sample}_contigs.fasta",
    output:
        blast="{sample}/blast/{sample}_blast.asn",
        tabular="{sample}/blast/{sample}_blast.tsv",
        sam="{sample}/blast/{sample}_blast.sam",
        pairs="{sample}/blast/{sample}_blast.txt",
        
    params:
        panel=config['panel'],
    threads:
        config['threads_sample']
    message:
        "[{wildcards.sample}] aligning contigs to panel with BLAST"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}/blast.log"
    shell:
        """
        exec 2> {log}
        
        blastn -query {input.contigs} -out {output.blast} -db {params.panel} \
            -task megablast \
            -outfmt 11 -num_threads {threads}
        
        blast_formatter -archive {output.blast} -out {output.tabular} \
            -outfmt "6 qseqid sseqid slen qstart qend sstart send bitscore length sstrand pident mismatch "
        echo "\n" >> {output.tabular}
        sed -i '1 i\query_id\tsubject_id\tsubjectlength\tquery_start\tquery_end\tsubject_start\tsubject_end\tbitscore\tlength\tstrand\tidentity\tmismatch' {output.tabular}
                
        blast_formatter -archive {output.blast} -out {output.sam} -outfmt 17
        
        blast_formatter -archive {output.blast} -out {output.pairs} -outfmt 0
        """


rule find_events:
    input:
        contigs="{sample}/assembly/{sample}_contigs.fasta",
    output:
        blast="{sample}/blast/{sample}_events.asn",
        tabular="{sample}/blast/{sample}_events.tsv",
        sam="{sample}/blast/{sample}_events.sam",
        pairs="{sample}/blast/{sample}_events.txt",
    params:
        events=config['event_detection']
    message:
        "[{wildcards.sample}] detecting event-specific sequences with BLAST"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}/event_blast.log"
    shell:
        """
        exec 2> {log}
        
        blastn -query {params.events} -out {output.blast} -subject {input.contigs} \
            -task blastn -outfmt 11 -subject_besthit -perc_identity 100
        
        blast_formatter -archive {output.blast} -out {output.tabular} \
            -outfmt "6 qseqid sseqid slen qlen qstart qend sstart send bitscore sstrand length pident" 
        echo "\n" >> {output.tabular}
        sed -i '1 i\query_id\tsubject_id\tsubjectlength\tquery_length\tquery_start\tquery_end\tsubject_start\tsubject_end\tbitscore\tstrand\tlength\tidentity' {output.tabular}
        
        blast_formatter -archive {output.blast} -out {output.sam} -outfmt 17
        
        blast_formatter -archive {output.blast} -out {output.pairs} -outfmt 0
        """

