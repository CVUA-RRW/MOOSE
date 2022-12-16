shell.executable("bash")


rule blast_scaffolds:
    input:
        scaffold="{sample}/assembly/spades/scaffolds.fasta",
    output:
        blast="{sample}/blast/{sample}_blast.tsv",
    params:
        panel=config['panel'],
        blast_id=config['blast_id'],
        blast_cov=config['blast_cov'],
    threads:
        config['threads_sample']
    message:
        "Aligning scaffold to panel for {wildcards.sample}"
    conda:
        "../envs/blast.yaml"
    log:
        "logs/{sample}/blast.log"
    shell:
        """
        exec 2> {log}
        
        blastn -query {input.scaffold} -out {output.blast} -db {params.panel} \
            -task megablast \
            -outfmt "6 qseqid sseqid slen qstart qend sstart send bitscore length pident mismatch " \
            -num_threads {threads}
        
        sed -i '1 i\query_id\tsubject_id\tsubjectlength\tquery_start\tquery_end\tsubject_start\tsubject_end\tbitscore\tlength\tidentity\tmismatch' {output.blast}
            
        """