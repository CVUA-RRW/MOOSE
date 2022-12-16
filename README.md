# MOOSE
Modified Organism Screening by Enrichment

## Prototyping plan

Loosely based on the Debode paper:

1. Quality trimming (fastp)
2. Align reads to reference genomes of big 6 (bwa)
3. Get unaligned reads
2b. Align reads to panel
3b. Get all reads where at least one mate maps
4a. Merge all reads
4b. Make contigs (Spades)
5. Find elements in contigs (BLAST)
6. Count reads in contigs (per found elements? Per contig?)
7. Compare host reads / transgene reads

## TODO

- recheck all the flags settings on samtools ops
- perform samtools ops on bam
- change hacky seq counts to idxstats
- convert counts to rpkm

## Ideas for future dev

- Compare bwa to simple kraken?
- If bwa upgrade to bwa2
- Test other assemblers (mira/masurca)
- Assemble multiple contigs sources (CISA?)
- Contig Blast viewer
- Automated event detection (API to BCH/Euginius)?

## Notes

### BWA reference indexing:

1-Get assmblies from e.g. Genbank (wget or use ncbi browser)

2- concatenate asemblies:

```
cd path/to/assemblies
cat *.fna > merged_refs.fa
```

3- Index genomes

```
bwa index merged_refs.fa
```

Then use merged_refs.fa as a reference

### BLAST database for panel elements:

```
makeblastdb -in panel.fa -dbtype nucl
```