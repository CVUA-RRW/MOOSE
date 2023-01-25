# MOOSE
Modified Organism Screening by Enrichment

## Prototyping plan

Loosely based on the Debode paper:

1. Quality trimming (fastp)
2. Align reads to reference genomes of big 6 (bwa)
3. Extract RPKM values for each host
4. reference alignement to elements (SAUTE) + polishing (GFA connector ?)
5. Find elements in contigs (BLAST)
6. Count reads in contigs (per found elements? Per contig?)
7. Try to detect event specific sequences (JRC-EURL Methods) in contigs

## TODO

- Create reports
- chop contigs and reblast part with no matches? 
	set a min length to chop (100bp?)
	use full BLAST db (local)
	hit selection? (best hit?)

## Ideas for future dev

- Compare bwa to simple kraken?
- If bwa upgrade to bwa2
- Automated event detection (API to BCH/Euginius)?

## Notes

### BWA reference indexing:

1-Get assemblies from e.g. Genbank (wget or use ncbi browser)

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

### Panel sequences
Multifasta format, unwrapped sequences!