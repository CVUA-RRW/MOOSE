# MOOSE
Modified Organism Screening by Enrichment

## Prototyping plan

Loosely based on the Debode paper:

1. Quality trimming (fastp)
2. Align reads to reference genomes of big 6 (bwa?)
3. Get unaligned reads
4. Make contigs (Spades?)
5. Find elements in contigs (BLAST)
6. Count reads in contigs (per found elements? Per contig?)
7. Compare host reads / transgene reads

## Ideas for future dev

- Test other assemblers (mira, 
- Assemble multiple contigs sources with CISA
- Contig Blast viewer
- Automated event detection (API to BCH/Euginius?)
