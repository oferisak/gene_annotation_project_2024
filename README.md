A project with scripts that produces premade annotation tables for all the genes to be used in downstream analysis
scripts include:

## Usage: 
When adding a new target / changing an annotation file like clinvar, take these steps to prepare it for the app
Step 1. generate the target tables:
- ClinVar coverage = produce_clinvar_vs_target_tables - just update the target list and run the analysis (it will generate for all the targets). output file = {target}.{clinvar_file}.clinvar_coverage.missing_variants.tsv
- target coverage = produce_genes_vs_target_tables - output file = {target}.refseq.incomplete_coverage.tsv
- GIAB stratifications = produce_giab_stratification_vs_gene_tables - you first need to select the GIAB strats you want and put them in the data/giab folder (hg19/grch37), then run the script - output file = {target}.giab_stratifications.non_zero.tsv

Step 2. after the tables are ready, copy them to the app/data/target_data folder


## GIAB stratifications

### Segmental duplications

#### File descriptions:

- GRCh3X_chainSelf.bed.gz
Describes non-trivial alignments of the genome reference to itself (excluding ALT loci). Non-trivial self-chains are regions where another part of the genome aligns to it because the sequences are similar (e.g., due to genomic duplication events). Further information on UCSC tracks can be found at https://genome.ucsc.edu/cgi-bin/hgTables.

- GRCh3X_segdups.bed.gz
The UCSC genomicSuperDups file was processed to remove all but chromosomes 1-22, X and Y and overlapping regions merged.

### Low complexity regions

#### File Descriptions:
All beds have 5bp slop added on each side to capture variants at the edge of the repeats (sometimes insertions were not captured properly before in stratifications).

- GRCh3X_SimpleRepeat*_slop5.bed.gz
perfect repeats of different unit sizes (i.e., homopolymers, and dinucleotide, trinucleotide, and quadnucleotide STRs) and different total repeat lengths (i.e., <=50bp, 51-200bp, or >200bp)

- GRCh3X_SimpleRepeat_imperfecthomopolgt*_slop5.bed.gz
perfect homopolymers >*p as well as imperfect homopolymers where a single base was repeated >10bp except for a 1bp interruption by a different base

- GRCh3X_AllTandemRepeats_*_slop5.bed.gz
union of SimpleRepeat dinucleotide, trinucleotide, and quadnucleotide STRs as well as UCSC Genome Brower RepeatMasker_LowComplexity, RepeatMasker_SimpleRepeats, RepeatMasker_Satellite, and TRF_SimpleRepeat.

- GRCh3X_AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz
union of all perfect homopolymers >6bp and imperfect homopolymers >10bp

- GRCh3X_AllTandemRepeatsandHomopolymers_slop5.bed.gz
union of AllTandemRepeats_* with AllHomopolymers_gt6bp_imperfectgt10bp_slop5.bed.gz

- GRCh3X_satellites_slop5.bed.gz
Centromeric and Pericentromeric Satellite Annotations

- GRCh3X_notin*_slop5.bed.gz
are non-overlapping complements of the stratification regions (i.e., genome after excluding the regions).

- GRCh3X_allTandemRepeats.bed.gz
union of all tandem repeats

