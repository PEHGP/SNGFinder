# SNGFinder
SNGFinder is a tool for identifying new genes based on the syntenic method. It requires a focal species (the species for which the new gene needs to be identified), species of the same lineage or genus as the focal species, and outgroup species. Syntenic regions are obtained based on whole genome alignment, and clusters of genes orthologous to the focal species genes are obtained within the same syntenic region. If a focal species gene has no orthologous genes in any outgroup in a orthologous cluster, then the gene is a new gene. In addition, SNGFinder can classify new genes into copy number variation genes, duplicate genes, and orphan genes based on the origin mechanism.
# Installation

# How to use
## 1. Repeat masker
If the genome has been repeat soft-masked, skip this step.  
The genomes downloaded from NCBI or ensemble usually have been repeat soft-masked.  
`RepeatMasker` is recommended, but `WindowMasker` is used here.  
```
windowmasker -mk_counts -in genome.fasta -out genome.counts
windowmasker -ustat genome.counts -in genome.fasta -out genome_masked.fasta -dust true -outfmt fasta
```
`genome_masked.fasta` is the repeat soft-masking fasta genome.
## 2. Whole genome alignment
### Alignment
```
cactus --binariesMode local --defaultMemory 100G --defaultDisk 30G job_store SpAll.txt SpAll.hal
```
`job_store` is used to store intermediate files.
`SpAll.txt` is a text file containing the locations of the input sequences as well as their phylogenetic tree. For a detailed explanation of this file, please refer to the [cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#interface). The file is formatted as follows:
```
(Sp1:0.23781,((Sp2:0.335387,Sp3:0.399047)N2:0.0755071,((Sp4:0.0245695,Sp5:0.0279553)N4:0.0525508,(Sp6:0.0735841,(Sp7:0.072269,(Sp8:0.0329604,Sp9:0.0467526)N7:0.0277942)N6:0.028579)N5:0.0253669)N3:0.281769)N1:0.23781)N0;
Sp1 Sp1_genome.fasta.gz
Sp2 Sp2_genome.fasta.gz
Sp3 Sp3_genome.fasta.gz
Sp4 Sp4_genome.fasta.gz
Sp5 Sp5_genome.fasta.gz
Sp6 Sp6_genome.fasta.gz
Sp7 Sp7_genome.fasta.gz
Sp8 Sp8_genome.fasta.gz
Sp9 Sp9_genome.fasta.gz
```
`SpAll.hal` is cactus output
[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) may be slow, so you can use a [cluster](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster) or [GPU](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#gpu-acceleration) to speed it up.
### MAF export
```
hal2maf --maxBlockLen 100000 --noAncestors --noDupes --onlyOrthologs --unique --refGenome Sp1 SpAll.hal stdout|taffy view|taffy add-gap-bases -a SpAll.hal|taffy norm -k -d |mafFilter -m - -N 0.95 |mafDuplicateFilter -m - -k |bgzip >SpAll_ref_Sp1.maf.gz
```
`SpAll.hal` is the input.
`--refGenome` a focal species for new gene identification.
`SpAll_ref_Sp1.maf.gz` is the maf format output.
For detailed instructions on MAF export, please refer to [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#maf-export)
## 3. Identification of new genes
# How to cite