# SNGFinder
SNGFinder is a tool for identifying new genes based on the syntenic method. It requires a focal species (the species for which the new gene needs to be identified), species of the same lineage or genus as the focal species, and outgroup species. Syntenic regions are obtained based on whole genome alignment, and clusters of genes orthologous to the focal species genes are obtained within the same syntenic region. If a focal species gene has no orthologous genes in any outgroup in a orthologous cluster, then the gene is a new gene. In addition, SNGFinder can classify new genes into copy number variation genes, duplicate genes, and orphan genes based on the origin mechanism.
# Installation
```
conda env create -f environment.yml -n SNGFinder
conda activate SNGFinder
pip install brotli
wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v3.1.4/cactus-bin-v3.1.4.tar.gz
tar xvzf cactus-bin-v3.1.4.tar.gz
cd cactus-bin-v3.1.4
python ./setup.py install
export PATH="$(pwd)/bin:$PATH" #important
cd ..
pip install git+https://github.com/PEHGP/SNGFinder.git
```
   
> [!NOTE]
> When the current Shell session is closed or exited, it is necessary to re-export PATH.
# How to use
Identifying new genes requires the following steps:  
1. [Repeat masker](#1-repeat-masker)
2. [Whole genome alignment](#2-whole-genome-alignment)
3. [Identification of new genes](#3-identification-of-new-genes)
## 1. Repeat masker
If the genome has been repeat soft-masked, skip this step.  
The genomes downloaded from NCBI or Ensembl usually have been repeat soft-masked.  
`RepeatMasker` is recommended, but `WindowMasker` is used here.  
```
windowmasker -mk_counts -in genome.fasta -out genome.counts
windowmasker -ustat genome.counts -in genome.fasta -out genome_masked.fasta -dust true -outfmt fasta
```
`genome_masked.fasta` is the repeat soft-masking fasta file for whole genome alignment.

## 2. Whole genome alignment
### 2.1 Alignment
```
cactus --binariesMode local --defaultMemory 100G --defaultDisk 30G job_store SpAll.txt SpAll.hal
```
`job_store` is used to store intermediate files.  
`SpAll.txt` is a text file containing the locations of the input sequences as well as their phylogenetic tree.  
  
For a detailed explanation of this file, please refer to the [cactus documentation](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#interface).  
`SpAll.hal` is cactus output  
  
**The file is formatted as follows:**  
```
(Sp1:0.23781,((Sp2:0.335387,Sp3:0.399047)N2:0.0755071,((Sp4:0.0245695,Sp5:0.0279553)N4:0.0525508,(Sp6:0.0735841,(Sp7:0.072269,(Sp8:0.0329604,Sp9:0.0467526)N7:0.0277942)N6:0.028579)N5:0.0253669)N3:0.281769)N1:0.23781)N0;
Sp1 Sp1_genome.fasta
Sp2 Sp2_genome.fasta
Sp3 Sp3_genome.fasta
Sp4 Sp4_genome.fasta
Sp5 Sp5_genome.fasta
Sp6 Sp6_genome.fasta
Sp7 Sp7_genome.fasta
Sp8 Sp8_genome.fasta
Sp9 Sp9_genome.fasta
```
> [!NOTE]
>The species name in the `SpAll.txt` can not include `.`.  
[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) may be slow, so you can use a [cluster](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster) or [GPU](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#gpu-acceleration) to speed it up.  
`Cactus` used the system temporary directory `/tmp` by default. If the prompt space is insufficient, you may need to use `--workdir` to specify the temporary directory.  
### 2.2 MAF export
```
hal2maf --maxBlockLen 100000 --noAncestors --noDupes --onlyOrthologs --unique --refGenome Sp1 SpAll.hal stdout|taffy view|taffy add-gap-bases -a SpAll.hal|taffy norm -k -d |mafFilter -m - -N 0.95 |mafDuplicateFilter -m - -k |bgzip >SpAll_ref_Sp1.maf.gz
```
`SpAll.hal` is the input.  
`--refGenome` a focal species for new gene identification.  
`SpAll_ref_Sp1.maf.gz` is the maf format output.  
  
For detailed instructions on MAF export, please refer to [here](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#maf-export).  
  
> [!NOTE]
>Errors may be reported on NFS file system. It is recommended to run on SSD.  
## 3. Identification of new genes
### Target file
**A `target` file is required, the format is as follows：**
```
Sp1 Sp1_genome.fasta.gz Sp1_protein.faa.gz  Sp1_gene.bed    Sp1_gene_to_protein.txt in
Sp2 Sp2_genome.fasta.gz Sp2_protein.faa.gz  Sp2_gene.bed    Sp2_gene_to_protein.txt in
Sp3 Sp3_genome.fasta.gz Sp3_protein.faa.gz  Sp3_gene.bed    Sp3_gene_to_protein.txt in
Sp4 Sp4_genome.fasta.gz Sp4_protein.faa.gz  Sp4_gene.bed    Sp4_gene_to_protein.txt in
Sp5 Sp5_genome.fasta.gz Sp5_protein.faa.gz  Sp5_gene.bed    Sp5_gene_to_protein.txt in
Sp6 Sp6_genome.fasta.gz Sp6_protein.faa.gz  Sp6_gene.bed    Sp6_gene_to_protein.txt in
Sp7 Sp7_genome.fasta.gz Sp7_protein.faa.gz  Sp7_gene.bed    Sp7_gene_to_protein.txt out
Sp8 Sp8_genome.fasta.gz Sp8_protein.faa.gz  Sp8_gene.bed    Sp8_gene_to_protein.txt out
Sp9 Sp9_genome.fasta.gz Sp9_protein.faa.gz  Sp9_gene.bed    Sp9_gene_to_protein.txt out
```
The `target` file has six columns, each separated by the `tab` key.  
The `first column` is species name.  
The `second column` is the fasta genome sequence.  
The `third column` is the fasta protein sequence.  
The `fourth column` is the genes bed file. The bed format only requires the first six columns, please refer to [here](https://grch37.ensembl.org/info/website/upload/bed.html) for details.  
The `fifth column` is the correspondence file between gene names and protein names. This file has two columns. The first column is the gene name and the second column is the protein name corresponding to the gene name.  
The `sixth column` is used to mark whether it is an outgroup species. If it is an outgroup species, fill in `out`. If not, fill in `in`.  
### SNGFinder
```
SNGFinder --prefix TestSp9 --maf SpAll_ref_Sp1.maf.gz --target target --refspecies Sp1 >test.log
```
`--prefix` is the result prefix name.  
`--maf` alignment file of maf format generated by cactus.  
`--target` target file.  
`--refspecies` the focus species for identifying new genes.  
  
**Detailed parameters**  
```
SNGFinder -h
usage: SNGFinder --prefix TestSp9 --maf SpAll_ref_Sp1.maf.gz --target target --refspecies Sp1

Identifying new genes based on the syntenic method.

options:
  -h, --help            show this help message and exit
  --prefix PREFIX       Result file prefix.
  --maf MAF             Whole-genome multiple sequence alignment file in maf format.
  --target TARGET       A target file.
  --blastpickle BLASTPICKLE
                        The pickle file that stores the best alignment results of blast. if this file is provided, skip the GetBlastd function.
  --proteinpickle PROTEINPICKLE
                        The pickle file that stores the protein sequence. if this file is provided, it is not necessary to read the protein sequence from the FASTA file.
  --genetotranspickle GENETOTRANSPICKLE
                        The pickle file that stores the correspondence between genes and transcripts.
  --transtogenepickle TRANSTOGENEPICKLE
                        The pickle file that stores the corresponding between transcripts and genes.
  --allgenepickle ALLGENEPICKLE
                        The pickle file that stores all gene names for each species.
  --blastout BLASTOUT   all vs. all blast results.
  --syntenicpickle SYNTENICPICKLE
                        The pickle file that stores the syntenic information between species. If this file is provided, skip the GetSyntenic function.
  --blockpickle BLOCKPICKLE
                        The pickle file that stores the syntenic gene block information between species. If this file is provided, skip the GetSyntenicGene function.
  --thpickle THPICKLE   The pickle file that stores the threshold for judging ortholog or paralog. If this file is provided, skip the GetRBH function.
  --chimericlistfile CHIMERICLISTFILE
                        The file that stores the list of chimeric genes. If this file is provided, skip the GetChimeras function.
  --refspecies REFSPECIES
                        Reference species name, which should be consistent with the species name in the target file.
  --outgroupfile OUTGROUPFILE
                        The file that stores the outgroup species, one species per line. If this file is provided, the outgroup species will be determined based on this file
                        instead of the target file.
  --dnaidentify DNAIDENTIFY
                        The threshold for judging ortholog based on DNA sequence identity. If the identity of the best blast hit in the syntenic region is greater than or equal
                        to this threshold, it will be considered as an ortholog.
  --inflation INFLATION
                        The inflation parameter for MCL clustering of paralog genes. The larger the value, the finer the clustering.
  --overlapth OVERLAPTH
                        Maximum overlap between blast hits allowed when judging chimeric genes.
  --leftpercentile LEFTPERCENTILE
                        The left percentile of the blast hits bit score distribution used to judge chimeric genes. If the blast hits bit score is less than this percentile, it
                        will be removed.
  --rightpercentile RIGHTPERCENTILE
                        The right percentile of the blast hits bit score distribution used to judge chimeric genes. If the blast hits bit score is greater than this percentile,
                        it will be removed.
  --nochimeras          If this option is set, the chimeric genes will not be identified.
  --nounannotated       If this option is set, the unannotated genes will not be identified.
```
### Output
`TestSp9_NewGenesOrphan.xls` Orphan genes in focus species. The other columns represent genes that are orthologous to the genes of the focal species.  
  
`TestSp9_NewGenesRepeat.xls` CNV genes in focus species. The other columns represent genes that are orthologous to the genes of the focal species. Which genes in the focus species are CNV with each other is not provided in this table. You can search for the information in the `TestSp9_RepeatList.txt` file.  
  
`TestSp9_NewGenesParalog.xls` Duplicated genes in focus species. The other columns represent genes that are orthologous to the genes of the focal species. Which genes in the focus species are duplicated with each other is not provided in this table. You can search for them in the `TestSp9_paralog_mcl.out` file, but don't forget to exclude the CNVs in the `TestSp9_RepeatList.txt` file.  
  
`TestSp9_OrthoGroupCluster.xls` Genes that are orthologous to the genes of the focus species.  
  
`TestSp9_ChimerasTrans.xls` Potential chimeric genes in the focal species. If the genes in the focus species are separated by commas, it is fission in the focus species; if there is only one gene in the focus species, it is fusion. The other columns represent the genes supporting this chimeric gene in other species.

In the table, genes containing the word `translocation` indicate that translocation may have occurred at that position, while genes containing `*_region*::*:*-*.pep*` indicate unannotated genes.  
# How to cite
