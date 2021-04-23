---
layout: default-overview
title: Training ab-initio predictor
exercises: 60
questions:
  - What do I need to train Augustus?
  - What are the step to train Augustus?
objectives:
  - Going through the steps
  - Understanding the different part of the training (it is complex so take your time!)
---

# Prerequisites

Setup the folder structure:

```bash
export data=/home/data/data_annotation/
export augustus_training_path=~/annotation/structural_annotation/augustus_training
export maker_evidence_path=~/annotation/structural_annotation/maker_evidence
cd
mkdir -p $augustus_training_path
cd $augustus_training_path
```

# Introduction

From the maker run evidence based, we can train our ab-initio predictors and then use them for the second run of annotation.
You will need a set of genomic sequences with gene structures (sequence coordinates of starts and ends of exons and genes) and the most important part is to select the right set of genes.
In many cases, or as a first step towards modeling complete genes, it is sufficient to have only the coding parts of the gene structure (CDS).
We will only train augustus today as it is one the best ab-initio predictor and one of the hardest to train.
Maker also support SNAP (Works good, easy to train, not as good as others ab-initio especially on longer intron genomes), GeneMark (Self training, no hints, buggy, not good for fragmented genomes or long introns), FGENESH (Works great, costs money even for training) and now EVM.


## Training Augustus

You will need to symlink the evidence-based annotation (the gff annotation file from the first run of maker) and the genome fasta sequence.

```bash
ln -s $maker_evidence_path/maker_evidence/maker_annotation.gff maker_evidence.gff
ln -s $data/genome/genome.fa
```

## Compile a set of training and test genes

* First of all, as we will generate lot of files we create structured folders to save them in a sorted way.
(you can check that you are in the proper folder by doing pwd, you should be in augustus_training/)

```bash
mkdir filter  
mkdir protein  
mkdir nonredundant  
mkdir blast_recursive  
mkdir gff2genbank  
```

* Then, we select only the protein coding genes from the maker.gff file and remove all other type of level2 features: tRNA, snoRNA, etc. ( :bulb: **Tips**: in this case we only have mRNA)

```bash
conda deactivate
conda activate agat
agat_sp_separate_by_record_type.pl -g maker_evidence.gff -o maker_results_noAbinitio_clean
ln -s maker_results_noAbinitio_clean/mrna.gff
```

* Next step, we need to filter the best genes we will use for the training, we need complete genes and a AED score under 0.3 (:bulb: **Tips**:those are our default parameters you can change them if you want to be more selective with an AED under 0.1, you can also set a distance between genes if you know your genes are further appart).

```
agat_sp_filter_feature_by_attribute_value.pl --gff mrna.gff --value 0.3 -a _AED -t ">" -o filter/codingGeneFeatures.filter.gff
```

* Now we filter mRNAs to keep only the longest isoform when several are avaialble.

```
agat_sp_keep_longest_isoform.pl -f filter/codingGeneFeatures.filter.gff -o filter/codingGeneFeatures.filter.longest_cds.gff
```

* We then check that the gene models are complete (has a start and stop codon), and remove the incomplete ones.

```
agat_sp_filter_incomplete_gene_coding_models.pl --gff filter/codingGeneFeatures.filter.longest_cds.gff -f genome.fa -o filter/codingGeneFeatures.filter.longest_cds.complete.gff
```

* We may also filter the gene models by distance from neighboring genes in order to be sure to have nice intergenic region wihtout genes in it in order to train properly intergenic parameters. (default 500 bp)

```
agat_sp_filter_by_locus_distance.pl --gff filter/codingGeneFeatures.filter.longest_cds.complete.gff -o filter/codingGeneFeatures.filter.longest_cds.complete.good_distance.gff
```

* To avoid to bias the training and give an exhaustive view of the diversity of gene we have to remove the ones that are too similar to each other. In order to do so, we translate our coding genes into proteins, format the protein fasta file to be able to run a recursive blast and then filter them.

```
agat_sp_extract_sequences.pl -o protein/codingGeneFeatures.filter.longest_cds.complete.good_distance_proteins.fa -f genome.fa -p -cfs -cis -ct 1 --g filter/codingGeneFeatures.filter.longest_cds.complete.good_distance.gff
```

```
conda deactivate
conda activate blast
makeblastdb -in protein/codingGeneFeatures.filter.longest_cds.complete.good_distance_proteins.fa -dbtype prot -out protein/codingGeneFeatures.filter.longest_cds.complete.good_distance.proteins


blastp -query protein/codingGeneFeatures.filter.longest_cds.complete.good_distance_proteins.fa -db protein/codingGeneFeatures.filter.longest_cds.complete.good_distance.proteins -outfmt 6 -out blast_recursive/codingGeneFeatures.filter.longest_cds.complete.good_distance.proteins.fa.blast_recursive
```
```
conda deactivate
conda activate agat
agat_sp_filter_by_mrnaBlastValue.pl --gff filter/codingGeneFeatures.filter.longest_cds.complete.good_distance.gff --blast blast_recursive/codingGeneFeatures.filter.longest_cds.complete.good_distance.proteins.fa.blast_recursive --outfile nonredundant/codingGeneFeatures.nr.gff

```

* Sequences need to be converted to a simple genbank format.

```
conda deactivate
conda activate bioingo
/opt/anaconda3/envs/bioinfo/scripts/gff2gbSmallDNA.pl nonredundant/codingGeneFeatures.nr.gff $data/genome/genome.fa 500 gff2genbank/codingGeneFeatures.nr.gbk
```

* You should split the set of gene structures randomly.
In theory in order for the test accuracy to be statistically meaningful the test set should also be large enough (100-200 genes).

```bash
/opt/anaconda3/envs/bioinfo/scripts/randomSplit.pl gff2genbank/codingGeneFeatures.nr.gbk 100
```

:question:What happened? how can you solve it? what might be the consequences of it?

<details>
<summary>:key: Click to see the solution .</summary>
There are not 100 genes in the file, because we are using only the chr4 of drosophila.
The training will probably not be good!
</details>

## Train Augustus

Now that you have created a set of gene to train augustus, let's train it!

You need to copy the config folder of augustus in your folder (to be able to write in it)

```
cp -r /opt/anaconda3/envs/bioinfo/config/ .
export AUGUSTUS_CONFIG_PATH=~/annotation/structural_annotation/augustus_training/config
```

:bangbang: **Replace <$USER> by your username.** 

```bash
/opt/anaconda3/envs/bioinfo/scripts/new_species.pl --species=dmel_$USER

etraining --species=dmel_$USER gff2genbank/codingGeneFeatures.nr.gbk.train

augustus --species=dmel_$USER gff2genbank/codingGeneFeatures.nr.gbk.test | tee run.log 
```

- Look at the accuracy report, what does it mean? why? see [Training Augustus](http://www.vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html)
