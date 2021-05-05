---
layout: default-overview
title: Abinitio with augustus
exercises: 45
questions:
  - How to run augustus?
  - How visualise you results?
objectives:
  - Run augustus
  - Open and load data in webapollo
---


# Prerequisites

Setup the folder structure:

```bash
export data=/home/data/data_annotation/
export abinitio_augustus_path=~/annotation/structural_annotation/abinitio_augustus
cd
mkdir -p $abinitio_augustus_path
cd $abinitio_augustus_path
```

# Running an ab initio gene finder

<u><strong>Ab initio gene finders:</strong></u> These methods have been around for a very long time, and there are many different programs to try. We will in this exercise focus on the gene finder Augustus. These gene finders use likelihoods to find the most likely genes in the genome. They are aware of start and stop codons and splice sites, and will only try to predict genes that follow these rules. The most important factor here is that the gene finder needs to be trained on the organism you are running the program on, otherwise the probabilities for introns, exons, etc. will not be correct. Luckily, these training files are available for Drosophila.

:mortar_board: **Augustus:**

Augustus comes with a list of species that already have a trained hmm model available. You can have a look at this list like that:  

```
conda activate bioinfo
augustus --species=help
```

:question:Did you see the appropriate model for Drosophila Melanogaster?
<details>
<summary>:key: Click to see the solution .</summary>
<code>
fly                                      | Drosophila melanogaster
</code>
</details>


So, let's now launch Augustus on our genome with the `fly` model.

```
augustus --species=fly $data/genome/genome.fa --gff3=yes --progress=true > augustus_drosophila.gff
```

if you wish to annotate isoforms too, use the following command:

```
augustus --species=fly $data/genome/genome.fa --gff3=yes --progress=true --alternatives-from-sampling=true > augustus_drosophila_isoform.gff
```

Take a look at the result file using ‘less augustus\_drosophila.gff’. What kinds of features have been annotated? Does it tell you anything about UTRs?

<details>
<summary>:key: Click to see the solution .</summary>
You have annotated genes, transcripts, exons, CDS, introns...<br/>
No UTR for drosophila melanogaster but can do it for other species see <A href="http://bioinf.uni-greifswald.de/augustus/"> augustus</A>
</details>


The gff-format of Augustus is non-standard (looks like gtf) so to view it in a genome browser you need to convert it. You can do this using the following command line:

```
conda deactivate
conda activate agat
agat_convert_sp_gxf2gxf.pl -g augustus_drosophila.gff -o augustus_drosophila.gff3
```
To better understand what contains your gff file you may use a script that will provide you some statistics like this one:
```
agat_sp_statistics.pl --gff augustus_drosophila.gff3
```
:question:How many genes have you annotated?
<details>
<summary>:key: Click to see the solution .</summary>

Compute transcript with isoforms if any <br/>

Number of genes                              60<br/>
Number of transcripts                        60<br/>
Number of mrnas with utr both sides          60<br/>
Number of mrnas with at least one utr        60<br/>
Number of cdss                               60<br/>

</details>


Transfer the augustus\_drosophila.gff3 to your computer using scp from a local terminal:    
```
scp -P <port_number> <login_name>@<host>:~/annotation/structural_annotation/abinitio_augustus/augustus_drosophila.gff3 .  
```
Load the file in [Webapollo (Here find the instruction)](webapollo_usage)
<br/>The official Ensembl annotation is available in the genome browser.  
:question: How does the Augustus annotation compare with the Ensembl annotation? Are they identical?

<details>
<summary>:key: Click to see the solution .</summary>
No they are not!<br/>
You can notice that augustus has less genes than Ensembl and also that often genes from augustus are longer and contain several genes from Ensembl, this shows the importance of manual curation and having external evidence for annotation.
</details>

:mortar_board: **Augustus with yeast models:**  
Run augustus on the same genome file but using settings for yeast instead (change species to Saccharomyces).

<details>
<summary>:key: Click to see the solution .</summary>
<pre class="code">
conda deactivate
conda activate bioinfo
augustus --species=saccharomyces $data/genome/genome.fa --gff3=on > augustus_saccharomyces.gff
</pre>
</details>

Have a look at the statistics to have a first impression of what are the differences compared to the previous annotation (like you did previously).

<details>
<summary>:key: Click to see the solution .</summary>
<pre class="code">
conda deactivate
conda activate agat
agat_sp_statistics.pl --gff augustus_saccharomyces.gff
</pre>
</details>


Load this result file into Webapollo and compare with your earlier results.  

:question: Can you based on this draw any conclusions about how a typical yeast gene differs from a typical Drosophila gene?

<details>
<summary>:key: Click to see the solution .</summary>
You can see that yeast genes are smaller, often one exon genes, if they have several exons the introns are smaller than for drosophila.
</details>

# Closing remarks

We have seen how to launch a quick annotation using the abinitio tool Augustus.
We have also seen the importance to use a species specific hmm model into the ab initio tool. Thus, the limitation of this approach is linked to the pre-trained species that are available.
