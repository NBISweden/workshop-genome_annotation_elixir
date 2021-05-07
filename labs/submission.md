---
layout: default-overview
title: Submission
exercises: 20
questions:
  - How to prepare annotation file for submission?
  - How to submit an annotation to INSDC database?
objectives:
  - be able to prepare an annotation file for submission
---

# Prerequisites
For this exercise you need to use command line.

Setup the folder structure:

```bash
conda activate agat
export data=/home/data/data_annotation/
export submission_path=~/annotation/submission
mkdir -p $submission_path
cd $submission_path
ln -s $data/genome/genome.fa
ln -s ~/annotation/functional_annotation/final_annotation.gff

```

# Submission to public repository (creation of an EMBL file)

Once your are satisfied by the wonderful annotation you have done, it is important to submit it to a public repostiroy. Fisrt, you will be applaused by the community because you share your nice work. Second, this is often mandatory if you wish to publish some work related to this annotation.

Current state-of-the-art genome annotation tools use the GFF3 format as output, while this format is not accepted as submission format by the International Nucleotide Sequence Database Collaboration (INSDC) databases. Converting the GFF3 format to a format accepted by one of the three INSDC databases is a key step in the achievement of genome annotation projects. However, the flexibility existing in the GFF3 format makes this conversion task difficult to perform.

In order to submit to **NCBI**, the use of a tool like [GAG](https://genomeannotation.github.io/GAG/) will save you lot time.  
In order to submit to **EBI**, the use of a tool like [EMBLmyGFF3](https://github.com/NBISweden/EMBLmyGFF3) will be your best choice.

In real life, prior to a submission to ENA, you need to create an account and create a project asking a locus_tag for your annotation. You have also to fill lot of metada information related to the assembly and so on. We will skip those tasks using fake information.

### Data preparation for submission to ENA (EBI)

First you need polish your annotation to filter or flag suprious cases (e.g short intron < 10 bp) otherwise the submission might fail :

```bash
agat_sp_flag_short_introns.pl --gff final_annotation.gff -o maker_final_short_intron_flagged.gff
agat_sp_fix_features_locations_duplicated.pl --gff maker_final_short_intron_flagged.gff -o maker_final_short_intron_flagged_duplicated_location_fixed.gff
```

### Data convertion

Then you will run EMBLmyGFF3 but you need first to get rid of exons that are often source of problem. Anyway the exon inforamtion will be redundant because is stored within the mRNA features.

```bash
conda activate embl
EMBLmyGFF3 --expose_translations
```
:bangbang: If it didn't work you can do
```
cp /opt/anaconda3/envs/embl/lib/python3.8/site-packages/EMBLmyGFF3/modules/translation_gff_feature_to_embl_feature.json .
```
Then modify translation_gff_feature_to_embl_feature.json to get rid of exons during the convertion.

```bash
nano translation_gff_feature_to_embl_feature.json
```
<details>
<summary>:key: Click here to see the expected translation_gff_feature_to_embl_feature.json.</summary>
  ...  
 "exon": {<br>   
   "remove": true <br>
 }<br>   
  ...   
</details>   


You can now run the conversion:

```bash
EMBLmyGFF3 maker_final_short_intron_flagged_duplicated_location_fixed.gff genome.fa -o my_annotation_ready_to_submit.embl \
        --topology linear \
        --molecule_type 'genomic DNA' \
        --transl_table 1  \
        --species 'Drosophila melanogaster' \
        --locus_tag LOCUSTAG \
        --project_id PRJXXXXXXX \
        -o my_annotation_ready_to_submit.embl
```

### Check the sanity of your embl file

If you use the Webin-CLI program (Command Line Submissions) from ENA, it contains embl-api-validator that will check the sanity of your EMBL file automatically as first step.  
If you don't use the Webin-CLI program (Interactive Submissions, Programmatic Submissions) or you just want to check the sanity of your file you can use directly the embl-api-validator. You can find it at the [ena repository](https://github.com/enasequence/sequencetools).  

Download the validator and validate your file:

```bash
wget https://repo1.maven.org/maven2/uk/ac/ebi/ena/sequence/embl-api-validator/1.1.265/embl-api-validator-1.1.265.jar
java -jar embl-api-validator-1.1.265.jar -r my_annotation_ready_to_submit.embl
```


If the file is validated, you now have a EMBL flat file ready to submit. In theory to finsish the submission, you will have to send this archived file to their ftp server and finish the submission process in the website side too.
But we will not go further. We are done. CONGRATULATIONS! You know most of the secrets needed to understand the annotations and perform your own!  
