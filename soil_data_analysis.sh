#!/bin/bash
## This script is used to analyze 16S data generated in Jones et al. 2019 paper. 
## Author Zhen Fu (fuzhen327@gmail.com)
## Illumina HiSeq was used to sequence 85 soil samples in vegetable farms in the western USA, including organic, integrated with livestock, and conventional farms.
## Qiime1 was used to analyze the data.

## The sequences can be downloaded from NCBI. However, NCBI-SRA does merge the paired-end reads into one file, which is difficult to separate.
## Therefore, we highly recommend downloading the orginal sequences from ENA (https://www.ebi.ac.uk/ena). Simply search "SRP126862", a list of all 85 
## samples(runs) should show up. There is a "Bulk Download File" icon, and a java-based tool will help users to download.
## I personally have tried to use a browser on Linux system and Windows, and I found that Windows worked
## better for unkonwn reason. If you have any question about getting sequneces, please feel free to contact the author. 

## The follow loop does three tasks: 1). unzip the fastq files; 2) rename the forward and reverse sequences; 3). join the forward and reverse reads using Qiime. 
## The loop takes care of the naming as far as "Run_table_soil_samples" is under the working directory. "Run_table_soil_samples" is in GitHub. 

while read LINE; do

  a=$(echo "$LINE" |cut -f6)
  b=$(echo "$LINE" |cut -f8)
  echo "${a}"
  forward_seq=$a'_1.fastq.gz'
  reverse_seq=$a'_1.fastq.gz'
  forward_seq=$a'_1.fastq'
  reverse_seq=$a'_2.fastq'

  new_name_f=$b'_1.fq'
  new_name_r=$b'_2.fq'

gunzip ${forward_seq} && mv ${forward_seq} ${new_name_f}
gunzip ${reverse_seq} && mv ${reverse_seq} ${new_name_r}

join_paired_ends.py -f ${new_name_f} -r ${new_name_r} -o ${b} -m fastq-join

  done < Run_table_soil_samples


## Next, we will use split_libraries_fastq.py. Orignally, this script in Qiime is used to split the libraries based on barcodes. Here, we use this to 
## geneate "seqs.fna", which is a fasta formatted file where each sequence is renamed. Additionally, the advantage of this script is it also trims
## and checks the quality of our merged fastq sequences. 

split_libraries_fastq.py -i \
CA.1.ORG.2014.join.fq,CA.2.ORG.2014.join.fq,CA.1.INT.2014.join.fq,CA.2.INT.2014.join.fq,CA.3.INT.2014.join.fq,CA.3.ORG.2014.join.fq,CA.4.INT.2014.join.fq,CA.4.ORG.2014.join.fq,\
CA.5.ORG.2014.join.fq,CA.5.INT.2014.join.fq,OR.1.ORG.2014.join.fq,OR.1.INT.2014.join.fq,OR.2.INT.2014.join.fq,OR.2.CONV.2014.join.fq,WA.1.INT.2014.join.fq,OR.3.CONV.2014.join.fq,\
OR.3.ORG.2014.join.fq,OR.4.INT.2014.join.fq,OR.4.CONV.2014.join.fq,OR.4.ORG.2014.join.fq,OR.5.CONV.2014.join.fq,OR.5.INT.2014.join.fq,OR.5.ORG.2014.join.fq,WA.1.CONV.2014.join.fq,\
OR.3.INT.2014.join.fq,WA.1.ORG.2014.join.fq,WA.2.INT.2014.join.fq,WA.3.INT.2014.join.fq,WA.3.ORG.2014.join.fq,WA.4.CONV.2014.join.fq,WA.4.INT.2014.join.fq,WA.4.ORG.2014.join.fq,\
WA.5.ORG.2014.join.fq,WA.5.INT.2014.join.fq,WA.1.CONV.2015.join.fq,WA.1.INT.2015.join.fq,WA.1.ORG.2015.join.fq,WA.3.INT.2015.join.fq,WA.3.ORG.2015.join.fq,WA.4.CONV.2015.join.fq,\
WA.4.INT.2015.join.fq,WA.4.ORG.2015.join.fq,WA.5.INT.2015.join.fq,WA.5.ORG.2015.join.fq,OR.1.INT.2015.join.fq,OR.1.ORG.2015.join.fq,OR.5.ORG.2015.join.fq,OR.2.CONV.2015.join.fq,\
OR.2.INT.2015.join.fq,OR.3.CONV.2015.join.fq,OR.3.INT.2015.join.fq,OR.3.ORG.2015.join.fq,OR.4.CONV.2015.join.fq,OR.4.INT.2015.join.fq,OR.4.ORG.2015.join.fq,OR.5.CONV.2015.join.fq,\
OR.5.INT.2015.join.fq,CA.2.CONV.2015.join.fq,CA.2.INT.2015.join.fq,CA.2.ORG.2015.join.fq,CA.3.CONV.2015.join.fq,CA.3.INT.2015join.fq,CA.3.ORG.2015join.fq,CA.4.INT.2015.join.fq,\
CA.4.ORG.2015.join.fq,CA.5.INT.2015.join.fq,CA.5.ORG.2015.join.fq,CA.1.CONV.2015.join.fq,CA.1.ORG.2015.join.fq,CA.2.CONV.2016.join.fq,CA.2.ORG.2016.join.fq,CA.4.INT.2016.join.fq,\
CA.3.CONV.2016.join.fq,CA.3.INT.2016.join.fq,CA.1.ORG.2016.join.fq,OR.4.ORG.2016.join.fq,OR.1.ORG.2016.join.fq,OR.3.INT.2016.join.fq,OR.5.CONV.2016.join.fq,OR.4.INT.2016join.fq,\
OR.2.CONV.2016join.fq,WA.4.ORG.2016join.fq,WA.4.CONV.2016.join.fq,WA.3.INT.2016.join.fq,WA.5.INT.2016.join.fq \
-o split_library --sample_ids \
CA.1.ORG.2014,CA.2.ORG.2014,CA.1.INT.2014,CA.2.INT.2014,CA.3.INT.2014,CA.3.ORG.2014,CA.4.INT.2014,CA.4.ORG.2014,CA.5.ORG.2014,CA.5.INT.2014,OR.1.ORG.2014,OR.1.INT.2014,\
OR.2.INT.2014,OR.2.CONV.2014,WA.1.INT.2014,OR.3.CONV.2014,OR.3.ORG.2014,OR.4.INT.2014,OR.4.CONV.2014,OR.4.ORG.2014,OR.5.CONV.2014,OR.5.INT.2014,OR.5.ORG.2014,WA.1.CONV.2014,\
OR.3.INT.2014,WA.1.ORG.2014,WA.2.INT.2014,WA.3.INT.2014,WA.3.ORG.2014,WA.4.CONV.2014,WA.4.INT.2014,WA.4.ORG.2014,WA.5.ORG.2014,WA.5.INT.2014,WA.1.CONV.2015,WA.1.INT.2015,\
WA.1.ORG.2015,WA.3.INT.2015,WA.3.ORG.2015,WA.4.CONV.2015,WA.4.INT.2015,WA.4.ORG.2015,WA.5.INT.2015,WA.5.ORG.2015,OR.1.INT.2015,OR.1.ORG.2015,OR.5.ORG.2015,OR.2.CONV.2015,\
OR.2.INT.2015,OR.3.CONV.2015,OR.3.INT.2015,OR.3.ORG.2015,OR.4.CONV.2015,OR.4.INT.2015,OR.4.ORG.2015,OR.5.CONV.2015,OR.5.INT.2015,CA.2.CONV.2015,CA.2.INT.2015,CA.2.ORG.2015,\
CA.3.CONV.2015,CA.3.INT.2015,CA.3.ORG.2015,CA.4.INT.2015,CA.4.ORG.2015,CA.5.INT.2015,CA.5.ORG.2015,CA.1.CONV.2015,CA.1.ORG.2015,CA.2.CONV.2016,CA.2.ORG.2016,CA.4.INT.2016,\
CA.3.CONV.2016,CA.3.INT.2016,CA.1.ORG.2016,OR.4.ORG.2016,OR.1.ORG.2016,OR.3.INT.2016,OR.5.CONV.2016,OR.4.INT.2016,OR.2.CONV.2016,WA.4.ORG.2016,WA.4.CONV.2016,WA.3.INT.2016,\
WA.5.INT.2016 \
--barcode_type 'not-barcoded' --phred_offset 33 -q 19


## The above command will create a giant file called "seqs.fna", which is a fasta file. 
## next we can pick OTU based on this giant fasta file: 
pick_otus.py -i ./split_library/seqs.fna -o picked_otus_default

## Pick a representative sequence for each OTU
pick_rep_set.py -i ./picked_otus_default/seqs_otus.txt -f ./split_library/seqs.fna -o rep_seqs_pickOTU_default

## Assign taxonomy to OTU representative sequences: 
assign_taxonomy.py -i rep_seqs_pickOTU_default -o taxonomy_results

## Align OTU representative sequences:
align_seqs.py -i rep_seqs_pickOTU_default -o alignment/

## Filter the alignment:
filter_alignment.py -i ./alignment/rep_seqs_pickOTU_default_aligned.fasta -o filtered_alignment/ --remove_outliers

## Build a phylogenetic tree:
make_phylogeny.py -i ./filtered_alignment/rep_seqs_pickOTU_default_aligned_pfiltered.fasta -o rep_set_tree

## Make the OTU table
make_otu_table.py -i ./picked_otus_default/seqs_otus.txt -t ./taxonomy_results/rep_seqs_pickOTU_default_tax_assignments.txt -o otu_table.biom

## From the step above, it generates a file called biom. This is a binary file, you will need program biom to convert this binrary file to a format similar to CSV/TSV.
## Make sure you have biom installed. 
biom convert -i otu_table.biom -o otu_table_tabseparated.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"

## Then you can generate a summary table
summarize_taxa.py -i otu_table.biom -o taxonomy_summaries
## This table records the relative abundance and taxonomy of each OTU across all your samples.


## Next, we can compute alpha diversity and generate alpha rarefaction
## Since alpha diversity includes many measurements, we can choose which ones we would like Qiime to estimate
## Here, this infomraiton is included in a file called "alpha_params.txt", which has one line of information as below:
## alpha_diversity:metrics chao1,chao1_ci,dominance,equitability,observed_otus,observed_species,simpson_reciprocal,shannon,simpson,simpson_e,PD_whole_tree
## Last file you will need is a mapping file; it recorded information of our samples, such as year, state, and farming types. This file in on Github
alpha_rarefaction.py -i otu_table.biom -m mapping_soil.tsv -o alpha -p alpha_params.txt -t rep_set.tre

## Rarefaction, since the lowest coverage for our samples is ~58,000 reads. So we will do a rarefaction curve starting from 10,000 and goes up to 58,000
alpha_rarefaction.py -i otu_table.biom -m pre_mapping.tsv -o alpha_rare -p alpha_params.txt --min_rare_depth 10000 --max_rare_depth 58000 --num_steps 5 -t rep_set.tre














