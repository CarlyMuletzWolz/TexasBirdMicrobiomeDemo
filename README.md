# TexasBirdDemo
 Tutorial to practice microbiome processing and analyses


This pipeline is always a work in progress, so do not consider it complete, but a great place to start. I have plenty, plenty more code for various derivatives of analyses (there are many!). Please let me know if you run into issues along the way because I will continue to share this code with others, so feedback is critical and helpful for everyone involved. 
 
Please find attached files for your use:

 
### Sequence quality filtering, pair-end reads, assign taxonomy using dada2 in R

File name: Dada2_TxBird.R

I do not have access to the sequence data at this point, but have included the file after processing seqtab_TxBird.rds for you to look at and complete processing if useful.

At some point I will get them into this repository to go through my dada2 pipeline, but I think just seeing the code at least should be helpful. Please reach out to me if you would like to have those data to run the dada2 pipeline.


### Pre-processing step. Remove singletons, non-bacterial taxa, contaminants and estimate alpha diversity

File name: TxBird_preProcess.Rmd
Knitted file: TxBird_preProcess.html. You can download the .html to see what all the code looks like when it has been run.

_ _All the files you need for this generated in Step 1: TxBird_feature_table.csv (ASV/OTU table), TxBird_meta.csv (metadata file), TxBird_taxonomy.csv (taxonomy file, TxBird_feature_DNAsequences.fasta (DNA sequences of the ASVs)


### Analysis step

File name: TxBird_InitAnalysis.Rmd
Knitted file: TxBird_InitAnalysis.html. You can download the .html to see what all the code looks like when it has been run.

_ _All the files you need for this: TxBird_feature_tableCLEAN.csv (ASV/OTU table after step 2), TxBird_meta_alpha.csv (metadata file after step 2), TxBird_taxonomy.csv (taxonomy file, TxBird_feature_DNAsequences.fasta (DNA sequences of the ASVs)

