
### By: Carly Muletz Wolz

### DADA2 pipeline and creating feature table, taxonomy table and sequence file for later analyses


#############  INITIAL processing of files in terminal to get this into dada2 format #############

## From basespace, the files are downloaded with each sample having a folder 
# and within that folder are the forward and reverse reads

## To make it compatible with dada2 and R we need to copy all the fastq files to one folder 
# Then you can delete the 

## You need to navigate to the project folder (most likely called FASTQ_Generation...) in terminal 
# cd Documents/BaseSpace/TxBird-140691582/FASTQ_Generation_2019-09-07_07_01_43Z-192355164/
# 1st make directory, second move all the files within the folders to the new directory
# then remove the folders with nothing in them now

# mkdir TxBird
# mv -v *_L001*/* TxBird/
# rm -r *_L001*/


### YOU MUST UNZIP all of the files

# gunzip TxBird/*_L001*


###############  INSTALL DADA2 if you don't already have it.  ################
### Follow tutorial on how to install (follow 1. and 2.) https://benjjneb.github.io/dada2/dada-installation.html

library(dada2)
packageVersion("dada2")
## received warning:no function found corresponding to methods exports from ‘GenomicAlignments’ for: ‘concatenateObjects’


########## DADA2 tutorial is very helpful: https://benjjneb.github.io/dada2/tutorial.html ##########


## Run data is here

setwd("/Users/Carly/Documents/BaseSpace/TxBird-140691582/FASTQ_Generation_2019-09-07_07_01_43Z-192355164/TxBird/")

path <- "/Users/Carly/Documents/BaseSpace/TxBird-140691582/FASTQ_Generation_2019-09-07_07_01_43Z-192355164/TxBird/"
list.files(path)



#FILTER AND TRIM 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#INSPECTION OF QUALITY PROFILES

plotQualityProfile(fnFs[8:16])

plotQualityProfile(fnRs[8:16])


#FILTERING AND TRIMMING 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


## parameters for filtering data
out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,180),
                     maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23, 
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out1)
str(out1)
mean(out1[,2])
mean(out1[,1])

## on average 65% with this 
mean(out1[,2])/mean(out1[,1])


## NEED to run whichever parameter last so that filtFs and filtRs are from this
out3 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,190),
                      maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23,
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)
head(out3)
mean(out3[,2])

## on average 65% with this 
mean(out3[,2])/mean(out3[,1])

## GOING with out1


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

setwd("/Users/Carly/Dropbox (Smithsonian)/SI_projects/TX_migBirds/Analyses/")

## Used out 1 parameters
saveRDS(seqtab, "seqtab_TxBird.rds")

## If want to read back in 
##seqtab <- readRDS("seqtab_run1.rds")

dim(seqtab)

## Distribution of amplicon sizes in bp
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)


## https://benjjneb.github.io/dada2/training.html
## Assign taxonomy is up three directories so that I can use these files for multiple projects
taxa <- assignTaxonomy(seqtab.nochim, "/Users/Carly/Documents/BaseSpace/Assign_Taxonomy_Dada2/rdp_train_set_16.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/Carly/Documents/BaseSpace/Assign_Taxonomy_Dada2/rdp_species_assignment_16.fa.gz")

#  inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


library(phyloseq); packageVersion("phyloseq")

## combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps

## rename ASVs to numbers
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format

library(seqRFLP)
## should work, but, below can also convert .csv file to fasta
## doesn't seem to work for merged data
## seq_data <- dataframe2fas(st.all, file = 'SalAMP_DNAsequences.fasta')


## convert feature table to matrix
site_species <-as(otu_table(ps), "matrix")

## need to change this to match mapping file later
rownames(site_species)

# not sure what this is doing. needed it for another project. Might be useful for underscores?
samples.out <- rownames(site_species)
rownames(site_species) <- sapply(strsplit(samples.out, "f"), `[`, 1)
rownames(site_species)

## transpose to make a species by site matrix

species_site <- t(site_species)

# taxon table 
tax <- as(tax_table(ps), "matrix")
seqs


getwd()
setwd("/Users/Carly/Dropbox (Smithsonian)/SI_projects/TX_migBirds/Analyses/")

## Write this file out and look at it. Determine what to do with ASVs in Neg Controls
## select whole worksheet and filter by reads in neg controls, 
## MAKE SURE you have selected all columns/rows before you filter

## Remove ASVs if they are in all negative controls or in one negative control, but almost all samples
## pretty subjective, no guidelines really on what to do here
## Once done making neg control decisions, remove negative controls from file 
## and rename _final.csv
write.csv(species_site, "TxBird_feature_table.csv")

write.csv(tax, "TxBird_taxonomy.csv")
write.csv(seqs, 'TxBird_feature_DNAsequences.csv')

library(seqRFLP)
## you want to generate a fasta file of sequences for reading back in 
seq_data <- dataframe2fas(seqs, file = "TxBird_feature_DNAsequence.fasta")


