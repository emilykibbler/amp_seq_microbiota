## Intro and package loading ---------------------------

setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/crypto_data") 

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/load_data.R")
# source("AVS554_packages.R")
# install_necessary()
# install_optional()

library(tidyverse)
library(devtools)
library(Rqc)
library(BiocParallel)
library(Rqc)
library(Biostrings)
library(dada2); packageVersion('dada2')
library(beepr)
library(phyloseq)
library(decontam)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(ggtext)
library(lme4)
library(lmerTest)
library(emmeans)
library(remotes)
library(PerformanceAnalytics)
library(asbio)
library(microbiome)
library(ggpubr)


setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/crypto_data") 

## Lab 2: import and filter ----------------------

# "Types" of data could be forward and reverse, batches, etc.
# Code expects each "type" to have a corresponding folder with format type_raw in the wd

file_info <- metadata_import(types = c("forward", "reverse"))


# One plot for each data type will be in the wd in .png format
# Default data type is raw
qualplots(file_info); beep()
# Forward starts with 3610171 reads
# Reverse starts with 3610171 reads as well

# this didn't work, come back later to troubleshoot
for (i in 1:nrow(file_info)) {
  rqc(path = file_info$raw_data_dir[i][1], 
      sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
      file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())
}; beep()

# F and R both look good. I will trim both at 250 (300-50)

filtFs <- file.path(file_info$filt_dir[1], paste0(file_info$sample_names[[1]], "_F_filt.fastq.gz"))
filtRs <- file.path(file_info$filt_dir[2], paste0(file_info$sample_names[[2]], "_R_filt.fastq.gz"))


# filter that sample pair and output as new same names, run this whole chunk from filtoutput through verbose=TRUE)	  
filtoutput <- filterAndTrim( #filter command
  file_info$file_names_full[[1]], filtFs, file_info$file_names_full[[2]], filtRs, # specifies in/out variables
  trimLeft = c(10, 10), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
  trimRight = c(50,50), # cuts of last XX bases from F,R reads, or hash this out and use truncLen
  #truncLen=c(330, 330), # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
  maxEE = c(2,3), ## The maxEE parameter sets the maximum number of expected errors allowed in a read. Always >1
  verbose = TRUE); beep() # verbose = print to screen
saveRDS(filtoutput, "filtoutput.rds")

qualplots(file_info, type = "filtered"); beep()

ggplot(filtoutput) +
  geom_point(aes(row.names(filtoutput), reads.in, color = "Raw reads")) +
  geom_point(aes(row.names(filtoutput), reads.out, color = "Filtered reads")) +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  theme(legend.title = element_blank()) +
  ylab("Reads") +
  xlab("Sample") +
  ggtitle("Crypto data, raw reads and filtered reads")
ggsave("refiltered_crypto_reads_in_out.png")


DNA_to_read <- readDNAStringSet(file_info$file_names_full[[1]][1], format = "fastq", with.qualities = TRUE)
# Pull out a sequence I can BLAST
toString(DNA_to_read[1])
# Coming up as cow mitochondrial?

## DADA2 pick sequence variants (Lab 3) ----------------------

file_info <- readRDS("file_info.rds")

create_and_plot_error_profile(file_info, bases = 1e8) ; beep()

errF <- readRDS("error_profile_forward.RDS")
errR <- readRDS("error_profile_reverse.RDS")

dadaFs <- dada(list.files(file_info$filt_dir[1], full.names = TRUE),
               err = errF, 
               multithread = TRUE, verbose = FALSE) #; beep()
dadaRs <- dada(list.files(file_info$filt_dir[2], full.names = TRUE),
               err = errR, 
               multithread = TRUE, verbose = FALSE); beep()




mergers <- mergePairs(dadaF = dadaFs,
                      derepF = list.files(file_info$filt_dir[1], full.names = TRUE),
                      dadaR = dadaRs,
                      derepR = list.files(file_info$filt_dir[2], full.names = TRUE),
                      verbose = FALSE) ; beep()
seqtab <- makeSequenceTable(mergers)

# 3. Get the dimensions of your table, which is a great way to count how many samples are and SVs included.
dim(seqtab) 


saveRDS(seqtab, '1e8_merged_seqtab.rds') 
# Ran again with 1e6 as bases variable
saveRDS(seqtab, 'merged_seqtab_1e6bases.rds') 

summary(str_length(colnames(seqtab)))

# for (i in 1:nrow(seqtab)) {
#   print(rownames(seqtab)[i])
#   print(summary(seqtab[i,]))
# }


seqtab_1e6 <- readRDS("merged_seqtab_1e6bases.rds")
seqtab_1e8 <- readRDS("1e8_merged_seqtab.rds")


# Side quest: Check if the error calculation difference makes any difference in reads per sample ----------------
total_reads_1e6 <- rowSums(seqtab_1e6)
total_reads_1e8 <- rowSums(seqtab_1e8)
df <- as.data.frame(cbind(total_reads_1e6, total_reads_1e8)) 
df$sample <- row.names(df)
df %>% ggplot() +
  geom_point(aes(x = sample, y = total_reads_1e6), color = "purple", position = "jitter") +
  geom_point(aes(x = sample, y = total_reads_1e8), color = "pink", position = "jitter") +
  ylab("Total reads") +
  ggtitle("Comparing read totals with different error calculations")
# No difference so that's good

SVs_found_by_sample_1e6 <- list()
SVs_found_by_sample_1e8 <- list()

for (i in 1:nrow(seqtab_1e6)) {
  print(paste("Number of SVs found in", rownames(seqtab_1e6)[i], ":"))
  print(length(seqtab_1e6[i,][seqtab_1e6[i,] > 0]))
  SVs_found_by_sample_1e6 <- c(SVs_found_by_sample_1e6, length(seqtab_1e6[i,][seqtab_1e6[i,] > 0]))
}
for (i in 1:nrow(seqtab_1e8)) {
  print(paste("Number of SVs found in", rownames(seqtab_1e8)[i], ":"))
  print(length(seqtab_1e8[i,][seqtab_1e8[i,] > 0]))
  SVs_found_by_sample_1e8 <- c(SVs_found_by_sample_1e8, length(seqtab_1e8[i,][seqtab_1e8[i,] > 0]))
}

df <- as.data.frame(cbind(SVs_found_by_sample_1e6, SVs_found_by_sample_1e8))
df$sample <- row.names(seqtab_1e6)
df$SVs_found_by_sample_1e6 <- as.numeric(df$SVs_found_by_sample_1e6)
df$SVs_found_by_sample_1e8 <- as.numeric(df$SVs_found_by_sample_1e8)

df %>% ggplot() +
  geom_point(aes(x = sample, y = SVs_found_by_sample_1e6), color = "purple", position = "jitter") +
  geom_point(aes(x = sample, y = SVs_found_by_sample_1e8), color = "pink", position = "jitter")
# maybe a little higher with 1e8


# df %>% ggplot() +
#   geom_boxplot(aes(y = SVs_found_by_sample_1e6), color = "purple") +
#   geom_boxplot(aes(y = SVs_found_by_sample_1e8), color = "pink")

View(df)

df <- data.frame(matrix(ncol = 3, nrow = 23))
colnames(df) <- c("sample", "err_type", "SVs")
df$sample <- row.names(seqtab_1e6)
df$err_type <- "1e6"
df$SVs <- unlist(SVs_found_by_sample_1e6)

temp <- data.frame(matrix(ncol = 3, nrow = 23))
colnames(temp) <- c("sample", "err_type", "SVs")
temp$sample <- row.names(seqtab_1e8)
temp$err_type <- "1e8"
temp$SVs <- unlist(SVs_found_by_sample_1e8)
identical(df$sample, temp$sample)

df <- rbind(df, temp)
df %>% ggplot() +
  geom_boxplot(aes(x = err_type, y = SVs))

dim(seqtab_1e6)
dim(seqtab_1e8)
summary(unlist(SVs_found_by_sample_1e6))
summary(unlist(SVs_found_by_sample_1e8))

# End of side quest: Bases used for error calculation does not make a big difference  ----------------
# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/crypto_data")

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")

seqtab <- readRDS("merged_seqtab_1e6bases.rds")

SVs_found_by_sample <- make_SV_summary(seqtab)
View(SVs_found_by_sample)
# Most SVs: 	47_R1_F_filt.fastq.gz -- 1758
# Most reads: 40_R1_F_filt.fastq.gz -- 221825
# Fewest SVs: 37_R1_F_filt.fastq.gz -- 167
# Fewest reads: 58_R1_F_filt.fastq.gz -- 27324


# DADA2 Remove chimeras from seqtab (Lab 4) ----------------------------	

# Chimeras are accidentally created sequences during lab protocols. Remove them.	

# Necessary packages -- should be loaded from the top of this script
    # library(dada2); packageVersion('dada2')
    # library(tidyverse)


# 1. Remove chimeras. Leave the method as consensus. multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE) # a few moments to run
# These samples have 9121 bimeras out of 12845 input sequences.

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab.nochim) 
dim(seqtab)
# 23 samples and 3724 SVs (dowm from 12845)

# Calculate the percentage of chimeras identified out of the total
sum(seqtab.nochim)/sum(seqtab) # 0.732 retained
1 - sum(seqtab.nochim)/sum(seqtab) # meaning 0.248 retained

# 3. Save your chimera-free sequence table.
saveRDS(seqtab.nochim, 'seqtab.nochim.rds')

#reload if needed
# seqtab.nochim <- readRDS('seqtab.nochim.rds')


## Workflow verification steps  (Lab 4) --------------------------
# This section can be optional, but as you are learning, or playing with new data, it is helpful to run some internal checks to assess whether you like how your quality control steps worked.

# 1. Track reads through the analysis here to see how many were lost at each QC step, in case you were too stringent. Only need to count F reads, even if you did paired F/R
#set function to count the unique reads
getN <- function(x) {
  sum(getUniques(x)) 
}

# load the filtered output file if necessary
# filtoutput <- readRDS("filtoutput.rds")

# 2. Load the metadata file that goes with your sequencing data so you can match factors to seq data. This uses a datatable made in Excel and saved as a .csv file
meta <- read.csv("Crypto_bacteria_metadata.csv", 
                 header = TRUE, # specifies that the top row is column names
                 nrow = 23, 
                 row.names = 1) # row.names specifies which column to call the samples by (sample name), doesn't have to be the first column
meta <- rename(meta, "Sample_source" = "Sample_type")



# 2.5 check the dimensions of the three data files you need for this to make sure the number of rows matches in each. If they do not, you may need to add/remove rows from your metadata file in case samples were removed/retained from your dataset.
dim(filtoutput) # 23 x 2
dim(seqtab.nochim) # 23 x 3724
dim(meta) # 23 x 11


head(row.names(filtoutput))
head(row.names(seqtab.nochim))
head(row.names(meta))

# Need to strip extra labeling

row.names(filtoutput) <- str_remove_all(row.names(filtoutput), "_R1.fastq.gz")
row.names(seqtab.nochim) <- str_remove_all(row.names(seqtab.nochim), "_R1_F_filt.fastq.gz")
saveRDS(seqtab.nochim,"seqtab_nochim.rds")
identical(row.names(filtoutput), row.names(seqtab.nochim)) # true, yay
identical(row.names(meta), row.names(filtoutput)) # true


# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab.nochim), meta$Sample_type) 
# 4. Assign column names to that new variable
colnames(track) <- c("reads.in","filtered", "nonchimeras", "Sample_type") 

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'tracked_seqs.RDS') 

# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))


# plot with Sample_type along the X axis
ggplot(plotData,aes(x = Sample_type, y = as.numeric(totals))) + geom_jitter(aes(color = type)) + 
  ylab("Sequences") + 
  xlab("Sample_type") +
  theme(axis.text.x = element_text(angle = 25, size = 8, hjust = 1)) +
  ggtitle("Reads by sample type")
ggsave("reads_filt_chim_sample_type.png")

# or, plot with QA stage along the X axis
ggplot(plotData,aes(x = type,y = as.numeric(totals))) + geom_jitter(aes(color = Sample_type)) + 
  ylab("Sequences") + 
  xlab("QA stage") +
  ggtitle("Reads by filtering step")
ggsave("reads_sample_type_QC_status.png")








# 9. Randomly select a negative control sample to see what's in it
unqs.NC <- seqtab.nochim["59",]

# Sort by # seqs and drop SVs absent in the negative control sample
unqs.NC <- sort(unqs.NC[unqs.NC > 0], decreasing = TRUE) 

# Print out how many reads are inferred in the negative control
cat("DADA2 inferred", length(unqs.NC), "sequence variants present in the selected sample.\n")
# DADA2 inferred 245 sequence variants present in the selected sample.

# Plot the number of sequences in the SVs found in the negative control sample.
plot(unqs.NC, ylab = "Number of seqs/SV, Neg Control", xlab = "SVs", main = "Number of sequences/SV in Negative Control") 

# OPTIONAL: add taxonomy and see what's in the negative control, for example if you want to know if you have the same contaminant showing up in negative controls over time.
# taxa.NC <- assignTaxonomy(unqs.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
# write.csv(taxa.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_NC_EXAMPLE.csv') 



# 10. Randomly select a positive control to see what's in it
unqs.PC <- seqtab.nochim["47",] #CHANGE ME, update with sample name of a positive control, such as Mock_S192

# Sort by # seqs and drop SVs absent in the positive control sample
unqs.PC <- sort(unqs.PC[unqs.PC > 0], decreasing = TRUE) 

# Print out how many reads are inferred in the positive control
cat("DADA2 inferred", length(unqs.PC), "sequence variants present in the selected sample.\n")
# DADA2 inferred 342 sequence variants present in the selected sample.

# Plot the number of sequences in the SVs found in the positive control
plot(unqs.PC, ylab = "Number of seqs/SV, Pos Control", xlab = "SVs", main = "Number of sequences/SV in one small intestine sample") 

# OPTIONAL: add taxonomy and see what's in the positive control, for example if you want to know if you only got back what you put in.
# taxa.PC <- assignTaxonomy(unqs.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
# write.csv(taxa.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_PC_EXAMPLE.csv') 






# Lab 5 (Phyloseq)------------------------------
library(phyloseq)

# load seqtab as needed
seqtab.nochim <- readRDS('seqtab_nochim.RDS') 


# load metadata table as needed
meta <- read.csv("Crypto_bacteria_metadata.csv", 
                 header = TRUE, 
                 nrow = 23, # weird trailing empty rows in this csv
                 row.names = 1) 
meta <- rename(meta, "Sample_source" = "Sample_type")
saveRDS(meta, "meta.rds")


# Check the sample sames and see if they match
row.names(meta)
row.names(seqtab.nochim)
identical(row.names(meta), row.names(seqtab.nochim))


## create a phyloseq object with all samples (subset out later)
EX_ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), # even though it's called an OTU table, it will use the SVs from my seqtab
                  sample_data(meta))

EX_ps
saveRDS(EX_ps, "crypto_ex_ps.rds")
# how many samples made it this far? 23 samples, 11 variables, 3724 taxa



# If you have sequenced controls, it's a good idea to look at your data in comparison to them.

# Alpha diversity peek 
plot_richness(EX_ps, x = "Sample_type", # CHANGE the x-axis to a factor of your choice
              measures = c("Observed","Chao1", "Shannon"), # these are some of the alpha diversity measures you can use 
              color = "Sample_source") + # CHANGE the color to a factor of your choice
  theme_bw() + # a ggplot theme to make the graph look nice
  theme(axis.text.x = element_text(angle = 75, hjust = 1))


ggsave("alpha_diversity_plot.png")


## Note: hopefully, diversity/richness is lower in the Negative controls than real samples
## SUMMARY: Generally, negative controls sort of are lower



# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
plot(sort(taxa_sums(EX_ps), TRUE), 
     type = "h", 
     ylim = c(0, 20)) #limit the y-axis to better see the long tail of rare taxa

# optional to save a copy of this figure


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
EX.ord <- ordinate(EX_ps, #calculate similarities
                   method = "PCoA", #ordination type
                   "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)
saveRDS(EX.ord, "crypto_ex_ord.rds")

plot_ordination(EX_ps, EX.ord, type = "samples", color = "Sample_type", shape = "Host_source") +
  ggtitle("Ordination plot")
ggsave("crypto_ordination_plot.png")
# Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues (see next section)
# Horse-shoe patterns indicate underlying patterns to the data
# The negative controls are supposed to be separate from the experimental samples but mine are not

# save this graph for later


# Removing contamination using negative controls (Lab 5)-----------------------------------
### This section can only be accomplished if you have sequenced negative controls and/or DNA quantification data from when the pool was created.  If you have these, CHOOSE ONE method to follow.  If you don't, skip this section and consider including sequencing controls in your next run.

# Dr. Ishaq's method: https://github.com/SueIshaq/Examples-DADA2-Phyloseq
# decontam in R: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

## ------ Dr Ishaq's method -------------
# Dr. Ishaq's method creates vectors out of the SV table data for negative controls, and subtracts those SVs from the sample data.  Depending on the type of negative control, these are removed from the whole data set or from subsets of batches. Remove PCR and sampling materials negative control SVs fully from all samples, and remove extraction kit SVs fully from each dna_extraction_batch, respectively.
phylo <- readRDS("crypto_ex_ps.rds")

# Look at individual controls to inspect for any trends
  
  phylo
  fresh_DMEM <- subset_samples(phylo, Host_source == "fresh_DMEM")
  fresh_DMEM <- prune_taxa(taxa_sums(fresh_DMEM) > 0, fresh_DMEM)
  # fresh_DMEM
  DMEM_with_antibiotics <- subset_samples(phylo, Host_source == "DMEM_with_antibiotics") 
  DMEM_with_antibiotics <- prune_taxa(taxa_sums(DMEM_with_antibiotics) > 0, DMEM_with_antibiotics)
  # DMEM_with_antibiotics
  Negative_control <- subset_samples(phylo, Host_source == "Negative_control") 
  Negative_control <- prune_taxa(taxa_sums(Negative_control) > 0, Negative_control)
  fresh_DMEM
  DMEM_with_antibiotics
  Negative_control

cleaned <- clean_phylo_data(phylo, neg_con = "Negative_control")

# Start with 3724 taxa and 23 samples
# 426 taxa found in 3 negative control samples
# End with 3298 taxa and 20 samples
saveRDS(cleaned, "clean_data_Ishaq.rds")

# 10. Did it work? Check your ordination again
cleaned_ord <- ordinate(cleaned, #calculate similarities
                           method = "PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(cleaned, cleaned_ord, type = "samples", color = "Sample_type", title = "Crypto after cleaning out negative controls")
ggsave("ord_after_clean.png")

# save this graph for later

### Summary of decontamination ----
# CHANGE ME: Which types of negative controls did you use and why? 
# Used all available negative controls. No batch controls necessary -- one extraction batch

# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?
# No obvious clusters before cleaning but maybe a little bit better spread after



## ------ Decontam method -------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam) # This packages identify contaminants by frequency of SVs.

# The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?

# 1. Look for contaminants using DNA quantification data from your metadata file
contamdf.freq <- isContaminant(EX_ps, method="frequency", conc="quant_reading") #CHANGE conc= "quant_reading" to the column heading that holds the concentration information in your metadata

# 2. See what it looks like 
head(contamdf.freq) 

# 3. Make it into a table
table(contamdf.freq$contaminant)

# 4. See which SVs are categorized as contaminants based on frequency
head(which(contamdf.freq$contaminant)) 

# 5. Get rid of the contaminants 
EX_ps.noncontamfreq <- prune_taxa(!contamdf.freq$contaminant, EX_ps)

# 6. How much is left
EX_ps.noncontamfreg
# CHANGE ME: number of samples and SVs left 



# The second contaminant identification method uses prevalence (presence/absence) across samples as compared to the prevalence in negative controls. Essentially, is it in most of your samples but not most of your controls?


# 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
sample_data(EX_ps)$is.neg <- sample_data(EX_ps)$Sample_or_Control == "Control Sample"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(EX_ps, method="prevalence", neg="is.neg", threshold=0.5)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))

# 11. get rid of the contaminants 
EX_ps_NC_decontam_clean <- prune_taxa(!contamdf.prev05$contaminant, EX_ps.noncontamfreq) # remove them from the original phyloseq object, or the one with cleaned out SVs by frequency

# 12. Check how many are left
EX_ps_NC_decontam_clean
# CHANGE ME: number of samples and SVs left


# 13. Did it work? Check your ordination again
EX_cleaner.ord <- ordinate(EX_ps_NC_decontam_clean, #calculate similarities
                           method ="PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(EX_ps_NC_decontam_clean, EX_cleaner.ord, type="samples", color="dna_extraction_batch", title="After cleaning out negative controls")

# save this graph for later

### Summary of decontamination ----
# CHANGE ME: Which types of negative controls did you use and why? 
# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?



# Lab 6_A: DADA2 assigning taxonomy from a reference database file (Lab 6) ---------------------------------


# Be sure to have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 

# The Silva file shown below is used for 16S rRNA (prokaryotes) and nicely formatted versions can be downloaded from the DADA2 website, which also contains some (18S) eukaryotic versions.  Let me know if you need a custom one for your dataset.
#
# 1. take the SVs (otu table) from the phyloseq object and make it into a dataframe again

clean_data <- readRDS("clean_data_Ishaq.rds")

seqtab.nochim.decontam = as(otu_table(clean_data), "matrix")
seqtab.nochim.decontam = as.matrix(seqtab.nochim.decontam)


## 2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops.
clean_all.taxa <- assignTaxonomy(seqtab.nochim.decontam, 
                                 'silva_nr99_v138.1_train_set.fa.gz',        # CHANGE file path
                                 tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                 minBoot = 75, verbose = TRUE) ; beep("treasure") 
# took 7 min for me


saveRDS(clean_all.taxa, 'clean_all_taxa.rds')
write.csv(clean_all.taxa, 'clean_all_taxa.csv')
# clean_all.taxa

# OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.
# clean_all.taxa <- readRDS("clean_all_taxa.rds")


clean_all.taxa.species <- addSpecies(clean_all.taxa, 'silva_species_assignment_v138.1.fa.gz', allowMultiple = FALSE, verbose = FALSE) ; beep("treasure") 
# this took 15 min

saveRDS(clean_all.taxa.species, 'clean_all.taxa.species.rds') 
write.csv(clean_all.taxa.species, 'clean_all.taxa.species.csv') 





## Remake the phyloseq object with the new taxonomy file ----------------------

# reload metadata table as needed
meta <- readRDS("meta.rds")
clean_data <- readRDS("clean_data_Ishaq.rds")

#reload taxa table as needed
clean_all.taxa.species <- readRDS('clean_all.taxa.species.rds')


otu_t <- otu_table(clean_data)

## create a phyloseq object with all samples
phylo_clean_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                     sample_data(meta),
                                     tax_table(clean_all.taxa.species))
saveRDS(phylo_clean_with_species, "phylo_clean_with_species.rds")




# Lab 6_B: Clean out unwanted taxa ------------------------------------
# Regardless of whether you included sequenced negative controls, you can remove taxa which are of no interest to you

# clean out chloroplast and mitochondria. can also elect to remove other contaminating domains or kingdoms as needed. 
#== it means to keep that, or set as equal to it. if you have != it means not that (not matching). & adds more sections to this.
# require(dplyr)

phylo_clean_with_species <- readRDS("phylo_clean_with_species.rds")

# Optional: explore your taxonomy before filtering. Use the tax table you made
df <- as.data.frame(EX_ps) #example
df <- as.data.frame(phylo_clean_with_species) #not working
# I think this is right
df <- as.data.frame(phylo_clean_with_species@tax_table)
table(df$Kingdom)
table(df$Phylum)
table(df$Class)
table(df$Order)
table(df$Family)

phylo_clean_wanted_species  <- phylo_clean_with_species  %>% # CHANGE ME to your phyloseq object name 
  subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast") # CHANGE ME to taxa you want to remove

## clean out taxa/SV columns that are no longer present
phylo_clean_wanted_species <- prune_taxa(taxa_sums(phylo_clean_wanted_species) > 0, phylo_clean_wanted_species)
phylo_clean_wanted_species <- prune_samples(sample_sums(phylo_clean_wanted_species) > 0, phylo_clean_wanted_species)
phylo_clean_wanted_species
# 20 samples and 2459 SVs left 
phylo_clean_with_species
# Started with 20 samples and 3298 taxa

# Save your clean phyloseq object
saveRDS(phylo_clean_wanted_species, 'phylo_clean_wanted_species.rds') 

# Reload as needed
# phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds') 

# Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y

# Reload as needed
# EX_ps_clean <- readRDS('EX_ps_clean_phyloseq_object.RDS') 
# Reload as needed
# phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds')
load_lab_seven_partA()


# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/


tab <- otu_table(phylo_clean_wanted_species)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows

r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) # better
df <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE, tidy = TRUE) # this could be used to re-plot with ggplot
# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_clean_wanted_species)))
# smallest reads (sequences) in a sample ______ Sample 57, 1250 reads
# largest in a sample ______ Sample 51, 87328 reads
# number of samples with <5000 ______ 7/20

# Lamb data showed not much difference between rarification with or without replace
# Just to with replace to simplify, at least for now

  # phylo_no_replace_rar <- rarefy_even_depth(phylo_clean_wanted_species, 
  #                                           sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
  #                                           replace = FALSE, #sampling with or without replacement
  #                                           trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
  #                                           rngseed = 711, 
  #                                           verbose = TRUE)

  # 
  # saveRDS(phylo_no_replace_rar, 'clean_phylo_no_replace_rarified.rds')


phylo_rarified <- rarefy_even_depth(phylo_clean_wanted_species, 
                                       sample.size = 3000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                       replace = TRUE, #sampling with or without replacement
                                       trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                       rngseed = 711, 
                                       verbose = TRUE)

# set.seed(711) was used
# 3 samples removed. Set read threshold at 3k instead of 5k to maintain 17 of the 20 samples
# 1176 OTUs removed

saveRDS(phylo_rarified, 'clean_phylo_replace_rarified.rds')
write.csv(otu_table(phylo_rarified), 'clean_phylo_replace_rarified.csv')



# Alpha diversity (Lab 7 and 8) ---------------------------------
## Alpha diversity graphics (Lab 7) ----------------------


# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
# sample_data(phylo_rarified)$Diet <- factor(sample_data(phylo_replace_rar)$Diet, levels = c("LooseAlfalfa", "PelletedAlfalfa"), labels = c("Loose Alfalfa", "Pelleted Alfalfa")) #CHANGE ME

# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")



sample_data(phylo_rarified)$Treatment <- paste(sample_data(phylo_rarified)$Sample_type, sample_data(phylo_rarified)$Crypto_added, sep = "_")
sample_data(phylo_rarified)$Sample_type <- factor(sample_data(phylo_rarified)$Sample_type, 
                                                  levels = c("large_intestine_biopsy_culture", "small_intestine_biopsy_culture"), 
                                                  labels = c("Large intestine", "Small intestine"))
sample_data(phylo_rarified)$Treatment2 <- paste(sample_data(phylo_rarified)$Treatment, sample_data(phylo_rarified)$Antibiotics, sep = "_")
sample_data(phylo_rarified)$Antibiotics <- factor(sample_data(phylo_rarified)$Antibiotics, 
                                                  levels = c("no", "yes"), 
                                                  labels = c("No", "Yes"))

plot_richness(phylo_rarified, 
              x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  geom_violin(trim = TRUE, aes(fill = Crypto_added)) +
  scale_fill_discrete(name = "Crypto added") + 
  geom_boxplot(width = 0.1, aes(group = Treatment, color = Crypto_added)) + 
  scale_color_hue(guide = "none") +
  # facet_grid(.~Sample_type, switch = "x", space = "free", scales = "free") +
  facet_grid(rows = vars(Sample_type), cols = vars(Antibiotics), space = "free", switch = "x") +
  ylab("Observed Bacterial Richness (SVs)") +
  xlab("Antibiotics added") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.text.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),) +
  labs(title = "Observed SV diversity",
       subtitle = "Crypto data<br>Rarefied")
ggsave("crypto_diversity_plot_offset.png")




# two different richness plots grouped

plot1 <- plot_richness(phylo_rarified, 
                       x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) +
  geom_violin(trim = TRUE, aes(fill = Crypto_added)) + 
  scale_fill_discrete(name = "Crypto added") +
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group box plots
  facet_grid(.~Sample_type, switch = "x", space = "free", scales = "free") +
  theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)")

plot1

plot2 <- plot_richness(phylo_rarified, 
                       x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) + #example, get rid of x axis labels
  geom_violin(trim = TRUE, aes(fill = Crypto_added)) + #optional. CHANGE ME, A is whatever factor to color violins
  scale_fill_discrete(name = "Crypto added") +
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group box plots
  facet_grid(.~Sample_type, switch = "x", space = "free", scales = "free") +
  theme(legend.position = "none") + #use to get rid of your legend
  ylab("Shannon Bacterial Richness (SVs)")

# plot together
p <- ggarrange(plot1,
          ggarrange(plot2, nrow = 1, labels = c("B")),
          ncol = 2, labels = "A", common.legend = TRUE, legend = "bottom")
annotate_figure(p, top = text_grob("Alpha diversity", face = "bold", size = 15))

ggsave("different_diversity_plots_panels.png")

    

# use phyloseq to measure alpha diversity
phylo_clean.rar.rich <- estimate_richness(phylo_rarified, measure = c("Observed", "Shannon")) #change to whatever measures you want

# OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# install.packages("remotes")
# remotes::install_github("twbattaglia/btools")
# phylo_faith <- estimate_pd(clean_phylo_replace_rarified)

# OPTIONAL measure evenness for each SV individually
phylo_clean.rar.even_SV <- evenness(otu_table(phylo_rarified))

# measure evenness for each sample
phylo_clean.rar.even <- phylo_clean.rar.rich$Shannon/log(phylo_clean.rar.rich$Observed)


# Coerce to data.frame and add the metadata for these samples
phylo_clean.rar.sd = as(sample_data(phylo_rarified), "matrix")
phylo_clean.rar.sd = as.data.frame(phylo_clean.rar.sd)
phylo_clean.rar.rich.df <- cbind(phylo_clean.rar.rich, phylo_clean.rar.even, phylo_clean.rar.sd)
saveRDS(phylo_clean.rar.rich.df, "phylo_clean.rar.rich.df.rds")

# make a histogram to look at the shape of the data (bell curve? skew?). You can save this graph for your own benefit if you want.
hist(phylo_clean.rar.rich.df$Observed)
hist(phylo_clean.rar.rich.df$Observed, breaks = 14)

# Want to know how much skew there is? Measure Kurtosis (a.k.a. tailedness). You can save this graph for your own benefit if you want.

kurtosis(phylo_clean.rar.rich.df$Observed)
# mine is 0.2228839
# Red flag is as it approaches +/- 2. High value means a lot of positive outliers


head(phylo_clean.rar.rich.df)

#check the distribution of your data, which will change your approach
shapiro.test(phylo_clean.rar.rich.df$Shannon)
# W = 0.93056, p-value = 0.2224: my shannon diversity metric is #FIXME distributed (p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution)

shapiro.test(phylo_clean.rar.rich.df$Observed)
# W = 0.96498, p-value = 0.7261, my observed data are #FIXME ?not normal. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution

shapiro.test(phylo_clean.rar.rich.df$phylo_clean.rar.even)
# W = 0.96669, p-value = 0.7583, #FIXME not normally distributed. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution


# DONE UP TO HERE ------------------------------------