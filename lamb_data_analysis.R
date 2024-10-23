## Intro and package loading ---------------------------

setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/load_data.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/AVS554_packages.R")
# install_necessary()
# install_optional()

load_libraries()


  ## Load raw data (Lab 2) ---------------------------

# When you first get a new dataset, you should play around and find out how many sequences, samples, and sequence lengths you have, as well as look at the quality score of the dataset. This will help you plan your quality control steps. 
# This only needs to be done once, after which you will use the Filtered Data you create in Lab 1.
# Note: If you are using data from Roche 454 pyrosequencing you may want to use Pyronoise or "shh.flows" in MOTHUR to denoise your raw data and create cleaner fastq files.  This is not needed for Illumina data.



# 1. Specify the folder path to the raw (unfiltered) sequence files by assigning the folder path to a variable.  Use fastq.qz (preferred) or fastq formats.  For organization, I keep forwards (R1) and reverses (R2) in separate folders. If you have  "Undetermined" sequence files, remove this from the folder - these are sequences that your sequencer could not identify the barcode and assign it to a sample.
setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")


file_info <- metadata_import(types = c("Forward", "Reverse")) # saves a copy as an rds

# Look at the whole dataset to aggregate all samples onto one graph. 
# One plot for each data type will be in the wd in .png format
qualplots(file_info); beep() # Default data type is raw
# qualplots(subset(file_info, data_subset == "Forward"))

# run this command to create a quality control report
for (i in 1:nrow(file_info)) {
  rqc(path = file_info$raw_data_dir[i][1], 
      sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
      file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())
}; beep()


# Optional: search for strings in your data
search_strings("AAAAAAAAAAA", file_info$file_names_full[[1]][[9]]) # pick a random one to search, as an example


DNA_to_read <- readDNAStringSet(file_info$file_names_full[[1]][1], format = "fastq", with.qualities = TRUE)
# Pull out a sequence I can BLAST
toString(DNA_to_read[1])
# 16S

# saveRDS(file_info, "file_info.rds") # this happens inside the "metadata import" function
file_info <- readRDS("file_info.rds")

### Summary of quality ----
# This dataset contains 34 number of samples, and 2560717 forward reads which are 300 bases long
# The F reads quality goes below 35 at cycle ~250.
# The R reads quality goes below 35 at cycle 20. The quality not is good enough to use


## Filter and trim raw sequences (Lab 2) ---------------------------
# Filer and trim sequences based on various parameters.  There are two ways to do this, on the Forward, R1 reads only (they tend to be better quality than Reverse, R2s), or you can do it on Forward/Reverse at the same time with the intent to merge them together into contigs later.

subset(file_info, data_subset == "Forward") %>%
  filtering(trunclen = 250) # ; beep("treasure")
# I retained 1866909 reads

qualplots(file_info, type = "filtered"); beep("treasure")

filtoutput <- readRDS("lamb_data_Forward_filtered_output.rds")

# check the dimensions of the variable you created, outputs two numbers: rows (# samples), columns (# of info it added)
dim(filtoutput) 
# take a look at the counts
View(filtoutput)
# get a sum total for raw reads in and filtered reads out
colSums(filtoutput)
summary(as.data.frame(filtoutput)$reads.in)
summary(as.data.frame(filtoutput)$reads.out)

### Summary of filtering ----
# Which parameters did you use and why? 
    # I threw out reverse reads due to poor quality. I trimmed at 250, where the quality was starting to dip.
# What was the largest and small number of raw reads?
    # Largest: 316288
    # Smallest: 4058
# What was the largest and small number of filtered reads?
    # Negative controls were generally tossed out, with 3/4 negatives ending up with 0 reads
    # Largest: 270381
# How many total raw reads?
    # 2560717 raw reads (forward)
# How many total filtered reads?
    # I retained 1866909 reads



## Lab 3---------------------------


source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")
load_lab_three() # reads file_info, errF into global environment

# file_info <- readRDS("file_info.rds")
create_and_plot_error_profile(file_info, bases = 1e6) ; beep("treasure") # writes the error profile RDS files
# 4473120 total bases in 18638 reads from 1 samples will be used for learning the error rates.
# This data set is only analyzing forward
errF <- readRDS("error_profile_Forward.RDS")


# to just plot
  # plotErrors(err, nominalQ = TRUE) +
  #   ggtitle("")
  # ggsave(".png")


dadaFs <- dada(list.files(file_info$filt_dir[1], full.names = TRUE),
               err = errF, 
               multithread = TRUE, verbose = FALSE) ; beep("treasure")

seqtab <- makeSequenceTable(dadaFs); beep("treasure")

dim(seqtab) 
# 34 samples  and 12425 of SVs

saveRDS(seqtab, 'seqtab.rds') 

SVs_found_by_sample <- make_SV_summary(seqtab)
View(SVs_found_by_sample)
# Fewest SVs: PelletedAlfalfa_wk0_21_for_filt.fastq.gz -- 30
# Most SVs: PelletedAlfalfa_wk2_93_for_filt.fastq.gz -- 3386


# Lab 4_A: DADA2 Remove chimeras from seqtab ----------------------------	

# Chimeras are accidentally created sequences during lab protocols. Remove them.	


# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")
# seqtab <- readRDS("seqtab.rds")
load_lab_four_partA() # reads seqtab into global environment

# 1. Remove chimeras. Leave the method as consensus. multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE) # a few moments to run
# My samples have 1128 bimeras out of 12425 total sequences

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab.nochim) 
# 34 samples and 11297 SVs, down from 12425

# Calculate the percentage of chimeras identified out of the total
round(sum(seqtab.nochim)/sum(seqtab), digits = 3) # 0.955 retained
round(1 - sum(seqtab.nochim)/sum(seqtab), digits = 3) # meaning 0.045 thrown out

# 3. Save your chimera-free sequence table.
saveRDS(seqtab.nochim, 'seqtab.nochim.rds')


## Lab 4_B: Workflow verification steps  --------------------------
# This section can be optional, but as you are learning, or playing with new data, it is helpful to run some internal checks to assess whether you like how your quality control steps worked.
load_lab_four_partB() # seqtab.nochim, filtoutput


    # 1. Track reads through the analysis here to see how many were lost at each QC step, in case you were too stringent. Only need to count F reads, even if you did paired F/R
    #set function to count the unique reads
    getN <- function(x) {
      sum(getUniques(x)) 
    }


# 2. Load the metadata file that goes with your sequencing data so you can match factors to seq data. This uses a datatable made in Excel and saved as a .csv file
meta <- read.csv("Particle_Size_metadata.csv", 
                 header = TRUE, # specifies that the top row is column names
                 row.names = 1) # row.names specifies which column to call the samples by (sample name), doesn't have to be the first column


# 2.5 check the dimensions of the three data files you need for this to make sure the number of rows matches in each. If they do not, you may need to add/remove rows from your metadata file in case samples were removed/retained from your dataset.
dim(filtoutput) # 34 x 2
dim(seqtab.nochim) # 34 x 30
dim(meta) # 34 x 30


head(row.names(filtoutput))
head(row.names(seqtab.nochim))
head(row.names(meta))

# Need to strip extra labeling

row.names(filtoutput) <- str_remove_all(row.names(filtoutput), "_for.fastq.gz")
row.names(seqtab.nochim) <- str_remove_all(row.names(seqtab.nochim), "_for_filt.fastq.gz")
identical(row.names(filtoutput), row.names(seqtab.nochim)) # true
identical(row.names(filtoutput), row.names(meta)) # false
meta <- meta[order(row.names(meta)),]
identical(row.names(filtoutput), row.names(meta)) #true now

saveRDS(meta, "meta.rds")
saveRDS(filtoutput, "filtoutput.rds")
saveRDS(seqtab.nochim, "seqtab.nochim.rds")

# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab.nochim), meta$Sample_type) #Sample_type should be the column name in your metadata file of the groups in your data

# 4. Assign column names to that new variable
colnames(track) <- c("reads.in", "filtered", "nonchimeras", "Sample_type") # CHANGE "Sample_type" to be the column name in your metadata file of the groups in your data

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'tracked_seqs.RDS') #CHANGE file name to whatever you want

# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))


# plot with Sample_type along the X axis
ggplot(plotData,aes(x = Sample_type, y = as.numeric(totals))) + geom_jitter(aes(color = type)) + 
  ylab("Sequences") + 
  xlab("Sample_type") +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1)) +
  ggtitle("Reads by sample type")
ggsave("reads_filt_chim_sample_type.png")

# or, plot with QA stage along the X axis
ggplot(plotData,aes(x = type,y = as.numeric(totals))) + geom_jitter(aes(color = Sample_type)) + 
  ylab("Sequences") + 
  xlab("QA stage") +
  ggtitle("Reads by filtering step")
ggsave("reads_sample_type_QC_status.png")








# 9. Randomly select a negative control sample to see what's in it
unqs.NC <- seqtab.nochim["Mock_S192",]

# Sort by # seqs and drop SVs absent in the negative control sample
unqs.NC <- sort(unqs.NC[unqs.NC > 0], decreasing = TRUE) 

# Print out how many reads are inferred in the negative control
cat("DADA2 inferred", length(unqs.NC), "sequence variants present in the selected sample.\n")
# DADA2 inferred 148 sequence variants present in the selected sample.

# Plot the number of sequences in the SVs found in the negative control sample.
plot(unqs.NC, ylab = "Number of seqs/SV, Neg Control", xlab = "SVs", main = "Number of sequences/SV in Negative Control - Mock") #CHANGE ME

# OPTIONAL: add taxonomy and see what's in the negative control, for example if you want to know if you have the same contaminant showing up in negative controls over time.
# taxa.NC <- assignTaxonomy(unqs.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
# write.csv(taxa.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_NC_EXAMPLE.csv') 



# 10. Randomly select a positive control to see what's in it
unqs.PC <- seqtab.nochim["LooseAlfalfa_wk0_13",] #CHANGE ME, update with sample name of a positive control, such as Mock_S192

# Sort by # seqs and drop SVs absent in the positive control sample
unqs.PC <- sort(unqs.PC[unqs.PC > 0], decreasing = TRUE) 

# Print out how many reads are inferred in the positive control
cat("DADA2 inferred", length(unqs.PC), "sequence variants present in the selected sample.\n")
# DADA2 inferred 342 sequence variants present in the selected sample.

# Plot the number of sequences in the SVs found in the positive control
plot(unqs.PC, ylab = "Number of seqs/SV, Pos Control", xlab = "SVs", main = "Number of sequences/SV in one experimental sample") 

# OPTIONAL: add taxonomy and see what's in the positive control, for example if you want to know if you only got back what you put in.
# taxa.PC <- assignTaxonomy(unqs.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
# write.csv(taxa.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_PC_EXAMPLE.csv') 



# Phyloseq First look (Lab 5_A) ------------------------------
library(phyloseq)
setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")
load_lab_five_partA() # seqtab.nochim, meta

  # # load seqtab as needed
  # seqtab.nochim <- readRDS('seqtab_nochim.RDS') 
  # meta <- readRDS("meta.rds")


# Check the sample sames and see if they match
identical(row.names(meta), row.names(seqtab.nochim))


## create a phyloseq object with all samples (subset out later)
phylo <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), # even though it's called an OTU table, it will use the SVs from my seqtab
                  sample_data(meta))

phylo
saveRDS(phylo, "phylo.rds")
# how many samples made it this far? 34 samples, 11297 taxa. Taxa are the SV columns



# If you have sequenced controls, it's a good idea to look at your data in comparison to them.

# Alpha diversity peek 
plot_richness(phylo, x = "Sample_type", # CHANGE the x-axis to a factor of your choice
              measures = c("Observed","Chao1", "Shannon"), # these are some of the alpha diversity measures you can use 
              color = "Treatment") + # CHANGE the color to a factor of your choice
  theme_bw() + # a ggplot theme to make the graph look nice
  theme(axis.text.x = element_text(angle = 75, hjust = 1))
# Ignore the warning that prints -- it would be relevant for OTUs instead of SVs
ggsave("alpha_diversity_plot.png")


## Note: hopefully, diversity/richness is lower in the Negative controls than real samples
## SUMMARY: Generally yes, negative controls are lower


# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
plot(sort(taxa_sums(phylo), TRUE), 
     type = "h", 
     ylim = c(0, 20)) #limit the y-axis to better see the long tail of rare taxa


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
phylo_ord <- ordinate(phylo, #calculate similarities
                   method = "PCoA", #ordination type
                   "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)


plot_ordination(phylo, phylo_ord, type = "samples", shape = "Sample_type", color = "DNA_extraction_batch") +
  ggtitle("Lamb, ordination plot, pre-clean")
ggsave("lamb_before_ordination_plot.png") # save this graph for later
plot_ordination(phylo, phylo_ord, type = "samples", color = "Sample_type", shape = "DNA_extraction_batch") +
  ggtitle("Lamb, ordination plot, pre-clean")
ggsave("lamb_before_ordination_plot2.png")
# Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues (see next section)
# Horse-shoe patterns indicate underlying patterns to the data
# The negative controls are supposed to be separate from the experimental samples and they are




# Lab 5_B: Removing contamination using negative controls-----------------------------------
### This section can only be accomplished if you have sequenced negative controls and/or DNA quantification data from when the pool was created.  If you have these, CHOOSE ONE method to follow.  If you don't, skip this section and consider including sequencing controls in your next run.

setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/lamb_data")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
# phylo <- readRDS("phylo.rds") # read in phyloseq dataset if necessary
# meta <- readRDS("meta.rds")
load_lab_five_partB() # phylo and meta


# Dr. Ishaq's method: https://github.com/SueIshaq/Examples-DADA2-Phyloseq

## ------ Dr Ishaq's method -------------

# Dr. Ishaq's method creates vectors out of the SV table data for negative controls, and subtracts those SVs from the sample data.  Depending on the type of negative control, these are removed from the whole data set or from subsets of batches. Remove PCR and sampling materials negative control SVs fully from all samples, and remove extraction kit SVs fully from each dna_extraction_batch, respectively.
# With modifications by Emily

clean_data <-  clean_phylo_data(phylo, neg_con = "NegCon_PCR", batch_con = "NegCon_kit")
clean_data
# Ends up 11162 and 30
saveRDS(clean_data, "clean_data_Ishaq.rds")

# Did it work? Check your ordination again
clean_ord <- ordinate(clean_data, #calculate similarities
                           method = "PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(clean_data, clean_ord, 
                type = "samples", color = "Diet", 
                shape = "DNA_extraction_batch", title = "Lamb, after clean")
ggsave("lamb_after_NC_clean_ordination_plot.png") # save this graph for later



### Summary of decontamination ----
# Which types of negative controls did you use and why? 
# I used the NegCon_PCR as a control to remove any contaminants introduced at the PCR step
# I used NegCon_kit as a control to remove any contaminants introduced at the extraction step
# I left the mock sample in as a pipeline control
# The mock was successfully removed by the data cleaning using the other negative controls!

# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?
# The clustering was by sample type (negative controls clustered separately), not extraction batch
# After clean, no distinct clusters



## ------ Decontam method -------------
# decontam in R: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html




# The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?

# 1. Look for contaminants using DNA quantification data from your metadata file
contamdf.freq <- isContaminant(phylo, method = "frequency", conc = "DNA_extract_concentration") #CHANGE conc= "quant_reading" to the column heading that holds the concentration information in your metadata
saveRDS(contamdf.freq, "contamdf_freq.rds")


# 2. See what it looks like 
head(contamdf.freq) 

# 3. Make it into a table
table(contamdf.freq$contaminant)

# 4. See which SVs are categorized as contaminants based on frequency
head(which(contamdf.freq$contaminant)) 

# 5. Get rid of the contaminants 
phylo.noncontamfreq <- prune_taxa(!contamdf.freq$contaminant, phylo)

# 6. How much is left
phylo.noncontamfreq
# 11248 taxa, 34 samples. Started with 11297



# The second contaminant identification method uses prevalence (presence/absence) across samples as compared to the prevalence in negative controls. Essentially, is it in most of your samples but not most of your controls?


# 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
sample_data(phylo)$is.neg <- sample_data(phylo)$Sample_type == "NegCon_kit"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(phylo, method = "prevalence", neg = "is.neg", threshold = 0.5)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))

## Section copied from bioconductor vignettes and modified for my data
contamdf.prev <- isContaminant(phylo, method = "prevalence", neg = "is.neg")
table(contamdf.prev$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phylo, function(abund) 1*(abund > 0))
ps.pa
# 11297 taxa and 34 samples, so the same as start
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_type == "NegCon_kit", ps.pa)
ps.pa.neg
# 11297 taxa, 2 samples

ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_type == "Experimental", ps.pa)
ps.pa.pos
# 11297 taxa, 30 samples
unique(ps.pa.pos@sam_data[["Sample_type"]])
# took out everything except Experimental

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos = taxa_sums(ps.pa.pos), pa.neg = axa_sums(ps.pa.neg),
                    contaminant = contamdf.prev$contaminant)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True Samples)") +
  ggtitle("Prevalence of contaminants, by decontam method")
ggsave("prevalence_contams_by_decontam.png")



# 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
sample_data(phylo)$is.neg <- sample_data(phylo)$Sample_type == "NegCon_kit"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(phylo, method = "prevalence", neg = "is.neg", threshold = 0.5)
names(contamdf.prev05)
head(contamdf.prev05)
dim(contamdf.prev05)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))
# 11. get rid of the contaminants 
# FIXME
# prune_taxa(!contamdf.prev05$contaminant, phylo.noncontamfreq) # remove them from the original phyloseq object, or the one with cleaned out SVs by frequency
# head(!contamdf.prev05$contaminant)


# phylo_NC_decontam_clean <- prune_taxa(names(!contamdf.prev05$contaminant), phylo.noncontamfreq) # remove them from the original phyloseq object, or the one with cleaned out SVs by frequency
# 12. Check how many are left
# phylo_NC_decontam_clean
# 11248 taxa, 34 samples
# but this plot looks bad


ESK_phylo_NC_decontam_clean <- prune_taxa(row.names(subset(contamdf.prev05, contaminant == FALSE)), phylo.noncontamfreq)
ESK_phylo_NC_decontam_clean
# 11180 taxa, 34 samples
# 13. Did it work? Check your ordination again
phylo_decontam_cleaner.ord <- ordinate(ESK_phylo_NC_decontam_clean, #calculate similarities
                           method = "PCoA", #ordination type
                           "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(ESK_phylo_NC_decontam_clean, phylo_decontam_cleaner.ord, 
                type = "samples", 
                color = "DNA_extraction_batch", 
                shape = "Sample_type",
                title = "Lamb, after cleaning out negative controls by decon")
ggsave("lamb_decontam_after.png")
# looks ok I think?

plot_ordination(ESK_phylo_NC_decontam_clean, phylo_decontam_cleaner.ord, 
                type = "samples", 
                color = "Treatment", 
                shape = "Sample_type",
                title = "Lamb, after cleaning out negative controls by decon")
ggsave("lamb_decontam_after2.png")


### Summary of decontamination ----
# FIXME
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
                           minBoot = 75, verbose = TRUE) # ; beep("treasure") 



saveRDS(clean_all.taxa, 'clean_all_taxa.rds')
write.csv(clean_all.taxa, 'clean_all_taxa.csv')
# clean_all.taxa

# OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.
clean_all.taxa <- readRDS("clean_all_taxa.rds")

# first_half <- clean_all.taxa[1:5579,]
# second_half <- clean_all.taxa[5580:nrow(clean_all.taxa),]
first_half <- clean_all.taxa[1:round(nrow(clean_all.taxa)/2, digits = 0),]
second_half <- clean_all.taxa[round((nrow(clean_all.taxa)/2) + 1, digits = 0):nrow(clean_all.taxa),]

first_half.species <- addSpecies(first_half, 'silva_species_assignment_v138.1.fa.gz', allowMultiple = FALSE, verbose = FALSE)# ; beep("treasure") 
# Takes about 45 minutes
second_half.species <- addSpecies(second_half, 'silva_species_assignment_v138.1.fa.gz', allowMultiple = FALSE, verbose = FALSE) # ; beep("treasure") 
# Also takes about 45 minutes

clean_all.taxa.species <- rbind(first_half.species, second_half.species) # ; beep("treasure")

saveRDS(clean_all.taxa.species, 'clean_all.taxa.species.rds') 
write.csv(clean_all.taxa.species, 'clean_all.taxa.species.csv') 





## Remake the phyloseq object with the new taxonomy file ----------------------

# reload metadata table as needed
meta <- readRDS("meta.rds")
clean_data <- readRDS("clean_data_Ishaq.rds")

#reload taxa table as needed
clean_all.taxa.species <- readRDS('clean_all.taxa.species.rds')

class(otu_table(clean_data))
class(meta)
class(tax_table(clean_all.taxa.species))
class(clean_all.taxa.species)

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
# 30 samples and 7820 SVs left 
phylo_clean_with_species
# Started with 30 samples and 11158 taxa

# Save your clean phyloseq object
saveRDS(phylo_clean_wanted_species, 'phylo_clean_wanted_species.rds') 

# Reload as needed
  # phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds') 




# Write out a description of experimental design (Homework)

# Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y

# Reload as needed
# EX_ps_clean <- readRDS('EX_ps_clean_phyloseq_object.RDS') 
# Reload as needed
# phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds')
load_lab_seven_partA()


# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
# From Dr. Ishaq's example:
# EX-rarec <- rarecurve(otu_table(EX_ps_clean), step = 10, cex=0.5, label = TRUE) # step is major tick marks on x axis (x 1000), cex is text size, label is sample name
# rarec <- rarecurve(otu_table(phylo_clean_wanted_species), step = 10, cex = 0.5, label = TRUE) # step is major tick marks on x axis (x 1000), cex is text size, label is sample name
# # didn't work
# rarec <- rarecurve(otu_table(phylo_clean_wanted_species, taxa_are_rows = FALSE), step = 10, cex = 0.5, label = TRUE) # step is major tick marks on x axis (x 1000), cex is text size, label is sample name
# # still doesn't work

# this is from stack overflow
tab <- otu_table(phylo_clean_wanted_species)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
# rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)

r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = TRUE) 
r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) # better?
df <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE, tidy = TRUE) # idk what to do with this really

# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_clean_wanted_species)))
# smallest reads (sequences) in a sample ______ PelletedAlfalfa_wk0_29 3696
# largest in a sample ______ PelletedAlfalfa_wk2_93
# number of samples with >5000 ______ 26/30

phylo_no_replace_rar <- rarefy_even_depth(phylo_clean_wanted_species, 
                                     sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                     replace = FALSE, #sampling with or without replacement
                                     trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                     rngseed = 711, 
                                     verbose = TRUE)
# set.seed(711) was used
# 4 samples removed
# 1579 OTUs removed

saveRDS(phylo_no_replace_rar, 'clean_phylo_no_replace_rarified.rds')

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(phylo_no_replace_rar), 'clean_phylo_no_replace_rarified.csv')

phylo_replace_rar <- rarefy_even_depth(phylo_clean_wanted_species, 
                                          sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                          replace = TRUE, #sampling with or without replacement
                                          trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                          rngseed = 711, 
                                          verbose = TRUE)

# set.seed(711) was used
# 4 samples removed
# 1630 OTUs removed

saveRDS(phylo_replace_rar, 'clean_phylo_replace_rarified.rds')
write.csv(otu_table(phylo_replace_rar), 'clean_phylo_replace_rarified.csv')


# Pick up the statistical analysis here (Lab 7+) -------------------------------------
# Once you are satisfied with your QC, you can start from here with your phyloseq object as needed


load_lab_seven_partB() # Loads meta, phylo_clean_wanted_species, phylo_no_replace_rar, and phylo_replace_rar


# phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds')
# phylo_no_replace_rar <- readRDS('clean_phylo_no_replace_rarified.rds')

# # Rarefied version
# EX_ps_clean.rar <- readRDS('EX_ps_clean_rarefied_phyloseq_object.RDS')





# Alpha diversity (Lab 7 and 8) ---------------------------------
## Alpha diversity graphics (Lab 7) ----------------------
load_lab_seven_partB() # Phylo_clean_wanted_species; meta; phylo_no_replace_rar, phylo_replace_rar
# FIXME idk why this isn't working

phylo_clean_wanted_species <- readRDS('phylo_clean_wanted_species.rds')
meta <- readRDS("meta.rds")
# phylo_no_replace_rar <<- readRDS('clean_phylo_no_replace_rarified.rds')
phylo_replace_rar <- readRDS("clean_phylo_replace_rarified.rds")

# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
sample_data(phylo_replace_rar)$Diet <- factor(sample_data(phylo_replace_rar)$Diet, levels = c("LooseAlfalfa", "PelletedAlfalfa"), labels = c("Loose Alfalfa", "Pelleted Alfalfa")) #CHANGE ME
saveRDS(phylo_replace_rar, "clean_phylo_replace_rarified.rds")


# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

plot_richness(phylo_replace_rar, 
              x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) + #example, get rid of x axis labels
  geom_violin(trim = TRUE, aes(fill = Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group box plots
  facet_grid(.~Week, switch = "x") + # facet_grid(.~Week, space = "free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  xlab("Week") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  labs(title = "SVs compared to diet, per week",
       subtitle = "Rarefied<br>With replacement")
ggsave("richness_plot_rarefied_w_replacement.png")

plot_richness(phylo_no_replace_rar, 
              x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) + #example, get rid of x axis labels
  geom_violin(trim = TRUE, aes(fill = Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group box plots
  facet_grid(.~Week, switch = "x", space = "free", scales = "free") + # facet_grid(.~Week, space = "free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  xlab("Week") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  labs(title = "SVs compared to diet, per week",
       subtitle = "Rarefied<br>Without replacement")
ggsave("richness_plot_rarefied_no_replacement.png")
# Difference is so minor, just go with the rareified + replacement data


# two different richness plots grouped

plot1 <- plot_richness(phylo_replace_rar, 
                       x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim = TRUE, aes(fill = Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  ylab("Observed Bacterial Richness (SVs)")
plot1

plot2 <- plot_richness(phylo_replace_rar, 
                       x = "Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim = TRUE, aes(fill = Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
  ylab("Bacterial Richness (SVs)")


# plot together
ggarrange(plot1,
          ggarrange(plot2, nrow = 1, labels = c("B")),
          ncol = 2, labels = "A", common.legend = TRUE)
ggsave("different_richness_plots_panels.png")


# DONE UP TO HERE ------------------------------------
# Skipped some stuff ------------------------------------

# richness plot with significance added



# richness plot observed SVs with lines to fit the view screen
plot_richness(EX_ps_clean.rar, 
              x="Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures="Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") + 
  ylim(0,1500) + #define the y axis min/max
  geom_segment(aes(x = 1, y = 1200, xend = 2, yend = 1200)) +  geom_text(x = 1.5, y = 1250, label = "***") # add a drawn in line and significance tag, adjusting the x and y coordinates till it fits where you want in the window.  Add another for each line to add. As written, this will fit the view window you have, if you adjust that your segments will not adjust with it.




## Alpha diversity plotted against other metadata -----

# use phyloseq to measure alpha diversity
phylo_clean.rare.rich <- estimate_richness(clean_phylo_replace_rarified, measure = c("Observed", "Shannon")) #change to whatever measures you want
  
  # OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
  # install.packages("remotes")
  # remotes::install_github("twbattaglia/btools")
  phylo_faith <- estimate_pd(clean_phylo_replace_rarified)
  
  # measure evenness for each sample
  EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)
  
  # Coerce to data.frame and add the metadata for these samples
  EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
  EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
  EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)
  
  dim(EX_ps_clean.rar.rich.df)
  
  # make a graph using that dataframe. this is a generic example, you will want to personalize this.
  ggplot(data=EX_ps_clean.rar.rich.df, aes(x=Diet, y=Observed)) + 
    theme_minimal() + 
    geom_point(aes(color=Treatment), size = 3) + # sets colors by group and a set point size 
    xlab("") + 
    ylab("Bacterial Richness (SVs)") + 
    theme(text = element_text(size = 20)) # increases font size
  
  
  # make that same graph but drop any samples that lack data for that FactorA (replace FactorA in the code with your factor name).
  ggplot(data=subset(EX_ps_clean.rar.rich.df, !is.na(FactorA)), aes(x=Temperature, y=Observed)) + 
    theme_minimal() + 
    geom_point(aes(color=Group), size = 3) + # sets colors by group and a set point size 
    xlab("Temperature of Ocean Water (C)") + 
    ylab("Bacterial Richness (SVs)") + 
    theme(text = element_text(size = 20)) # increases font size
  
  
  
  # Make an alpha diversity table: https://www.geeksforgeeks.org/create-table-from-dataframe-in-r/
  EX_table = as.table(table(EX_ps_clean.rar.rich.df$Diet, EX_ps_clean.rar.rich.df$Week)) 
  EX_table



# Lab 8: Alpha diversity metrics statistics--------------
# Phyloseq can measure and visualize alpha diversity: https://joey711.github.io/phyloseq/plot_richness-examples.html
# phyloseq doesn't do stats or more complex graphs
#FIXME
  # add a load data
  # clean_phylo_replace_rarified <- readRDS("clean_phylo_replace_rarified.rds")
  # Picked up here ------------------------------------

# use phyloseq to measure alpha diversity
phylo_clean.rar.rich <- estimate_richness(clean_phylo_replace_rarified, measure = c("Observed", "Shannon")) #change to whatever measures you want
saveRDS(phylo_clean.rar.rich, "phylo_clean.rar.rich.rds")
  
# use phyloseq to calculate Faith's Diversity metric (optional), https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# EX_faith <- estimate_pd(EX_ps_clean.rar)

# OPTIONAL measure evenness for each SV individually
    # library(asbio)
    # library(microbiome)
    # EX_ps_clean.rar.even_SV <- evenness(otu_table(EX_ps_clean.rar))

# measure evenness for each sample
phylo_clean.rar.even <- phylo_clean.rar.rich$Shannon/log(phylo_clean.rar.rich$Observed)


# Coerce to data.frame and add the metadata for these samples
phylo_clean.rar.sd = as(sample_data(clean_phylo_replace_rarified), "matrix")
phylo_clean.rar.sd = as.data.frame(phylo_clean.rar.sd)
phylo_clean.rar.rich.df <- cbind(phylo_clean.rar.rich, phylo_clean.rar.even, phylo_clean.rar.sd)
saveRDS(phylo_clean.rar.rich.df, "phylo_clean.rar.rich.df.rds")

# make a histogram to look at the shape of the data (bell curve? skew?). You can save this graph for your own benefit if you want.
hist(phylo_clean.rar.rich.df$Observed)
hist(phylo_clean.rar.rich.df$Observed, breaks = 12)

# Want to know how much skew there is? Measure Kurtosis (a.k.a. tailedness). You can save this graph for your own benefit if you want.
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
kurtosis(phylo_clean.rar.rich.df$Observed)
# mine is -0.1604563
# Red flag is as it approaches +/- 2. High value means a lot of positive outliers


head(phylo_clean.rar.rich.df)

#check the distribution of your data, which will change your approach
shapiro.test(phylo_clean.rar.rich.df$Shannon)
# W = 0.72199, p-value = 1.076e-05: my shannon diversity metric is #FIXME distributed (p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution)

shapiro.test(phylo_clean.rar.rich.df$Observed)
#W = 0.95095, p-value = 0.244, my observed data are #FIXME ?not normal. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution

shapiro.test(phylo_clean.rar.rich.df$phylo_clean.rar.even)
# W = 0.60161, p-value = 3.086e-07, #FIXME not normally distributed. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution





## If your alpha diversity is normally distributed (parametric)-----------------


# pick your anova model based on your study design: https://www.statmethods.net/stats/anova.html. Replace "y" with the alpha diversity measure you want, such as observed, evvenness, shannon, etc.
# Replace "A and B" with the factors you want from your own data, and and replace x with numerical data

### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html



# you might run each factor separately to find out if they are significant on their own. 
summary(lm(Observed ~ as.numeric(PHOS), data = EX_ps_clean.rar.rich.df)) # PHOS is the phosphorus values from the lamb rumen data

# if you have multiple groups, this will give you the pairwise comparisons
emmeans(lm,pairwise ~ Diet) 
# paste the output here

# what if you want to add log transformation to a numerical variable in your metadata? Check out the syntax here: https://rpubs.com/marvinlemos/log-transformation


### Wilcoxon (a.k.a. Mann-Whitney) for pairwise numeric or categorical factors -------
# Wilcoxon in R (a.k.a. Mann-Whitney): https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
# Two data samples are independent if they come from distinct populations and the samples do not affect each other. Using the Mann-Whitney-Wilcoxon Test, we can decide whether the population distributions are identical while assuming them to follow the normal distribution. We can use this on non-normal data if we have many more than 30 samples per group and large variance between groups. Can only do pairwise comparisons (your factor can only have two variables in it)

# Replace "y" with the alpha diversity measure you want, "A" with the factor you want (can only do one as a time),
wilcox.test(y ~ A, data=EX_ps_clean.rar.rich.df)







### ANOVA for multiple categorical factors -------
# ANOVA is used when your factor has 2+ variables in it. If you have 3+ variables, you should run a post-hoc test to correct the p-value for multiple group comparisons

# One Way Anova (Completely Randomized Design)
fit <- aov(y ~ A, data=EX_ps_clean.rar.rich.df)

# Randomized Block Design (B is the blocking factor)
fit <- aov(y ~ A + B, data=EX_ps_clean.rar.rich.df)

# Two Way Factorial Design
fit <- aov(y ~ A*B, data=EX_ps_clean.rar.rich.df) 

# Analysis of Covariance
fit <- aov(y ~ A + x, data=EX_ps_clean.rar.rich.df) # x is a series of numerical values, like pH)


#example of one way anova from class
example1 <- lm(Observed ~ Diet, data=EX_ps_clean.rar.rich.df)
example2 <- lm(Observed ~ Week, data=EX_ps_clean.rar.rich.df)

#which is model better? Find the AIC value
anova(example1,example2)




#example of two way from class
summary(aov(Observed ~ Diet * Week, data=EX_ps_clean.rar.rich.df))

#For within subjects designs, the data frame has to be rearranged so that each measurement on a subject is a separate observation. See R and Analysis of Variance.

# One Within Factor
fit <- aov(y~A+Error(Subject/A),data=EX_ps_clean.rar.rich.df) #Subject would be a Factor in your metadata, y is the value of interest, A is the Factor of interest

# Two Within Factors W1 W2, Two Between Factors B1 B2
fit <- aov(y~(W1*W2*B1*B2)+Error(Subject/(W1*W2))+(B1*B2),
           data=EX_ps_clean.rar.rich.df)


# once you have chosen a model, add a correction if you have multiple group comparisons (multiple factor levels)
# Tukey Honestly Significant Differences
TukeyHSD(fit) # where fit comes from aov()




## If your alpha diversity is not normally distributed (non-parametric) -----------------
### Kruskal-Wallis Test -----
# K-W is the non-parametric version of ANOVA: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
kruskal.test(Observed ~ Diet, data = EX_ps_clean.rar.rich.df)
# data:  Observed by Diet
# Kruskal-Wallis chi-squared = 2.2931, df = 1, p-value = 0.13


# Follow it up with a Conover Test if you have multiple comparisons. Note, code changed recently
# https://rdrr.io/cran/DescTools/man/ConoverTest.html 
install.packages("conover.test")
library(conover.test)

conover.test(EX_ps_clean.rar.rich.df$Observed, EX_ps_clean.rar.rich.df$Week,
             method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) # method of correction, OK to pick just one

# output will look like this:
# Kruskal-Wallis rank sum test
# data: x and group
# Kruskal-Wallis chi-squared = 5.2867, df = 2, p-value = 0.07
# 
# Comparison of x by group                            
# (No adjustment)                                
# Col Mean-|
# Row Mean |          0          1
---------+----------------------
  #       1 |  -1.818304
  #         |     0.0417
#         |
#       2 |  -2.467699  -0.726045
#         |    0.0111*     0.2379
# 
# alpha = 0.05
# Reject Ho if p <= alpha/2
