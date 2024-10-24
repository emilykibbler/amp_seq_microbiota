
## Lab 1: Intro and package loading ---------------------------

setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal")


source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/AVS554_packages.R")
# install_necessary()
# install_optional()

load_libraries()


## Prelim data processing, splitting by read 1 and read 2  ---------------------------
# Separated into a different file

  

## Lab 2: Import and assess quality  ---------------------------

file_info <- metadata_import(types = c("Forward", "Reverse")) # saves a copy as an rds



# Look at the whole dataset to aggregate all samples onto one graph. 
# One plot for each data type will be in the wd in .png format
qualplots(file_info); beep() # Default data type is raw
# Notes: 68 samples files aggregated for both forward and reverse
  # Raw reads: 7035457

comparison <- as.data.frame(matrix(nrow = 0, ncol = 3))
colnames(comparison) <- c("Metric", "ESK", "Paper")
comparison[nrow(comparison) + 1, ] <- c("Raw", "7035457", "n.r.")

# run this command to create a quality control report
# data needs to be subsetted, too memory intensive to run all
# come back if I have time
for (i in 1:nrow(file_info)) {
  rqc(path = file_info$raw_data_dir[i][1], 
      sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
      file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())
}; beep("treasure")


filtoutput <- filterAndTrim( 
  file_info$file_names_full[[1]], 
  file.path(file_info$filt_dir[1], paste0(file_info$sample_names[[1]], "_F_filt.fastq.gz")), 
  file_info$file_names_full[[2]], 
  file.path(file_info$filt_dir[2], paste0(file_info$sample_names[[1]], "_R_filt.fastq.gz")), 
  trimLeft = c(20, 20), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
  # trimRight = c(50,50), # cuts of last XX bases from F,R reads, or hash this out and use truncLen
  truncLen = c(275, 250), # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
  maxEE = c(5,5), ## The maxEE parameter sets the maximum number of expected errors allowed in a read. Always >1
  verbose = FALSE) # verbose = print to screen
saveRDS(filtoutput, "filtoutput.rds")

as.data.frame(filtoutput)$reads.out/as.data.frame(filtoutput)$reads.in

qualplots(file_info, type = "filtered")

View(filtoutput)
### Summary of filtering ----
# Which parameters did you use and why? 

# I am trimming off 20 from the beginning of the reads
# When I trimmed only 10 bases, I lost ~50% of my reads later
# Github forum suggest trimming more of the beginning of the reads
# If primer sequences aren't removed, it can cause too many to be considered bimeras
# Forum also suggested relaxing truncation length and maxEE:
# with insufficient overlap, dada2 can't merge F-R
# 
# What was the largest and small number of raw reads?
# Largest: B166_R1.fastq.gz -- 249713
# Smallest: Neg16_R1.fastq.gz-- 108
# What was the largest and small number of filtered reads?
# Largest: B166_R1.fastq.gz -- 208273
# Smallest: Neg16_R1.fastq.gz -- 57
# How many total raw reads?
# 7035457 raw reads
# How many total filtered reads?
# I retained 6210016 reads
comparison[nrow(comparison) + 1, ] <- c("Total filtered reads", 6210016, 4150123)
median(as.data.frame(filtoutput)$reads.out[1:76])
comparison[nrow(comparison) + 1, ] <- c("Median exp_samp filtered reads", 72591, NA)
comparison[nrow(comparison) + 1, ] <- c("Minimum exp_samp filtered reads", 12559, NA)

## Lab 3---------------------------
  # DADA2 learn error rates (Lab 3) ----------------------------

# source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")


# file_info <- readRDS("file_info.rds")
create_and_plot_error_profile(file_info, bases = 1e6) ; beep("treasure") # writes the error profile RDS files
# Foward: 16906150 total bases in 73505 reads from 1 samples will be used for learning the error rates
# Reverse: 16455240 total bases in 91418 reads from 1 samples will be used for learning the error rates.
errF <- readRDS("error_profile_Forward.rds")
errR <- readRDS("error_profile_Reverse.rds")


dadaFs <- dada(list.files(file_info$filt_dir[1], full.names = TRUE),
               err = errF, 
               multithread = TRUE, verbose = FALSE) #; beep()
dadaRs <- dada(list.files(file_info$filt_dir[2], full.names = TRUE),
               err = errR, 
               multithread = TRUE, verbose = FALSE)


  ## DADA2 pick sequence variants (Lab 3) ----------------------


mergers <- mergePairs(dadaF = dadaFs,
                      derepF = list.files(file_info$filt_dir[1], full.names = TRUE),
                      dadaR = dadaRs,
                      derepR = list.files(file_info$filt_dir[2], full.names = TRUE),
                      verbose = FALSE) 
seqtab <- makeSequenceTable(mergers)
rowSums(seqtab)
saveRDS(seqtab, "seqtab.rds")
seqtab <- readRDS("seqtab.rds")

SVs_found_by_sample <- make_SV_summary(seqtab)
View(SVs_found_by_sample)
# Fewest SVs: Neg13_R1_F_filt.fastq.gz -- 1
# Most SVs: B213_R1_F_filt.fastq.gz -- 532

summary(SVs_found_by_sample$SVs[1:76])
comparison[nrow(comparison) + 1, ] <- c("Median ASV per exp_samp", 144.5, 80)



## Remove chimeras. Leave the method as consensus. -------------
# multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE) 
# My samples have Identified 875 bimeras out of 6019 input sequences.

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab_nochim) 
# 93 samples (row count) and 5144 SVs (column count)

# Calculate the percentage of chimeras identified out of the total
round(sum(seqtab_nochim)/sum(seqtab), digits = 3) # 0.996
sort(rowSums(seqtab_nochim))
# Experimental sample with fewest reads is B136, 6123 reads
median(rowSums(seqtab_nochim[1:76,]))
comparison[nrow(comparison) + 1, ] <- c("Median reads per exp_samp after dechim", 67171, NA)
summary(rowSums(seqtab_nochim[1:76,]))
comparison[nrow(comparison) + 1, ] <- c("Minimum reads per exp_samp after dechim", 6123, NA)

# 3. Save your chimera-free sequence table.
saveRDS(seqtab_nochim, 'seqtab_nochim.rds')



## Lab 4: Workflow verification  --------------------------
# This section can be optional, but as you are learning, or playing with new data, it is helpful to run some internal checks to assess whether you like how your quality control steps worked.

# 1. Track reads through the analysis here to see how many were lost at each QC step, in case you were too stringent. Only need to count F reads, even if you did paired F/R

# load the filtered output file if necessary
# filtoutput <- readRDS("filtoutput.rds")

# 2. Load the metadata file that goes with your sequencing data so you can match factors to seq data. This uses a datatable made in Excel and saved as a .csv file
meta <- read_xlsx("metadata.xlsx")
meta <- as.data.frame(meta)
head(meta)
row.names(meta) <- meta$SampleID
meta$Sample_type <- NA
for (i in 1:nrow(meta)) {
  if (meta$Group[i] == "controlneg") {
    meta$Sample_type[i] <- "negative"
  } else {meta$Sample_type[i] <- "experimental"}
}
saveRDS(meta, "meta.rds")


# 2.5 check the dimensions of the three data files you need for this to make sure the number of rows matches in each. If they do not, you may need to add/remove rows from your metadata file in case samples were removed/retained from your dataset.
dim(filtoutput) # 93 x 2
dim(seqtab_nochim) # 93 x 5144
dim(meta) # 98 x 46


head(row.names(filtoutput))
head(row.names(seqtab_nochim))
head(row.names(meta))

# Need to strip extra labeling

row.names(filtoutput) <- str_remove_all(row.names(filtoutput), "_R1.fastq.gz")
row.names(seqtab_nochim) <- str_remove_all(row.names(seqtab_nochim), "_R1_F_filt.fastq.gz")
saveRDS(seqtab_nochim,"seqtab_nochim.rds")
identical(row.names(filtoutput), row.names(seqtab_nochim)) # true, yay
identical(row.names(meta), row.names(filtoutput)) #false
length(row.names(meta)[!(row.names(meta) %in% row.names(filtoutput))]) # 5 more samples in metadata than in data
# Look at the sample metadata for what I will need to discard -- no match in dataset
View(meta[row.names(meta)[!(row.names(meta) %in% row.names(filtoutput))],])
# Look at the sample metadata for what I will keep
View(meta[row.names(meta)[row.names(meta) %in% row.names(filtoutput)],])


meta_cleaned <- meta[row.names(meta)[row.names(meta) %in% row.names(filtoutput)],]
identical(row.names(meta_cleaned), row.names(filtoutput)) # False
match(row.names(meta_cleaned), row.names(filtoutput)) # they do match, but in different orders
meta_cleaned <- meta_cleaned[order(row.names(meta_cleaned)),]
identical(row.names(meta_cleaned), row.names(filtoutput)) # true
identical(row.names(meta_cleaned), row.names(seqtab_nochim)) # true, good to go
saveRDS(meta_cleaned, "meta.rds")
saveRDS(filtoutput, "filtoutput.rds")
saveRDS(seqtab_nochim, "seqtab_nochim.rds")
# Write over meta df
meta <- meta_cleaned


# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab_nochim), meta_cleaned$Sample_type) 
# 4. Assign column names to that new variable
colnames(track) <- c("reads.in","filtered", "nonchimeras", "Sample_type") 

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'tracked_seqs.rds') 

# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))

### Plot QC steps ------------------

# plot with Sample_type along the X axis
ggplot(plotData,aes(x = Sample_type, y = as.numeric(totals))) + geom_jitter(aes(color = type)) + 
  ylab("Sequences") + 
  xlab("Sample_type") +
  # theme(axis.text.x = element_text(angle = 0, size = 12, vjust = -0.5)) +
  ggtitle("Reads by sample type")
ggsave("reads_filt_chim_sample_type.png")

# or, plot with QA stage along the X axis
ggplot(plotData,aes(x = type,y = as.numeric(totals))) + geom_jitter(aes(color = Sample_type)) + 
  ylab("Sequences") + 
  xlab("QA stage") +
  ggtitle("Reads by filtering step")
ggsave("reads_sample_type_QC_status.png")





# Lab 5 :Phyloseq First look (Lab 5_A) ------------------------------

# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal")
# load_lab_five_partA() # seqtab.nochim, meta

# # load data if necessary
seqtab_nochim <- readRDS('seqtab_nochim.rds')
meta <- readRDS("meta.rds")


# Check the sample sames and see if they match
identical(row.names(meta), row.names(seqtab_nochim))


## create a phyloseq object with all samples (controls will be subset out later)
phylo <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), # even though it's called an OTU table, it will use the SVs from my seqtab
                  sample_data(meta))

phylo
saveRDS(phylo, "phylo.rds")
# phylo <- readRDS("phylo.rds")
# how many samples made it this far? 93 samples, 46 variables, 5144 taxa. Taxa are the SV columns



# If you have sequenced controls, it's a good idea to look at your data in comparison to them.

# Alpha diversity peek 
plot_richness(phylo, x = "Sample_type", # CHANGE the x-axis to a factor of your choice
              measures = c("Observed","Chao1", "Shannon"), # these are some of the alpha diversity measures you can use 
              color = "Syndrome") + 
  geom_jitter() +
  theme_bw() + # a ggplot theme to make the graph look nice
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),) +
  ggtitle("Alpha diversity, before cleaning")

ggsave("initial_alpha_diversity_plot.png")


## Note: hopefully, diversity/richness is lower in the Negative controls than real samples



# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
plot(sort(taxa_sums(phylo), TRUE), 
     type = "h", 
     ylim = c(0, 20)) #limit the y-axis to better see the long tail of rare taxa


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
phylo_ord <- ordinate(phylo, #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

# If all the negatives are eliminated by filter and trim:
phylo_ord <- ordinate(prune_samples(sample_sums(phylo) > 0, phylo), #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)


plot_ordination(phylo, phylo_ord, type = "samples", color = "Group") +
  ggtitle("Ordination plot, pre-clean")
ggsave("ordination_before_plot.png") # save this graph for later

# Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues (see next section)
# Horse-shoe patterns indicate underlying patterns to the data
# The negative controls are supposed to be separate from the experimental samples and they are
# Out of 17 negatives controls provided, 13 have been filtered out so far



# Lab 5_B: Removing contamination using negative controls-----------------------------------
### This section can only be accomplished if you have sequenced negative controls and/or DNA quantification data from when the pool was created.  If you have these, CHOOSE ONE method to follow.  If you don't, skip this section and consider including sequencing controls in your next run.

phylo <- readRDS("phylo.rds") # read in phyloseq dataset if necessary
meta <- readRDS("meta.rds")
# load_lab_five_partB() # phylo and meta



## ------ Decontam method -------------
phylo <- readRDS("phylo.rds")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam) # This packages identify contaminants by frequency of SVs.

#### skipped, no quant data ----------
# The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?
# 

# The second contaminant identification method uses prevalence (presence/absence) across samples as compared to the prevalence in negative controls. Essentially, is it in most of your samples but not most of your controls?

# 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
sample_data(phylo)$is.neg <- sample_data(phylo)$Group == "controlneg"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(phylo, method = "prevalence", neg = "is.neg", threshold = 0.5)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))

# 11. get rid of the contaminants 
# remove them from the original phyloseq object, or the one with cleaned out SVs by frequency
phylo_clean_decontam_decontam <- prune_taxa(!contamdf.prev05$contaminant, phylo) 

# 12. Check how many are left
phylo_clean_decontam_decontam
# 93 samples and 5074 SVs left
saveRDS(phylo_clean_decontam_decontam, "phylo_clean_decontam_decontam.rds")

# 13. Did it work? Check your ordination again
phylo_clean_ord <- ordinate(phylo_clean_decontam_decontam, #calculate similarities
                            method = "PCoA", #ordination type
                            "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(phylo_clean_decontam_decontam, 
                phylo_clean_ord, 
                type = "samples", 
                color = "Group",
                title = "After cleaning with decontam")
ggsave("ord_plot_after_decontam.png")
# save this graph for later


median(rowSums(phylo_clean_decontam_decontam@otu_table[1:76,]))
comparison[nrow(comparison) + 1, ] <- c("Median reads per exp_samp after decontam method", 66853, 47323)
summary(rowSums(phylo_clean_decontam_decontam@otu_table[1:76,]))
comparison[nrow(comparison) + 1, ] <- c("Minimum reads per exp_samp after decontam method", 6112, NA)
median(rowSums(phylo_clean_decontam_decontam@otu_table[77:93,]))
comparison[nrow(comparison) + 1, ] <- c("Median reads per neg_ct after decontam method", 249, NA)
summary(rowSums(phylo_clean_decontam_decontam@otu_table[77:93,]))
comparison[nrow(comparison) + 1, ] <- c("Minimum reads per exp_samp after decontam method", 28, NA)



### Summary of decontamination ----
# CHANGE ME: Which types of negative controls did you use and why? 
# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?




# Lab 6_A: DADA2 assigning taxonomy from a reference database file (Lab 6) ---------------------------------


# Be sure to have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 

# The Silva file shown below is used for 16S rRNA (prokaryotes) and nicely formatted versions can be downloaded from the DADA2 website, which also contains some (18S) eukaryotic versions.  Let me know if you need a custom one for your dataset.
#
# 1. take the SVs (otu table) from the phyloseq object and make it into a dataframe again


phylo_clean_decontam_decontam <- readRDS("phylo_clean_decontam_decontam.rds")

seqtab_nochim_decontam = as(otu_table(phylo_clean_decontam_decontam), "matrix")
seqtab_nochim_decontam = as.matrix(seqtab_nochim_decontam)

## 2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops.
decontam_all_taxa <- assignTaxonomy(seqtab_nochim_decontam, 
                                    'silva_nr99_v138.1_train_set.fa.gz',        # CHANGE file path
                                    tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                    minBoot = 75, verbose = TRUE) # ; beep("treasure") 



saveRDS(decontam_all_taxa, 'decontam_all_taxa.rds')
write.csv(decontam_all_taxa, 'decontam_all_taxa.csv')
# decontam_all_taxa

# OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.
# decontam_all_taxa <- readRDS("decontam_all_taxa.rds")


decontam_all_taxa_species <- addSpecies(decontam_all_taxa, 
                                        'silva_species_assignment_v138.1.fa.gz', 
                                        allowMultiple = FALSE, 
                                        verbose = FALSE) ; beep("treasure") 

saveRDS(decontam_all_taxa_species, 'decontam_all_taxa_species.rds') 
write.csv(decontam_all_taxa_species, 'decontam_all_taxa_species.csv') 

### Remake the phyloseq object with the new taxonomy file ----------------------

# reload metadata table as needed
meta <- readRDS("meta.rds")
phylo_clean_decontam_decontam <- readRDS("phylo_clean_decontam_decontam.rds")

#reload taxa table as needed
decontam_all_taxa_species <- readRDS('decontam_all_taxa_species.rds')


otu_t <- otu_table(phylo_clean_decontam_decontam)

## create a phyloseq object with all samples
phylo_decontam_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                        sample_data(meta),
                                        tax_table(decontam_all_taxa_species))


saveRDS(phylo_decontam_with_species, "phylo_decontam_with_species.rds")



# Lab 6_B: Clean out unwanted taxa ------------------------------------
# Regardless of whether you included sequenced negative controls, you can remove taxa which are of no interest to you

# clean out chloroplast and mitochondria. can also elect to remove other contaminating domains or kingdoms as needed. 

# phylo_clean_with_species <- readRDS("phylo_clean_with_species.rds")
# Optional: explore your taxonomy before filtering. Use the tax table you made

df <- as.data.frame(phylo_decontam_with_species@tax_table)
table(df$Kingdom)
table(df$Phylum)
table(df$Class)
table(df$Order)
table(df$Family)

phylo_decontam_with_species <- phylo_decontam_with_species %>% # CHANGE ME to your phyloseq object name 
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast") # CHANGE ME to taxa you want to remove

saveRDS(phylo_decontam_with_species, "phylo_decontam_with_species.rds")

# Write out a description of experimental design (Homework)
# FIXME


# Lab 7: Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y



# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
#

tab <- otu_table(phylo_decontam_with_species)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
# rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)


r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_decontam_with_species)))
# smallest reads (sequences) in a sample ______ B136 -- 2922 
# largest in a sample ______ B166 197735
# number of experimental samples with <5000 ______ 1
# In the paper they rarified to 12k, I don't understand how I lost so many reads


phylo_decontam_rar <- rarefy_even_depth(phylo_decontam_with_species, 
                                        sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                        replace = TRUE, #sampling with or without replacement
                                        trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                        rngseed = 711, 
                                        verbose = TRUE)
# set.seed(711) was used
# 17 negative control samples removed and one experimental sample
# 363 OTUs removed


# this helps with plotting later
sample_data(phylo_decontam_rar)$Group <- factor(sample_data(phylo_decontam_rar)$Group, 
                                                levels = c("IPD_ATB", "IPD"), 
                                                labels = c("Antibiotics", "No antibiotics")) #CHANGE ME

saveRDS(phylo_decontam_rar, "decontam_phylo_rarified.rds")

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(phylo_decontam_rar), 'decontam_phylo_rarified.csv')

# Specifying taxa after rarification

phylo_decontam_strep_rar  <- phylo_decontam_rar  %>% 
  subset_taxa(Family == "Streptococcaceae") 

phylo_decontam_no_strep_rar   <- phylo_decontam_rar  %>% 
  subset_taxa(Family != "Streptococcaceae") 


## clean out taxa/SV columns that are no longer present
phylo_decontam_strep_rar <- prune_taxa(taxa_sums(phylo_decontam_strep_rar) > 0, phylo_decontam_strep_rar)
phylo_decontam_strep_rar <- prune_samples(sample_sums(phylo_decontam_strep_rar) > 0, phylo_decontam_strep_rar)
phylo_decontam_strep_rar
# 73 samples, 64 taxa

phylo_decontam_no_strep_rar <- prune_taxa(taxa_sums(phylo_decontam_no_strep_rar) > 0, phylo_decontam_no_strep_rar)
phylo_decontam_no_strep_rar <- prune_samples(sample_sums(phylo_decontam_no_strep_rar) > 0, phylo_decontam_no_strep_rar)
phylo_decontam_no_strep_rar
# 75 samples, 1457 taxa


saveRDS(phylo_decontam_strep_rar, 'phylo_decontam_strep_rar.rds') 
saveRDS(phylo_decontam_no_strep_rar, 'phylo_decontam_no_strep_rar.rds') 


# meta <- readRDS("meta.rds")


# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

### Alpha diversity plotting ------------------

# This script was getting too long and busy so I have split the plotting
# See separate .R file named nasal_aspirate_plots





## Alpha diversity plotted against other metadata -----

# use phyloseq to measure alpha diversity
phylo_rar_rich <- estimate_richness(phylo_decontam_rar, measures = c("Observed", "Shannon")) #change to whatever measures you want

# # OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
  # install.packages("remotes")
  # remotes::install_github("twbattaglia/btools")
  # EX_faith <- estimate_pd(EX_ps_clean.rar)

# measure evenness for each sample
phylo_rar_even <- phylo_rar_rich$Shannon/log(phylo_rar_rich$Observed)

# Coerce to data.frame and add the metadata for these samples
phylo_rar_sd = as(sample_data(phylo_decontam_rar), "matrix")
phylo_rar_sd = as.data.frame(phylo_rar_sd)
phylo_rar_rich_df <- cbind(phylo_rar_rich, phylo_rar_even, phylo_rar_sd)

dim(phylo_rar_rich_df) # 66 x 49

# make a graph using that dataframe. this is a generic example, you will want to personalize this.
ggplot(phylo_rar_rich_df, aes(x = Group, y = Observed)) + 
  theme_minimal() + 
  geom_boxplot(aes(color = Group)) +
  theme(legend.position = "none") +
  xlab("") + 
  ylab("Bacterial Richness (SVs)") + 
  # theme(text = element_text(size = 20)) # increases font size
  ggtitle("Bacterial richness")

### skipped, maybe come back later --------

# # make that same graph but drop any samples that lack data for that FactorA (replace FactorA in the code with your factor name).
# ggplot(data=subset(EX_ps_clean.rar.rich.df, !is.na(FactorA)), aes(x=Temperature, y=Observed)) + 
#   theme_minimal() + 
#   geom_point(aes(color=Group), size = 3) + # sets colors by group and a set point size 
#   xlab("Temperature of Ocean Water (C)") + 
#   ylab("Bacterial Richness (SVs)") + 
#   theme(text = element_text(size = 20)) # increases font size



# Make an alpha diversity table: https://www.geeksforgeeks.org/create-table-from-dataframe-in-r/
alpha_diversity_table = as.table(table(phylo_rar_rich_df$Group, phylo_rar_rich_df$Gestational_age)) 
alpha_diversity_table


# -----------
  # Write out your research questions, make them specific. 
  # Write a description of your experimental design, including notes about replication, repeated measures, nested design, or other factors which might impact your statistical models.
  
  
  
  
  
  # Alpha diversity metrics statistics (Lab 8)--------------
# Phyloseq can measure and visualize alpha diversity: https://joey711.github.io/phyloseq/plot_richness-examples.html
# phyloseq doesn't do stats or more complex graphs


# use phyloseq to measure alpha diversity

phylo_decontam_rar_rich <- estimate_richness(phylo_decontam_rar, measures = c("Observed", "Shannon")) #change to whatever measures you want

# use phyloseq to calculate Faith's Diversity metric (optional), https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# EX_faith <- estimate_pd(EX_ps_clean.rar)

# OPTIONAL measure evenness for each SV individually

phylo_decontam_rar_even_SV <- evenness(otu_table(phylo_decontam_rar, taxa_are_rows = FALSE))
  


test_phylo_decontam_rar_rich <- normal_stats(phylo_decontam_rar_rich, phylo_decontam_rar)
# saveRDS("phylo_decontam_rar_rich")

 ## If your alpha diversity is not normally distributed (non-parametric) -----------------
### Kruskal-Wallis Test -----
# K-W is the non-parametric version of ANOVA: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
kruskal.test(Observed ~ Group, data = phylo_decontam_rar_rich)
# data:  Observed by Diet
# Kruskal-Wallis chi-squared = 10.823, df = 1, p-value = 0.001002


# Follow it up with a Conover Test if you have multiple comparisons. Note, code changed recently
# https://rdrr.io/cran/DescTools/man/ConoverTest.html 


conover.test(phylo_decontam_rar_rich$Observed, phylo_decontam_rar_rich$Group,
             method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) # method of correction, OK to pick just one

# output will look like this:
# Kruskal-Wallis rank sum test
# data: x and group
# Kruskal-Wallis chi-squared = 10.823, df = 1, p-value = 0
# 
# 
# Comparison of x by group                            
# (No adjustment)                                
# Col Mean-|
#   Row Mean |   Antibiot
# ---------+-----------
#   No antib |   3.536360
# |    0.0004*
#   
#   alpha = 0.05
# Reject Ho if p <= alpha/2

### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html 



# you might run each factor separately to find out if they are significant on their own. 

numeric_columns <- c("Age", "Birth_weight", "Gestational_age", "House_surface", "Breastfeeding_time", "Pneumococcal_load", "Fever_time_before_sampling", "C_reactive_protein" , "Hemoglobin", "Leukocytes", "Hospitalization_days" )
df <- phylo_decontam_rar_rich
df[,numeric_columns] <- sapply(df[numeric_columns],as.numeric)
sapply(df, class)
numeric_columns <- colnames(df)[which(as.data.frame(sapply(df, class))[,1] == "numeric")]

summary(lm(Observed ~ ., data = df[, numeric_columns])) 
df_ab <- subset(df, Group == "Antibiotics")
summary(lm(Observed ~ ., data = df_ab[, numeric_columns])) 
df_no_ab <- subset(df, Group == "No antibiotics")
summary(lm(Observed ~ ., data = df_no_ab[, numeric_columns])) 

# what if you want to add log transformation to a numerical variable in your metadata? Check out the syntax here: https://rpubs.com/marvinlemos/log-transformation


### Linear mixed effects models for complicated experimental designs -----
# If you have more complicated experimental designs, you might need lmer or glmer models to accommodate that
# does not have a normal distribution, so compare means using glmer and poisson distribution (number of events within a time interval): https://towardsdatascience.com/the-poisson-distribution-and-poisson-process-explained-4e2cb17d459

# lm <- (glmer(Observed ~ A + C + (1|Year_collected), family=poisson, data=EX_ps_clean.rar.rich.df))
# 
# summary(lm)
# paste the output here


# if you have multiple groups, this will give you the pairwise comparisons
model <- lm(Observed ~ ., data = df[, numeric_columns])
emmeans(model, pairwise ~ Age) 
# paste the output here

emmeans(lm,pairwise ~ B) 
# paste the output here


