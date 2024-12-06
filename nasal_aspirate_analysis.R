
## Lab 1: Intro and package loading ---------------------------

setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal")


source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab8_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab9_functions.R")
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
summary(as.data.frame(filtoutput)$reads.out[1:76])
comparison[nrow(comparison) + 1, ] <- c("Median exp_samp filtered reads", 72591, NA)
comparison[nrow(comparison) + 1, ] <- c("Minimum exp_samp filtered reads", 12559, NA)

## Lab 3: DADA2 learn error rates (Lab 3) ----------------------------


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

SVs_found_by_sample <- make_SV_summary(seqtab) # this function is in lab3_functions
summary(SVs_found_by_sample$SVs[1:76])
View(SVs_found_by_sample)
# Fewest SVs: Neg13_R1_F_filt.fastq.gz -- 1
# Most SVs: B213_R1_F_filt.fastq.gz -- 532

summary(SVs_found_by_sample$SVs[1:76])
comparison[nrow(comparison) + 1, ] <- c("Median ASV per exp_samp", 144.5, 80)



## Remove chimeras.-------------
# Leave the method as consensus. 
# multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE) 
# My samples have Identified 875 bimeras out of 6019 input sequences.
  # View(make_SV_summary(seqtab_nochim))
  # sum(make_SV_summary(seqtab_nochim)$SVs)
  # summary(make_SV_summary(seqtab_nochim)$SVs[1:76])
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
#### Patient metadata analysis -----
pts_meta <- subset(meta, Group != "controlneg")
# single check to make sure I can run this function on categorical data
cor(as.numeric(as.factor(pts_meta[,"Group"])), as.numeric(as.factor(pts_meta[,"Gender"])), method = "kendall") # crushed it
numeric_columns <- c("Age", "Birth_weight", "Gestational_age", "House_surface", "Breastfeeding_time", "Pneumococcal_load", "Fever_time_before_sampling", "C_reactive_protein" , "Hemoglobin", "Leukocytes", "Hospitalization_days" , "Days_ATB_before_sampling")
categorical_pt_data <- pts_meta[, !(names(pts_meta) %in% numeric_columns)]
categorical_pt_data <- subset(categorical_pt_data, select = -Sample_type)
numeric_pt_data <- pts_meta[, names(pts_meta) %in% c("Group", numeric_columns)]

patient_data_correlation_summary <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(patient_data_correlation_summary) <- c(names(unlist(cor.test(as.numeric(as.factor(meta[,"Group"])), as.numeric(meta[,"Birth_weight"]), method = "pearson", na.action = "na.remove"))), "Variable")

for (i in 1:ncol(numeric_pt_data)) {
  if (colnames(numeric_pt_data)[i] != "Group") {
  temp <- subset(numeric_pt_data, !is.na(colnames(numeric_pt_data)[i]))
  patient_data_correlation_summary[nrow(patient_data_correlation_summary) + 1, ] <-
    c(unlist(cor.test(as.numeric(as.factor(temp[,"Group"])), as.numeric(temp[,i]), method = "pearson")), colnames(numeric_pt_data)[i])
  }
}

saveRDS(patient_data_correlation_summary, "numeric_patient_data_correlation_summary.rds")

# Chi square for categorical data

# chisq.test(as.factor(categorical_pt_data$Group), as.factor(categorical_pt_data$Ethnicity)) # example for me to check the format of the output
# unlist(chisq.test(as.factor(categorical_pt_data$Group), as.factor(categorical_pt_data$Ethnicity)))[1:4]

chisq_summary <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(chisq_summary) <- c("statistic.X-squared", "parameter.df", "p.value", "method", "variable")
for (i in 1:ncol(categorical_pt_data)) {
  temp <- subset(categorical_pt_data, !is.na(categorical_pt_data[,i]))
  chisq_summary[nrow(chisq_summary) + 1, ] <-
    c(unlist(chisq.test(as.factor(temp[,"Group"]), as.factor(temp[,i])))[1:4], colnames(categorical_pt_data)[i])
}

saveRDS(chisq_summary, "chisq_summary.rds")
### Filtering steps analysis ---------
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


meta <- meta[row.names(meta)[row.names(meta) %in% row.names(filtoutput)],]
identical(row.names(meta), row.names(filtoutput)) # False
match(row.names(meta), row.names(filtoutput)) # they do match, but in different orders
meta <- meta[order(row.names(meta)),]
identical(row.names(meta), row.names(filtoutput)) # true
identical(row.names(meta), row.names(seqtab_nochim)) # true, good to go
saveRDS(meta, "meta.rds")
saveRDS(filtoutput, "filtoutput.rds")
saveRDS(seqtab_nochim, "seqtab_nochim.rds")



# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab_nochim), meta$Sample_type) 
# 4. Assign column names to that new variable
colnames(track) <- c("reads.in","filtered", "nonchimeras", "Sample_type") 

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'track.rds') 

# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))

### Plot QC steps ------------------
# See plot script




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
phylo <- readRDS("phylo.rds")
# how many samples made it this far? 93 samples, 46 variables, 5144 taxa. Taxa are the SV columns



# If you have sequenced controls, it's a good idea to look at your data in comparison to them.

# Alpha diversity peek 
# see plot script


## Note: hopefully, diversity/richness is lower in the Negative controls than real samples. Check!



# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
plot(sort(taxa_sums(phylo), TRUE), 
     type = "h", 
     ylim = c(0, 20)) #limit the y-axis to better see the long tail of rare taxa


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
phylo_ord <- ordinate(phylo, #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

# See plotting script

# Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues (see next section)
# Horse-shoe patterns indicate underlying patterns to the data
# The negative controls are supposed to be separate from the experimental samples and they are




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

# skipped, no quant data 
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
summary(rowSums(phylo_clean_decontam_decontam@otu_table)[1:76])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6112   43586   66853   70767   92416  200156 

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
decontam_all_taxa <- readRDS("decontam_all_taxa.rds")


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


otu_t <- otu_table(phylo_clean_decontam_decontam, taxa_are_rows = FALSE)

## create a phyloseq object with all samples
phylo_decontam_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                        sample_data(meta),
                                        tax_table(decontam_all_taxa_species))


saveRDS(phylo_decontam_with_species, "phylo_decontam_with_species.rds")
summary(rowSums(phylo_decontam_with_species@otu_table)[1:76])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 6112   43586   66853   70767   92416  200156 

# Lab 6_B: Clean out unwanted taxa ------------------------------------
# Regardless of whether you included sequenced negative controls, you can remove taxa which are of no interest to you

# clean out chloroplast and mitochondria. can also elect to remove other contaminating domains or kingdoms as needed. 

phylo_decontam_with_species <- readRDS("phylo_decontam_with_species.rds")
# Optional: explore your taxonomy before filtering. Use the tax table you made

df <- as.data.frame(phylo_decontam_with_species@tax_table)
table(df$Kingdom) # 1 eukaryotic SV to be removed
table(df$Phylum)
table(df$Class)
table(df$Order) # 17 chloroplast to be removed
table(df$Family) # 724 mitochondria

phylo_decontam_with_species <- phylo_decontam_with_species %>% 
  subset_taxa(Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast")
summary(rowSums(phylo_decontam_with_species@otu_table)[1:76])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2922   40793   63951   67733   90951  197735 
saveRDS(phylo_decontam_with_species, "phylo_decontam_with_species.rds")

# Write out a description of experimental design (Homework)
# See manuscript


# Lab 7: Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y



# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
#

tab <- otu_table(phylo_decontam_with_species, taxa_are_rows = FALSE)
class(tab) <- "matrix" # tab <- as.matrix(tab) will do nothing for some reason
## you get a warning here, ignore, this is what we need to have
tab <- t(tab) # transpose observations to rows

r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_decontam_with_species)))
summary((rowSums(otu_table(phylo_decontam_with_species))[1:76]))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2922   40793   63951   67733   90951  197735 
# smallest reads (sequences) in a sample ______ B136 -- 2922 
# largest in a sample ______ B166 197735
# number of experimental samples with <5000 ______ 1
# In the paper they rarified to 12k, I don't understand how I lost so many reads

dim(phylo_decontam_with_species@otu_table)
# 93, 1884
summary(make_SV_summary(as.data.frame(phylo_decontam_with_species@otu_table))$SVs[1:76])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 5.00   19.00   39.50   59.76   83.50  265.00 

phylo_decontam_rar <- rarefy_even_depth(phylo_decontam_with_species, 
                                        sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                        replace = TRUE, #sampling with or without replacement
                                        trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                        rngseed = 711, 
                                        verbose = TRUE)
# set.seed(711) was used
# 17 negative control samples removed and one experimental sample
# 363 OTUs removed
phylo_decontam_rar
# 1521 taxa, 75 samples


# this helps with plotting later
sample_data(phylo_decontam_rar)$Group <- factor(sample_data(phylo_decontam_rar)$Group, 
                                                levels = c("IPD_ATB", "IPD"), 
                                                labels = c("Antibiotics", "No antibiotics")) 

saveRDS(phylo_decontam_rar, "phylo_decontam_rar.rds")

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(phylo_decontam_rar, taxa_are_rows = FALSE), 'decontam_phylo_rarified.csv')

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



# Write out your research questions, make them specific. 
  # Write a description of your experimental design, including notes about replication, repeated measures, nested design, or other factors which might impact your statistical models.
# See manuscript https://docs.google.com/document/d/1q4icwy0ZIGTJOfsY6cBDnvcO1PkW_Xpve7MGYT9Gw9k/edit?usp=sharing
# Document is currently only available to members of the UMaine system, as it is a work in progress
# The best parts of the report will be added to the GitHub readme
  
  
  
  
  # Lab 8: Alpha diversity metrics statistics (Lab 8)--------------
# Phyloseq can measure and visualize alpha diversity: https://joey711.github.io/phyloseq/plot_richness-examples.html
# phyloseq doesn't do stats or more complex graphs



# use phyloseq to measure alpha diversity

phylo_decontam_rar_rich <- estimate_richness(phylo_decontam_rar, measures = c("Chao1", "Observed", "Shannon")) #change to whatever measures you want
saveRDS(phylo_decontam_rar_rich, "phylo_decontam_rar_rich.rds")
# use phyloseq to calculate Faith's Diversity metric (optional), https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# EX_faith <- estimate_pd(EX_ps_clean.rar)

# OPTIONAL measure evenness for each SV individually
  phylo_decontam_rar_even_SV <- evenness(otu_table(phylo_decontam_rar, taxa_are_rows = FALSE))

phylo_decontam_rar_rich_df <- normal_stats(phylo_decontam_rar) # see lab8_functions.R
# Outputs:
  # Kurtosis: 1.334661
  # Shannon shapiro: W = 0.96103, p-value = 0.02104
  # Observed shapiro: W = 0.86331, p-value = 9.1e-07 -- I later update my richness metric from Observed to Chao1 to match the paper
  # Even shapiro: W = 0.96621, p-value = 0.04258

# Outputs:
  # Kurtosis: 1.940676
  # Shannon shapiro: W = 0.96103, p-value = 0.02104
  # Chao1 shapiro: W = 0.86163, p-value = 7.969e-07
  # Even shapiro: W = 0.96756, p-value = 0.05132


saveRDS(phylo_decontam_rar_rich_df, "phylo_decontam_rar_rich_df.rds")

 ## If your alpha diversity is not normally distributed (non-parametric) -----------------
### Kruskal-Wallis Test -----
# K-W is the non-parametric version of ANOVA: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
kruskal.test(Chao1 ~ Group, data = phylo_decontam_rar_rich_df)
# data:  Chao1 by Group
# Kruskal-Wallis chi-squared = 10.358, df = 1, p-value = 0.001289
kruskal.test(Shannon ~ Group, data = phylo_decontam_rar_rich_df)
# data:  Shannon by Group
# Kruskal-Wallis chi-squared = 4.7354, df = 1, p-value = 0.02955

# Follow it up with a Conover Test if you have multiple comparisons. Note, code changed recently
# https://rdrr.io/cran/DescTools/man/ConoverTest.html 


conover.test(phylo_decontam_rar_rich_df$Observed, phylo_decontam_rar_rich_df$Group,
             method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) # method of correction, OK to pick just one
# Output:
    # Kruskal-Wallis rank sum test
    # data: x and group
    # Kruskal-Wallis chi-squared = 10.823, df = 1, p-value = 0
    # 
    # 
    # Comparison of x by group                            
    # (No adjustment)                                
    # Col Mean-|
    #   Row Mean |   Antibiot
    #            +
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
phylo_decontam_rar_rich_df <- phylo_decontam_rar_rich_df
phylo_decontam_rar_rich_df[,numeric_columns] <- sapply(phylo_decontam_rar_rich_df[numeric_columns],as.numeric)
sapply(phylo_decontam_rar_rich_df, class)
numeric_columns <- colnames(phylo_decontam_rar_rich_df)[which(as.data.frame(sapply(phylo_decontam_rar_rich_df, class))[,1] == "numeric")]

summary(lm(Observed ~ ., data = phylo_decontam_rar_rich_df[, numeric_columns])) 
phylo_decontam_rar_rich_df_ab <- subset(phylo_decontam_rar_rich_df, Group == "Antibiotics")
summary(lm(Observed ~ ., data = df_ab[, numeric_columns])) 
df_no_ab <- subset(phylo_decontam_rar_rich_df, Group == "No antibiotics")
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
model <- lm(Observed ~ ., data = df[, c("Group", numeric_columns)])
emmeans(model, pairwise ~ Group) 
# paste the output here
  
  # emmeans
  # Group          emmean   SE df lower.CL upper.CL
  # Antibiotics      48.4 2.13 44     44.1     52.7
  # No antibiotics   51.6 3.30 44     44.9     58.2
  # 
  # Confidence level used: 0.95 
  # 
  # $contrasts
  # contrast                     estimate   SE df t.ratio p.value
  # Antibiotics - No antibiotics    -3.16 4.22 44  -0.749  0.4578

atb_df <- subset(phylo_decontam_rar_rich_df, Group == "Antibiotics")
no_atb_df <- subset(phylo_decontam_rar_rich_df, Group == "No antibiotics")
t.test(atb_df$Observed, no_atb_df$Observed) # different
t.test(atb_df$Shannon, no_atb_df$Shannon) # not significantly different


# Lab 9: More ways to graph alpha diversity  -------
# Resources for visualization
# https://www.data-to-viz.com/
# https://r-graph-gallery.com/index.html

### Heatmaps --------------
# base plot
plot_heatmap(phylo_decontam_rar)
# follow the tutorial to make prettier versions in phyloseq: https://joey711.github.io/phyloseq/plot_heatmap-examples.html
# go to the plot script to see code for a fancier heatmaps

# example from class
phylo_decontam_rar_glom = tax_glom(phylo_decontam_rar, "Phylum")




#example from class, but looks like trash without a lot of fancying up
# heatmap(otu_table(EX_ps_clean.rar))


### Correlogram --------------
# This makes a correlation matrix plot: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")
# phylo_decontam_rar_rich <- readRDS("phylo_decontam_rar_rich.rds") # this is the same as below but less
phylo_decontam_rar_rich_df <- readRDS("phylo_decontam_rar_rich_df.rds") # this contains the same info as above but sample metadata too

# corrplot requires one dataframe of metadata and/or seqtab data.  Any columns in the dataframe will be correlated against all others.  Too many columns makes the graph unreadable, try to keep it to <50.

# select to top 15 most abundant SVs from your phyloseq object and extract to a dataframe
#take out top 100 taxa based on abundance
top100_decontam_rar_corr_df <- create_top_N_corr_df(phylo_decontam_rar, 100, phylo_decontam_rar_rich_df) # see lab9_functions
# check header to make sure it looks good
head(top100_decontam_rar_corr_df)
# Top 100 don't all have reads after rare



# This is from browsing the metadata
metadata_numeric_columns <- c("Age", "Birth_weight", "Gestational_age", "House_surface", "Breastfeeding_time", "Pneumococcal_load", "Fever_time_before_sampling", "C_reactive_protein" , "Hemoglobin", "Leukocytes", "Hospitalization_days" )

top30_corr_df <- create_top_N_corr_df(phylo_decontam_rar, 30, phylo_decontam_rar_rich_df) # see lab9_functions
top30_corr_df[,metadata_numeric_columns] <- sapply(top30_corr_df[metadata_numeric_columns],as.numeric)
numeric_columns <- colnames(top30_corr_df)[which(as.data.frame(sapply(top30_corr_df, class))[,1] == "numeric")]
sapply(top30_corr_df[,numeric_columns], class)

# alternatively, load a premade dataframe containing your chosen SVs from your sequence table (otu)table in phyloseq) and the metadata you want to include
# EX_ps_clean.rar.corr.df <- read.csv("example_correlogram_dataframe.csv", check.names = FALSE, header=T, row.names=1)

# run correlations
corr_calc <- cor(top30_corr_df[,numeric_columns],
                    use = "complete.obs", # use=Specifies the handling of missing data. Options are all.obs (assumes no missing data - missing data will produce an error), complete.obs (listwise deletion), and pairwise.complete.obs (pairwise deletion)
                    method = "spearman") # correlation method=pearson, spearman, or kendall
# corr.mtest is a custom function written by Dr. Ishaq, see lab9_functions file
res1 <- cor.mtest(top30_corr_df[,numeric_columns], 0.95)

### Note, if you have too few samples, you may receive an error about too few observations. 
# You may ignore the error message, or remove than column


# Note: if you come up with a warning message, 
# it means that one or more of your columns generated correlations of value 0, which makes it angry.  
# Visualize the plot, and "?" will come up in those columns.  
# To fix, remove the column or add a data transformation, to the dataframe you created to make this.  
# Set the corr.mtest function and run corr.mtest again.

# res2 <- cor.mtest(EX_ps_clean.rar.corr.df,0.99)

colnames(res1[[1]]) <- colnames(top30_corr_df[,numeric_columns])
rownames(res1[[1]]) <- colnames(top30_corr_df[,numeric_columns])

## plot only those correlations with a significant p-value <0.05, using hierarchical clustering
corrplot(corr_calc, 
         type = "lower", #shape of the plot itself: full, upper, lower
         method = "circle", #shape of the data: circle, square, ellipse, number, shade, color, pie 
         order = "hclust", #how to cluster samples: AOE, hclust, FPC, alphabet, or leave blank
         p.mat = res1[[1]], #which correlations to use
         sig.level = 0.05, #sets significance cutoff
         insig = "blank", #leaves p > 0.05 blank
         tl.col = "black", # text color
         tl.cex = .6, #text size
         col = brewer.pal(n = 10, name = "RdYlBu")) #specify a color palette and number of shades needed






### Barplots------------------
# can add ggplot components to make it pretty
#don't recommend using genus here, it make crash R

# see plotting script



# example from class
  # EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")
  # plot_bar(EX_ps_clean.rar.glom, fill="Phylum") +   
  #   facet_grid(~Diet, space="free", scales="free") + 
  #   theme(legend.position = "bottom", axis.text.x = element_blank()) 

# make a stacked 100% bar chart in phyloseq

  # EX_ps_clean.rar.stacked = transform_sample_counts(EX_ps_clean.rar, function(x) x / sum(x) )
  # plot_bar(EX_ps_clean.rar.stacked, fill="Phylum") 

# to filter by abundance and pool low abundance groups: https://github.com/joey711/phyloseq/issues/901


## Core community members -----------
#code from this website: https://microbiome.github.io/tutorials/CoremicrobiotaAmplicon.html



# graph the abundance of those shared taxa, here are some example: https://microbiome.github.io/tutorials/Core.html
source("~/Desktop/projects/R/AVS_554/lab9_functions.R")

# test_phylo.coreW <- core_taxa_finder(phylo_decontam_rar, c(1/10000, 25/100))
# identical(test_phylo.coreW, phylo.coreW)
phylo.coreW_25 <- core_taxa_finder(phylo_decontam_rar, c(1/10000, 25/100))
saveRDS(phylo.coreW_25, "phylo.coreW_25.rds")



phylo.coreW_35 <- core_taxa_finder(phylo_decontam_rar, c(1/10000, 35/100))
saveRDS(phylo.coreW_35, "phylo.coreW_35.rds")

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)
atb_phylo.coreW_35 <- subset_samples(phylo.coreW_35, Group == "Antibiotics")
saveRDS(atb_phylo.coreW_35, "atb_phylo.coreW_35.rds")
no_atb_phylo.coreW_35 <- subset_samples(phylo.coreW_35, Group == "No antibiotics")
saveRDS(no_atb_phylo.coreW_35, "no_atb_phylo.coreW_35.rds")

# Heat maps and boxplots of the above data can be found in the plot script
phylo.coreW_35 <- readRDS("phylo.coreW_35.rds")
atb_phylo.coreW_35 <- readRDS("atb_phylo.coreW_35.rds")
no_atb_phylo.coreW_35 <- readRDS("no_atb_phylo.coreW_35.rds")

for (i in 1:ncol(phylo.coreW_35@otu_table)) {
  print(colnames(atb_phylo.coreW_35@otu_table)[i])
  print(t.test(atb_phylo.coreW_35@otu_table[,i], no_atb_phylo.coreW_35@otu_table[,i]))
}

SVs_sig_diff_on_t_test <- c("GGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGCTTATGGTTAATACCCATAAGCCCTGACGTTACCCACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTATTTAAGTCAGATGTGAAAGCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGGATAACTAGAGTAGGTGAGAGGGGAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTCCCTGGCATCATACTGACACTGAGGTGCGAAAGCGTGGGTAGCAAACAG",
                            "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                            "GGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACAAATGTGTAAGTAACTATGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG")
# in case I got confused and pasted something twice
SVs_sig_diff_on_t_test <- unique(SVs_sig_diff_on_t_test)

aligns_to_AY_in_paper <- "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"

t.test(atb_phylo.coreW_35@otu_table[,aligns_to_AY_in_paper], no_atb_phylo.coreW_35@otu_table[,aligns_to_AY_in_paper])

aligns_to_HG_in_paper <- "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
t.test(atb_phylo.coreW_35@otu_table[,aligns_to_HG_in_paper], no_atb_phylo.coreW_35@otu_table[,aligns_to_HG_in_paper])

### Simple method ------------

# this is the one with SVs as colnames, sample as row names, and read numbers as values
view(phylo_decontam_rar@otu_table)
# this is the one with SVs as row names, classifications as colnames, class names as values
view(phylo_decontam_rar@tax_table)

phylo_decontam_rar_atb <- subset_samples(phylo_decontam_rar, Group == "Antibiotics")
phylo_decontam_rar_no_atb <- subset_samples(phylo_decontam_rar, Group == "No antibiotics")

phylo_decontam_rar_atb_otu <- as.data.frame(phylo_decontam_rar_atb@otu_table)
dim(phylo_decontam_rar_atb_otu)
phylo_decontam_rar_no_atb_otu <- as.data.frame(phylo_decontam_rar_no_atb@otu_table)
dim(phylo_decontam_rar_no_atb_otu)

sig_SVs <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(sig_SVs) <- c("SV", "p-value")

# unlist(t.test(phylo_decontam_rar_atb_otu[,1], phylo_decontam_rar_no_atb_otu[,1]))

for (i in 1:ncol(phylo_decontam_rar_atb_otu)) {
  if (unlist(t.test(phylo_decontam_rar_atb_otu[,i], phylo_decontam_rar_no_atb_otu[,i]))[3] < 0.05) {
    sig_SVs[nrow(sig_SVs) + 1, ] <- c(colnames(phylo_decontam_rar_atb_otu)[i], unlist(t.test(phylo_decontam_rar_atb_otu[,i], phylo_decontam_rar_no_atb_otu[,i]))[3])
  }
}

dim(sig_SVs)
tax_key <- as.data.frame(phylo_decontam_rar@tax_table)
tax_key$SV <- row.names(tax_key)
sig_SVs <- left_join(sig_SVs, tax_key, by = "SV")

phylo_decontam_rar_abundance <- microbiome::transform(phylo_decontam_rar, "compositional")
phylo_decontam_rar_abundance_df <- as.data.frame(phylo_decontam_rar_abundance@otu_table)
abundance_means <- as.data.frame(cbind(colnames(phylo_decontam_rar_abundance_df), colMeans(phylo_decontam_rar_abundance_df)))
colnames(abundance_means) <- c("SV", "mean_abundance")

phylo_decontam_rar_abundance_df$SampleID <- row.names(phylo_decontam_rar_abundance_df)
phylo_decontam_rar_abundance_df <- pivot_longer(phylo_decontam_rar_abundance_df, 
                                                cols = colnames(phylo_decontam_rar_abundance_df)[1:length(phylo_decontam_rar_abundance_df) - 1],
                                                names_to = "SV",
                                                values_to = "Abundance")
phylo_decontam_rar_abundance_df <- left_join(phylo_decontam_rar_abundance_df, 
                                             as.data.frame(phylo_decontam_rar_abundance@sam_data)[,1:2],
                                             by = "SampleID")

sig_SVs <- left_join(sig_SVs, phylo_decontam_rar_abundance_df, by = "SV")
sig_SVs <- left_join(sig_SVs, abundance_means, by = "SV")
sig_SVs$mean_abundance <- as.numeric(sig_SVs$mean_abundance)
view(sig_SVs)
ASV3 <- "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
ASV3 %in% unique(sig_SVs$SV)
view(subset(sig_SVs, SV == ASV3))

my_t_test_sig_SV_sequences <- unique(sig_SVs$SV)

saveRDS(sig_SVs, "sig_SVs.rds")

# Lab 10: Comparing changes in taxonomy (Lab 10) ---------------------------------
## Focusing on a single taxon -------


# transform to relative abundance
relabun.ps <- transform_sample_counts(phylo_decontam_rar, function(x) x / sum(x)) 

# glom ASVs by genus
ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE) 

#subset by the genus of choice. Taxon is not present if you get error: "length of 'dimnames' [1] not equal to array extent" 
ps_genusP <- subset_taxa(ps_genus, Genus %in% c("Streptococcus"))                      

# melt the data into a different configuration
genus.df <- psmelt(ps_genusP) 

# grab that abundance data
MySummary <- genus.df %>% group_by(Group) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 

#check it
head(MySummary)

# graph it
ggplot(data = MySummary, aes(x = Group, y = mean_abund)) +
  geom_point(aes(color = Group)) +
  ylab("Mean relative abundance of reads")



## DESeq Differential Abundance  ---------------------------------
# DESeq only does pairwise comparisons. To make a multifactorial comparison and graph, use "DESeq_and_ternary_plot_example.R".  You can also subset your data

packageVersion("DESeq2")

phylo_decontam_with_species <- readRDS("phylo_decontam_with_species.rds")
# The negative controls are in this cleaned but un-rarefied data set
# take those out
phylo_decontam_with_species_exp_samples <- subset_and_trim(phylo_decontam_with_species, "Sample_type", "experimental")  

# sample_data(phylo_decontam_with_species_exp_samples)$Group <- factor(sample_data(phylo_decontam_with_species_exp_samples)$Group, 
#                                                 levels = c("IPD_ATB", "IPD"), 
#                                                 labels = c("Antibiotics", "No antibiotics")) 
sample_data(phylo_decontam_with_species_exp_samples)$Group <- factor(sample_data(phylo_decontam_with_species_exp_samples)$Group, 
                                                                     levels = c("IPD", "IPD_ATB"), 
                                                                     labels = c("No_antibiotics", "Antibiotics")) 
saveRDS(phylo_decontam_with_species_exp_samples, "phylo_decontam_with_species_exp_samples.rds")


# OPTIONAL if you need to subset 
  # atb_decontam <- subset_and_trim(phylo_decontam_with_species_exp_samples, "Group", "Antibiotics") # subset and trim is in lab5_functions
  # no_atb_decontam <- subset_and_trim(phylo_decontam_with_species_exp_samples, "Group", "No antibiotics")
  # # A little confused; I made these subsets but didn't use them
  # # Maybe this is for if you have more than one factor?


# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(phylo_decontam_with_species_exp_samples, ~ Group)

# calculate differential abundance
gm_mean =  function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = suppressMessages(estimateSizeFactors(diagdds, geoMeans = geoMeans))
diagdds = suppressMessages(DESeq(diagdds, fitType = "local"))

# calculate significance for those abundance calculations
res = suppressMessages(results(diagdds))
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo_decontam_with_species_exp_samples)[rownames(sigtab), ], "matrix")) #CHANGE ME if you didn't subset your data

# head(sigtab)

dim(sigtab) # 17 x 13
# number of SVs that were different between them -- 17


# calculate log changes and set
# sigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ] # If you want only positive (or negative) log changes
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")] #CHANGE ME add Order or Species - if you have it

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels = names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels = names(x))


## if the Genus is empty, replace with the Family NOTE: this used to work but the syntax is broken for unknown reasons
# sigtab.df$Genus[is.na(sigtab.df$Genus)] <- sigtab.df$Family[is.na(sigtab.df$Genus)]

# if the Genus is empty, replace with the Family NOTE: works as of Feb 2020
sigtab$Genus = ifelse(is.na(sigtab$Genus), paste(sigtab$Family),paste(sigtab$Genus)) # ;sigtab

# OPTIONAL bind Genus and Species together - only works if you had species-level taxonomy AND you added it a few steps prior
sigtab$Genus.species <- paste(sigtab$Genus, sigtab$Species) # I ran this but it doesn't change the plot

saveRDS(sigtab, "sigtab.rds")

## graph differential abundance
# in separate plot script




## Feature Prediction (Differential Abundance) with Random Forest ---------------------------------
#  https://rpubs.com/michberr/randomforestmicrobe
# https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
# https://www.rdocumentation.org/packages/randomForest/versions/4.6-12/topics/randomForest


# Make a dataframe of training data with OTUs as column and samples as rows, which is the phyloseq OTU table

predictors <- otu_table(phylo_decontam_with_species_exp_samples, taxa_are_rows = FALSE)

dim(predictors)
# 76 samples, 1884 SVs

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(phylo_decontam_with_species_exp_samples), stringsAsFactors = FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)


# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

# clean up some of the other taxa info
colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))

### start here when choosing factors, can reuse the above lines as needed. one example for factorial data, and one for numeric data is provided. select as needed.

# Make one column for our outcome/response variable. Choose which one applies to the thing you want to test, and then follow factorial or numeric through the rest of the code section.
response <- as.factor(sample_data(phylo_decontam_with_species_exp_samples)$Group) 
# response <- as.numeric(sample_data(EX_ps_clean)$Numeric_Data) #CHANGE ME

# Combine response and SVs into data frame
rf.data <- data.frame(response, predictors)


# set seed for random number generation reproducibility
set.seed(2)

# classify for factorial data
response.pf <- rfPermute(response ~. , data = rf.data, na.action = na.omit, ntree = 500, nrep = 100) #na.omit ignores NAs in data (not tolerated). ntrees is how many forests to build, nreps generates p-value

# or this way for numeric data
# response.pf <- rfPermute(as.numeric(response) ~. , data = rf.data, na.action = na.omit, ntree= 500, nrep = 100)

print(response.pf)
# paste the print out here, especially the OOB error. 1-(Out-of-the-box error) = accuracy of your model
# Antibiotics No antibiotics pct.correct LCI_0.95 UCI_0.95
# Antibiotics             52              2        96.3     87.3     99.5
# No antibiotics          16              6        27.3     10.7     50.2
# Overall                 NA             NA        76.3     65.2     85.3
saveRDS(response.pf, "response_pf.rds")

# grab which features were labeled "important"
imp <- importance(response.pf, scale = TRUE)

# Make a data frame with predictor names and their importance
imp.df <- data.frame(predictors = rownames(imp), imp) 


# For factorial data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, MeanDecreaseAccuracy.pval <= 0.05) 
print(dim(imp.sig))
# 24 x 9


# or For factorial data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$MeanDecreaseAccuracy),]



#create levels to the factor based on SV table
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top so many predictors (more than 50 is a crowded graph)
imp.top <- imp.sort[1:50, ]

# figure out what they are and name them, otherwise will just be named the full SV
otunames <- imp.top$predictors


# need to grab the abundance data out of the otu_table and make a new data.frame, and add the taxa names back in
# grab the column names from the otu table that match those in the forest set
pred.abun.colnum <- which(colnames(rf.data) %in% otunames)

# when you find a match, grad the abudnance data
pred.abun <- rf.data[,sort(c(pred.abun.colnum))]

# make this into a dataframe for manipulation
pred.abun.df <- data.frame(pred.abun, stringsAsFactors = FALSE)

# use the row.names (sample names) from the phyloseq object to name the samples in your forest
row.names(pred.abun.df) <- row.names(sample_data(phylo_decontam_with_species_exp_samples))


# set color palette  
col.bro <- (rainbow(6))

# add white to that color list
col.bro <- append(col.bro, "#ffffff", after = 6)

# add some factors that you can use to make your graph pretty, as many as you want
pred.abun.df$Sample <- row.names(sample_data(phylo_decontam_with_species_exp_samples)) #always grab the sample names
pred.abun.df$Group <- sample_data(phylo_decontam_with_species_exp_samples)$Group 
# pred.abun.df$FactorB <- sample_data(EX_ps_clean)$FactorB #CHANGE ME


# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
# pred.abun.df$Date <- factor(pred.abun.df$Date, levels=c("21_Apr", "12_May", "01_Jun", "22_Jun", "25_Jul"), labels=c("21 Apr", "12 May", "01 Jun", "22 Jun", "25 Jul")) #CHANGE ME

## Troubleshooting ------

# reload these packages in this order, because sometimes the ddply function breaks
library(plyr)
library(dplyr)
# head(pred.abun.df)
# melt and transform the data using ALL the factors you added
m <- melt(pred.abun.df, id.vars = c("Sample", "Group"))
# m$variable <- str_remove_all(m$variable, ".NA")
# m$variable <- str_replace_all(m$variable, ".NA.", "_") # come back and clean up names later
head(m)
saveRDS(m, "m.rds")
# FIXME
esk_m <- m
esk_m$rescale <- log(1 + as.numeric(m$value))
head(esk_m)
m <- ddply(m, .(variable), mutate, rescale = log(1 + as.numeric(value))) #NOTE: may need to change it to as.numeric(value) if it tells you an error about non-binary operators. if you go to graph it and it doesn't recognize 'rescale', it means this piece failed.
# Note to self: ddply and transform is kind of like pivot_longer in the tidyverse
head(m)
identical(esk_m, m)


# adjusted plot for factorial data, recommend using sample ID as 'factor A' (x value)
ggplot(esk_m, aes(as.factor(Sample), variable)) + 
  theme_minimal() + 
  facet_grid(.~Group, space = "free", scales = "free") + #set up graph facets to separate out levels of FactorA
  geom_tile(aes(fill = esk_rescale), color = "gray") + #add the heatmap coloring 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + #use the preset color palette
  labs(fill = "Log abundance") + #rename legend heading
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 6)) +
  ylab('Predictor Taxa') +
  xlab('Sample') 



## LEFSe ---------------------------------
# can be done within R using the lefse package

# can transform phyloseq data to lefse-ready format with the command that you create, phyloseq_to_lefse
# uses code by seahorse001x: https://github.com/seashore001x/Rrumen/blob/master/phyloseq2lefse.R

# require(dplyr)
# require(tibble)
# 
# # this script defines a function to convert phyloseq object into lefse recognized file format. no need to change anything in this bit.
# phyloseq_to_lefse <- function(physeq){
#   # aggregate to genus level
#   ps_lefse <- physeq %>% tax_glom(taxrank = 'Genus', NArm = F)
#   
#   # extract taxonomic data from phyloseq object and then stored in a vector called lefse_tax
#   lefse_tax <- ps_lefse %>% tax_table %>% data.frame(stringsAsFactors=FALSE)
#   lefse_tax <- replace(lefse_tax, is.na(lefse_tax), 'Unclassified')
#   lefse_tax <- lefse_tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% mutate(id = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")) %>% ungroup %>% pull(id)
#   
#   # extract otu abundance matrix from phyloseq object and annotated with tax information
#   lefse_matrix <- otu_table(ps_lefse) %>% data.frame(stringsAsFactors = F) %>% t %>% data.frame
#    
#   # bookmark this is what is throwing the error
# #  colnames(lefse_matrix) <- lefse_tax
#   
# #  row.names(lefse_matrix) <- lefse_tax
#   
#   # extract sample metadata and order sample same in lefse_matrix
#   lefse_sample <- sample_data(ps_lefse)
#   
#   
#   # convert factor in lefse_sample into character in order to combine sample and abundance data
#   lefse_sample_isfactor <- sapply(lefse_sample, is.factor)
#   lefse_sample[,lefse_sample_isfactor] <- lefse_sample[,lefse_sample_isfactor] %>% lapply(as.character)
#   lefse_sample <- lefse_sample %>% data.frame
#   
#   lefse_table <- full_join(rownames_to_column(lefse_sample), rownames_to_column(lefse_matrix), by = ("rowname" = "rowname")) %>% t
#   
#   return(lefse_table)
# }
# 
# 
# EX_clean_for_lefse <- phyloseq_to_lefse(EX_ps_clean)


# If not already in environment, rerun the lines of section starting "feature prediction" up to when colnames(predictors) change

### start here when choosing factors, can reuse the above lines as needed. one example for factorial data, and one for numeric data is provided. select as needed.

# Make one column for our outcome/response variable. Choose which one applies to the thing you want to test, and then follow factorial or numeric through the rest of the code section.
response.class <- as.factor(sample_data(phylo_decontam_with_species_exp_samples)$Group) #CHANGE ME to the main factor you want to test
# response.subclass <- as.factor(sample_data(EX_ps_clean)$Week)

subjects <- as.factor(sample_data(phylo_decontam_with_species_exp_samples)$SampleID) # change this to be ID of the individual if animals/environments have repeated measures sampling

# Combine response and SVs into data frame
# Ex_ps_for_lefse.df <- data.frame(response.class, response.subclass, subjects, predictors)
phylo_lefse_df <- data.frame(response.class, subjects, predictors)


# save this output and upload it to the LEFSe web version on Galaxy: http://huttenhower.sph.harvard.edu/galaxy/
write.table(phylo_lefse_df, file = "phylo_lefse_df.txt", 
            append = TRUE, sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

### Troubleshooting Galaxy -------




#  Lab 11-12: Beta diversity ordinations and stats (Lab 11 and 12) ---------------------------------
# ordinations in phyloseq: https://joey711.github.io/phyloseq/plot_ordination-examples.html

## PCA  -----------------
# currently, phyloseq doesn't run a PCA, but you can manually perform one using this tutorial: https://www.datacamp.com/community/tutorials/pca-analysis-r


# calculate the components
phylo_decontam_rar_pca <- prcomp(otu_table(phylo_decontam_rar, taxa_are_rows = FALSE), center = TRUE, scale = TRUE) #CHANGE ME

# take a look
summary(phylo_decontam_rar_pca) 

# graph it. add ggplot2 text to make it pretty
ggbiplot(phylo_decontam_rar_pca)
ggbiplot(phylo_decontam_rar_pca, varname.size = 0, varname.abbrev = TRUE)



## PCoA in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
# non-euclidian distance
uJ_pcoa <- ordinate(phylo_decontam_rar, #calculate similarities
                       method = "PCoA", #ordination type
                       "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.
saveRDS(uJ_pcoa, "uJ_pcoa.rds")
# simple ordination
plot_ordination(phylo_decontam_rar, uJ_pcoa, type="samples", color="Group") #CHANGE ME but leave type as "samples"

# Detailed plot in the plot script
  


## NMDS in phyloseq -----------------
# euclidian distance
# use phyloseq to calculate ordination for use in the plot
uJ_nmds <- ordinate(phylo_decontam_rar, #calculate similarities
                       method = "NMDS", #ordination type
                       "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

# report the amount of stress at the end of the NMDS calculation, as it is typical to report that in your manuscript. >0.2 is not great, and >0.3 plot is meaningless
# Paste output here: Run 20 stress 0.2627185
saveRDS(uJ_nmds, "uJ_nmds.rds")
# simple ordination
plot_ordination(phylo_decontam_rar, uJ_nmds, type = "samples", color = "Group") 

# fancy plot in plot script




### ----- permanova stats for ordinations --------

### Troubleshooting adonis -------
# Example basic permanova test using Adonis in vegan
  #FIXME
  # Other people in class are having the same error
adonis2(uJ_nmds ~ Group, as(sample_data(phylo_decontam_rar), "data.frame"), permutations = 9999, na.action = na.omit, by = "terms") #CHANGE ME so variable and factor reflect your data
# paste output here

# Example permanova output
#             Df  SumsOfSqs MeanSqs F.Model R2      Pr(>F)   
#   Diet       1    0.6479 0.64792  1.7560 0.07303 0.0047 **
#   Week       1    0.7136 0.71361  1.9340 0.08044 0.0019 **
#   Diet:Week  1    0.4997 0.49967  1.3542 0.05632 0.0495 * 
#   Residuals 19    7.0106 0.36898         0.79021          
#   Total     22    8.8718                 1.00000          
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# interpretation: The F-model/F statistic is a measure of importance, R2 may not be relevant for factorial data


# skipping the following part because my experimental design is not complicated
    # Or, run a more complicated permanova test using Adonis in vegan if you have a complicated experimental design
    # if want to add repeated sampling use this line
    ctrl <- with(as.data.frame(sample_data(EX_ps_clean.rar)), how(blocks = Sheep_ID, nperm = 999))
    adonis2(EX_uJ_nmds ~ Diet * Week * as.numeric(Weight_kg), data = data.frame(sample_data(EX_ps_clean.rar)), # CHANGE ME so 'Factor' is the name of your Factor variables
            permutations = ctrl, #for repeated measures
            # strata = as.data.frame(sample_data(EX_ps_clean.rar))$Blocking_factor #CHANGE ME for replicates/blocks
    ) 
    
    ### NOTE: if adonis is not working and says your data is of the wrong type:
    # 1) check to make sure you don't have NAs in what you are trying use
    # 2) update the data.table package, then reload the vegan package
    # 3) Repeat the ordination calculation by trying the block of code below:
    
    # Code to use Vegan to calculation ordination data and run stats, in case phyloseq and vegan aren't sharing code well.
    # function and tutorial from https://rpubs.com/DKCH2020/587758

veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}
# export data from phyloseq to vegan-compatible object
vegan_otu_table <- veganotu(phylo_decontam_rar)

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
sampledf <- data.frame(sample_data(phylo_decontam_rar))

# run the ordination calculation, change the variable name from EX_uJ to reflect the calculation method you choose (EX_wBC)
uJ <- vegdist(wisconsin(sqrt(vegan_otu_table)), method = "jaccard", binary = TRUE) #CHANGE ME to be the method and weighted (binary = false) or unweighted version (binary = true)

adonis2(uJ ~ Group, 
        sampledf, 
        permutations = 9999, 
        na.action = na.omit,
        by = "terms")

# Output: 

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999

# adonis2(formula = EX_uJ ~ Group, data = as(sample_data(phylo_decontam_rar), "data.frame"), permutations = 9999, na.action = na.omit, by = "terms")
#           Df SumOfSqs     R2      F Pr(>F)    
# Group     1    1.042 0.0326 2.4604  1e-04 ***
# Residual 73   30.909 0.9674                  
# Total    74   31.950 1.0000                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# from class chat
# adonis2(dune ~ Management*A1, data = dune.env, by = "terms")


#### run betadispersion test to see how tight/loose the clusters are ----
# betadisp_EX <- betadisper(EX_uJ_nmds, sampledf$Group) # CHANGE me so variable and factor reflect your data
betadisp <- betadisper(uJ, sampledf$Group)
betadisp
# paste output here
  # Homogeneity of multivariate dispersions
  # 
  # Call: betadisper(d = EX_uJ, group = sampledf$Group)
  # 
  # No. of Positive Eigenvalues: 74
  # No. of Negative Eigenvalues: 0
  # 
  # Average distance to median:
  #   Antibiotics No antibiotics 
  # 0.6530         0.6095 
  # 
  # Eigenvalues for PCoA axes:
  #   (Showing 8 of 74 eigenvalues)
  # PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
  # 2.6850 1.2362 1.1580 1.0571 0.9137 0.8641 0.7208 0.6437 


# run ANOVA to see if clusters overlap or not
anova(betadisp)
# paste output here
  # Analysis of Variance Table
  # 
  # Response: Distances
  # Df   Sum Sq   Mean Sq F value    Pr(>F)    
  # Groups     1 0.029482 0.0294823  14.329 0.0003123 ***
  #   Residuals 73 0.150205 0.0020576                      
  # ---
  #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


permutest(betadisp)
# paste output here
# 
  # Permutation test for homogeneity of multivariate dispersions
  # Permutation: free
  # Number of permutations: 999
  # 
  # Response: Distances
  #             Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
  #  Groups     1 0.029482 0.0294823 14.329    999  0.001 ***
  # Residuals 73 0.150205 0.0020576                         
  # ---
  #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#  correct for multiple comparisons
TukeyHSD(betadisp)
# paste output here
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = distances ~ group, data = df)
# 
# $group
#                             diff         lwr         upr      p adj
# No antibiotics-Antibiotics -0.0435474 -0.06647553 -0.02061927 0.0003123


# Lab 13: Beta diversity components (Lab 13) ---------------------------------

# Component analysis perform ordination calculations on your sample data (metadata or environmental data) and on your sequence table data (SVs). Having missing or NA in the datasets may throw an error. Having outliers may distort the plots. 



## CCA in phyloseq -------------------------------------------
# correspondence analysis, or optionally, constrained correspondence analysis (a.k.a. canonical correspondence analysis),
# usually used with bell-curve or unimodal relationships hypothesis-based testing does not need normally distributed data,
#but if you have a lot of outliers you may want to consider adding a data transformation.

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! - check
# If you need to alter something, see the "beta diversity components" section.

bray_not_na <- phyloseq::distance(physeq = phylo_decontam_rar, method = "bray") #CHANGE ME

cca_ord <- ordinate(
  physeq = phylo_decontam_rar, #CHANGE ME
  method = "CCA",
  distance = bray_not_na,
  formula = ~ Group +  # If you want to see interaction between variables (like Group and Syndrome), change this + to a *
    Condition(Syndrome) 
)

cca_ord
# paste the output here
  # Call: cca(formula = OTU ~ Group + Condition(Syndrome), data = data)
  # 
  # -- Model Summary --
  #   
  #   Inertia Proportion Rank
  # Total         21.36890    1.00000     
  # Conditional    1.60084    0.07491    5
  # Constrained    0.31292    0.01464    1
  # Unconstrained 19.45513    0.91044   68
  # 
  # Inertia is scaled Chi-square
  # 
  # -- Eigenvalues --
  #   
  #   Eigenvalues for constrained axes:
  #   CCA1 
  # 0.31292 
  # 
  # Eigenvalues for unconstrained axes:
  #   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
  # 0.8499 0.8377 0.8179 0.7584 0.7431 0.6975 0.6720 0.6563 
  # (Showing 8 of 68 unconstrained eigenvalues)


# NOTE: do you have CCA eigenvalues for each of the factors you put in? 
# no...
# CCA1 = x axis, and CCA2 = y axis, so if you don't have both of those it will use CA which is the unconstrained axis instead.

cca_ord <- ordinate(
  physeq = phylo_decontam_rar, #CHANGE ME
  method = "CCA",
  distance = bray_not_na,
  formula = ~ Group)
cca_ord
# output
    # Call: cca(formula = OTU ~ Group, data = data)
    # 
    # -- Model Summary --
    #   
    #   Inertia Proportion Rank
    # Total         21.36890    1.00000     
    # Constrained    0.33920    0.01587    1
    # Unconstrained 21.02970    0.98413   73
    # 
    # Inertia is scaled Chi-square
    # 
    # -- Eigenvalues --
    #   
    #   Eigenvalues for constrained axes:
    #   CCA1 
    # 0.3392 
    # 
    # Eigenvalues for unconstrained axes:
    #   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
    # 0.9844 0.9446 0.8419 0.7711 0.7456 0.7144 0.6947 0.6839 
    # (Showing 8 of 73 unconstrained eigenvalues)

# NOTE: do you have CCA eigenvalues for each of the factors you put in? 
# CCA1 = x axis, and CCA2 = y axis, so if you don't have both of those it will use CA which is the unconstrained axis instead.


# anova of whole model
# anova(cca_ord, permu = 1000)

# anova of the factors (terms) you specified
anova(cca_ord, by = "terms", permu = 1000)
# paste the output here   
    # Permutation test for cca under reduced model
    # Terms added sequentially (first to last)
    # Permutation: free
    # Number of permutations: 999
    # 
    # Model: cca(formula = OTU ~ Group, data = data)
    # Df ChiSquare      F Pr(>F)  
    # Group     1    0.3392 1.1774   0.07 .
    # Residual 73   21.0297                
    # ---
    #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1      


# cca plot with FACTOR A (for me, Group) in the model
cca_plot <- plot_ordination(
  physeq = phylo_decontam_rar, 
  ordination = cca_ord, 
  color = "Group",
  axes = c(1,2)) + 
  theme_minimal() +
  geom_point(aes(colour = Group), size = 4) +
  geom_point(colour = "grey90", size = 1.5) 
cca_plot


# above looks kinda crazy...
# below -- I don't have CCA2, where is it supposed to come from?
# answer: the CCA2 must not be significant from my anova
# comment out that line and it will use CA1 instead

# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cca_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CCA1,
                 # yend = CCA2, # I don't have a CCA2
                 x = 0,
                 y = 0,
                 shape = NULL,
                 color = NULL,
                 label = labels)

label_map <- aes(x = 1.3 * CCA1,
                 # y = 1.3 * CCA2,
                 shape = NULL,
                 color = NULL,
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))


# Make a new graphic
cca_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) 

# It looks crazy so it might not be the best method for this data
# Put in paper: CCA was not significant on anything other than the primary variable (group) so that's why I didn't include it
# The constrained ordinations are best for large data sets with many conditions



## RDA in phyloseq -------------------------------------------
# Redundancy analysis for linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = phylo_decontam_rar, method = "bray")

rda_ord <- ordinate(
  physeq = phylo_decontam_rar, #CHANGEME
  method = "RDA",
  distance = bray_not_na,
  formula = ~ Group
)

rda_ord
# paste the output here

# NOTE: do you have RDA eigenvalues for each of the factors you put in? RDA1 = x axis, and RDA2 = y axis, so if you don't have both of those it will use PC which is the unconstrained axis instead.
    # 
    # Call: rda(formula = OTU ~ Group, data = data)
    # 
    # -- Model Summary --
    #   
    #   Inertia Proportion Rank
    # Total         8.405e+06  1.000e+00     
    # Constrained   3.431e+05  4.082e-02    1
    # Unconstrained 8.062e+06  9.592e-01   73
    # 
    # Inertia is variance
    # 
    # -- Eigenvalues --
    #   
    #   Eigenvalues for constrained axes:
    #   RDA1 
    # 343098 
    # 
    # Eigenvalues for unconstrained axes:
    #   PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8 
    # 1942629 1221730  939907  598770  378089  345196  308192  274695 
    # (Showing 8 of 73 unconstrained eigenvalues)

# anova of whole model
anova(rda_ord, permu = 1000)
# paste the output here    
    # Permutation test for rda under reduced model
    # Permutation: free
    # Number of permutations: 999
    # 
    # Model: rda(formula = OTU ~ Group, data = data)
    # Df Variance      F Pr(>F)    
    # Model     1   343098 3.1069  0.001 ***
    #   Residual 73  8061538                  
    # ---
    #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# anova of the factors (terms) you specified
anova(rda_ord, by = "terms", permu = 1000)
# paste the output here              
    # Permutation test for rda under reduced model
    # Terms added sequentially (first to last)
    # Permutation: free
    # Number of permutations: 999
    # 
    # Model: rda(formula = OTU ~ Group, data = data)
    # Df Variance      F Pr(>F)   
    # Group     1   343098 3.1069  0.002 **
    #   Residual 73  8061538                 
    # ---
    #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# rda plot with FACTORA in the model
rda_plot <- plot_ordination(
  physeq = phylo_decontam_rar, #CHANGE ME
  ordination = rda_ord, 
  color = "Group", 
  axes = c(1,2)) + 
  theme_minimal() +
  # aes(shape = as.factor(Week)) + #CHANGE ME
  geom_point(aes(colour = Group), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 



# Now add the environmental variables as arrows
arrowmat <- vegan::scores(rda_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, 
                 yend = RDA2,
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * RDA1, 
                 y = 1.3 * RDA2,
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 



## dbRDA in phyloseq -------------------------------------------
# Constrained Analysis of Principal Coordinates or distance-based RDA, for non-linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = phylo_decontam_rar, method = "bray")

cap_ord <- ordinate(
  physeq = phylo_decontam_rar, #CHANGEME
  method = "CAP", #borrows capscale from the vegan package
  distance = bray_not_na,
  formula = ~ Group #CHANGE ME use for repeated measures
)

cap_ord
# paste the output here

    # Call: capscale(formula = distance ~ Group, data = data)
    # 
    # -- Model Summary --
    #   
    #   Inertia Proportion Rank
    # Total         30.46888                
    # RealTotal     31.80872    1.00000     
    # Constrained    1.11588    0.03508    1
    # Unconstrained 30.69285    0.96492   54
    # Imaginary     -1.33985                
    # 
    # Inertia is squared Bray distance
    # 
    # -- Eigenvalues --
    #   
    #   Eigenvalues for constrained axes:
    #   CAP1 
    # 1.1159 
    # 
    # Eigenvalues for unconstrained axes:
    #   MDS1  MDS2  MDS3  MDS4  MDS5  MDS6  MDS7  MDS8 
    # 4.405 3.961 2.692 1.934 1.784 1.394 1.297 0.960 
    # (Showing 8 of 54 unconstrained eigenvalues)


# NOTE: do you have CAP eigenvalues for each of the factors you put in? CAP1 = x axis, and CAP2 = y axis, so if you don't have both of those it will use MDS which is the unconstrained axis instead.

# anova of whole model
anova(cap_ord, permu = 1000)
# paste the output here    

# anova of the factors (terms) you specified
anova(cap_ord, by = "terms", permu = 1000)
# paste the output here              

    # Permutation test for capscale under reduced model
    # Terms added sequentially (first to last)
    # Permutation: free
    # Number of permutations: 999
    # 
    # Model: capscale(formula = distance ~ Group, data = data)
    # Df SumOfSqs     F Pr(>F)    
    # Group     1   1.1159 2.654  0.001 ***
    #   Residual 73  30.6928                 
    # ---
    #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# CAP plot with FACTORA in the model
cap_plot <- plot_ordination(
  physeq = phylo_decontam_rar, #CHANGE ME
  ordination = cap_ord, 
  axes = c(1,2)) + 
  theme_minimal() +
  # aes(shape = as.factor(Week)) + #CHANGE ME
  geom_point(aes(colour = Group), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 




# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 # yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 # y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 



## variance partitioning --------------
# How much did any one factor contribute to the sample clustering?
# http://www.hiercourse.com/docs/variationPartitioning.html 

library(vegan)

# create a dataframe from your SV table
EX.df <- as.data.frame(otu_table(EX_ps_clean.rar))

# extract your sample data from the phyloseq object
env.df <- as.data.frame(sample_data(EX_ps_clean.rar))


# calculate Principal Coordinates Of Neighborhood Matrix, diversity distance data transformed into rectangular format
EX.pcnm <- pcnm(dist(bray_not_na))

# environmental variables as predictors of community similarity
cap.env <- capscale(EX.df ~ ., data=env.df, distance='bray')

# calculate CCA 
cap.pcnm <- capscale(EX.df ~ ., data=as.data.frame(scores(EX.pcnm)), distance='bray')


# make a model using SV ordination and sample data
mod0.env <- capscale(EX.df ~ 1, data=env.df, distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# make a model using SV ordination and the component scores
mod0.pcnm <- capscale(EX.df ~ 1, data=as.data.frame(scores(EX.pcnm)), distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# select variables in each predictor table
step.env <- ordistep(mod0.env, scope=formula(cap.env))

# check variance inflation factors, higher number = data are redundant to another factor, 1= data are unique.  If factors are redundant (conflated) to each other (for example they basically report the same thing, go back and remove one and mention they were redundant/conflated in your manuscript Methods)
vif.cca(step.env)


step.pcnm <- ordistep(mod0.pcnm, scope=formula(cap.pcnm))

step.env$anova
# paste output here

step.pcnm$anova
# paste output here

EX.pcnm.sub <- scores(EX.pcnm, 
                      choices=c(1,3:16)) #CHANGE ME to include the significant pcnm numbers from previous command. If you want them all, write 1:16. to select only a few of them, write something like 1,3,6,7:9. 

# partition variation among four predictor tables:
EX.var <- varpart(EX.df, 
                  ~ FACTORA, #CHANGE ME
                  ~ FactorB, #CHANGE ME
                  EX.pcnm.sub, data=env.df)

#plot 
par(mfrow=c(1, 2))
showvarparts(4)
plot(EX.var)



EX.var
# paste output here

anova(rda(EX.df  ~ FACTORA + Condition(EX.pcnm.sub), data=env.df)) # add + Condition(env.df$SAMPLE_ID) for repeated measures
# paste output here

























## fin -----------
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
