# This script is to see what happens if I change the analysis order
# I'm trying to find what is different from my analysis and the paper
# Other script order:
# 1 - Filter, trim
# 2 - Learn errors and dechim
# 3 - Decontam
# 4 - Species assignment
# 5 - Rarifaction

# New order:
# 1 - Filter, trim (same)
# 2 - Learn errors and dechim (same)
# 3 - Species assignment
# 4 - Decontam
# 5 - Rarifaction

# Maybe try swapping new-4 and new-5

# I already have a phyloseq object that has followed new 1, 2, and 3 for the decontamination method side quest
# so start at 4

## Functions and package loading ---------------------------

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

phylo_dirty_with_species <- readRDS("phylo_dirty_with_species.rds")

head(phylo_dirty_with_species@sam_data)
# the "group" column was modified in other analyses; make that the same
sample_data(phylo_dirty_with_species)$Group <- factor(sample_data(phylo_dirty_with_species)$Group, 
                                                levels = c("IPD_ATB", "IPD", "controlneg"), 
                                                labels = c("Antibiotics", "No antibiotics", "Negative control")) 
# and just save that so I don't have to do that again
saveRDS(phylo_dirty_with_species, "phylo_dirty_with_species.rds")

head(phylo_dirty_with_species@otu_table, n = c(6L, 6L))
head(row.names(phylo_dirty_with_species@otu_table))
phylo_dirty_with_species # 5144 taxa, 93 samples
summary(make_SV_summary(as.data.frame(phylo_dirty_with_species@otu_table))$SVs[1:76])
# Experimental samples SV stats:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 9.00   50.75  126.00  163.78  250.75  505.00 
summary(make_SV_summary(as.data.frame(phylo_dirty_with_species@otu_table))$SVs[77:93])
# Negative control samples SV stats:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 1.00    7.00   17.00   17.65   22.00   40.00 
summary(rowSums(phylo_dirty_with_species@otu_table)[1:76])
# Experimental samples read stats:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 6123   43858   67171   71260   92426  200291
summary(rowSums(phylo_dirty_with_species@otu_table)[77:93])
# Negative control samples read stats:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.0   255.0   662.0   892.6  1104.0  2947.0 

phylo_dirty_with_desired_species <- phylo_dirty_with_species %>% 
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast")
summary(rowSums(phylo_dirty_with_desired_species@otu_table)[1:76])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2933   41902   65362   68449   91030  198797


## ------ Decontam method -------------



# skipped first part, no quant data 
# The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?
# 

# The second contaminant identification method uses prevalence (presence/absence) across samples as compared to the prevalence in negative controls. Essentially, is it in most of your samples but not most of your controls?

# 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
sample_data(phylo_dirty_with_desired_species)$is.neg <- sample_data(phylo_dirty_with_desired_species)$Group == "Negative control"

# 8. Calculate prevalence of SVs in samples versus controls
contamdf.prev05 <- isContaminant(phylo_dirty_with_desired_species, method = "prevalence", neg = "is.neg", threshold = 0.5)

# 9. Make a table
table(contamdf.prev05$contaminant)

# 10. Look at it
head(which(contamdf.prev05$contaminant))

# 11. get rid of the contaminants 
# remove them from the original phyloseq object, or the one with cleaned out SVs by frequency
phylo_species_then_decontam <- prune_taxa(!contamdf.prev05$contaminant, phylo_dirty_with_desired_species) 

# 12. Check how many are left
phylo_species_then_decontam
# 93 samples and 1884 SVs left
saveRDS(phylo_species_then_decontam, "phylo_species_then_decontam.rds")

phylo_species_then_decontam <- readRDS("phylo_species_then_decontam.rds")


summary(rowSums(phylo_species_then_decontam@otu_table[1:76,]))
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2922   41322   64258   67958   91020  198662 
# very similar to doing it the other way around
summary(make_SV_summary(as.data.frame(phylo_species_then_decontam@otu_table))$SVs[1:76])
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 5.00   19.00   40.00   60.26   84.50  266.00 
dim(phylo_species_then_decontam@otu_table)
# 93, 1884

# Lab 7: Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y



# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
#

tab <- otu_table(phylo_species_then_decontam, taxa_are_rows = FALSE)
class(tab) <- "matrix" # tab <- as.matrix(tab) will do nothing for some reason
## you get a warning here, ignore, this is what we need to have
tab <- t(tab) # transpose observations to rows

r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# optional to save this plot

# take a look at rowsums, or total sequences per sample

# Result from decontam then species:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2922   40793   63951   67733   90951  197735 
  # smallest reads (sequences) in a sample ______ B136 -- 2922 
  # largest in a sample ______ B166 197735
  # number of experimental samples with <5000 ______ 1
  # In the paper they rarified to 12k, I don't understand how I lost so many reads
sort(rowSums(otu_table(phylo_species_then_decontam))) # lab negs have the least; only one sample < 5k, three have < 12k
rowSums(otu_table(phylo_species_then_decontam)) # the last 17 are the lab neg cts
summary((rowSums(otu_table(phylo_species_then_decontam))[1:76]))
# Output:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 2922   41322   64258   67958   91020  198662

dim(phylo_species_then_decontam@otu_table)
# 93, 1884 -- same as decontam then species

# Output from decontam then species:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 5.00   19.00   39.50   59.76   83.50  265.00 
summary(make_SV_summary(as.data.frame(phylo_species_then_decontam@otu_table))$SVs[1:76])
# Output:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 5.00   19.00   40.00   60.26   84.50  266.00 


phylo_species_then_decontam_rar <- rarefy_even_depth(phylo_species_then_decontam, 
                                        sample.size = 5000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                        replace = TRUE, #sampling with or without replacement
                                        trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                        rngseed = 711, 
                                        verbose = TRUE)
# set.seed(711) was used
# 17 lab negative control samples removed and one experimental sample, same as decontam then species
# 383 OTUs removed -- 363 other way around


phylo_species_then_decontam_rar
# 1501 taxa and 75 samples
# other way around: # 1521 taxa, 75 samples


# this helps with plotting later
sample_data(phylo_species_then_decontam_rar)$Group <- factor(sample_data(phylo_species_then_decontam_rar)$Group, 
                                                levels = c("IPD_ATB", "IPD"), 
                                                labels = c("Antibiotics", "No antibiotics")) 

saveRDS(phylo_species_then_decontam_rar, "phylo_species_then_decontam_rar.rds")




### CONCLUSION --------
# really doesn't matter