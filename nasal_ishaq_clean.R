# At first I did the S. Iqshaq method of cleaning
# I was losing too many reads so I switched to the decontam R package
# Here is the first round of code for posterity

## Comparison to see what taxa were removed in this method ---------
phylo <- readRDS("phylo.rds")

seqtab_nochim = as(otu_table(phylo), "matrix")
seqtab_nochim = as.matrix(seqtab_nochim)

## 2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops.
dirty_all_taxa <- assignTaxonomy(seqtab_nochim, 
                                 'silva_nr99_v138.1_train_set.fa.gz',        # CHANGE file path
                                 tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                 minBoot = 75, verbose = TRUE) ; beep("treasure") 



saveRDS(dirty_all_taxa, 'dirty_all_taxa.rds')
write.csv(dirty_all_taxa, 'dirty_all_taxa.csv')


dirty_all_taxa_species <- addSpecies(dirty_all_taxa, 
                                        'silva_species_assignment_v138.1.fa.gz', 
                                        allowMultiple = FALSE, 
                                        verbose = FALSE) ; beep("treasure") 

saveRDS(dirty_all_taxa_species, 'dirty_all_taxa_species.rds') 
write.csv(dirty_all_taxa_species, 'dirty_all_taxa_species.csv')

meta <- readRDS("meta.rds")
# reload taxa table as needed
dirty_all_taxa_species <- readRDS('dirty_all_taxa_species.rds')

otu_t <- otu_table(phylo)

## create a phyloseq object with all samples
phylo_dirty_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                        sample_data(meta),
                                        tax_table(dirty_all_taxa_species))


saveRDS(phylo_dirty_with_species, "phylo_dirty_with_species.rds")


# Dr. Ishaq's method: https://github.com/SueIshaq/Examples-DADA2-Phyloseq

## ------ Dr Ishaq's method -------------

# Dr. Ishaq's method creates vectors out of the SV table data for negative controls, and subtracts those SVs from the sample data.  Depending on the type of negative control, these are removed from the whole data set or from subsets of batches. Remove PCR and sampling materials negative control SVs fully from all samples, and remove extraction kit SVs fully from each dna_extraction_batch, respectively.
# With modifications by Emily

clean_data <-  clean_phylo_data(phylo, neg_con = "negative")
clean_data
# Ends up with 4952/5144 taxa and 76/93 samples maintained
saveRDS(clean_data, "clean_data.rds")
# Loosen cutoff
relaxed_clean_data <- relaxed_clean_phylo_data(phylo, neg_con = "negative", threshold = 5)
# only added in 2 more taxa, removed 15 out of 17 negative controls
# probably not worth chasing unless those two taxa are super important for some reason
summary(rowSums(relaxed_clean_data@otu_table)) # median and min did not change at all

# Did it work? Check your ordination again
clean_ord <- ordinate(clean_data, #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(clean_data, clean_ord, 
                type = "samples", color = "Group", 
                title = "Ordination plot, after Ishaq clean")
ggsave("after_clean_ordination_plot.png") # save this graph for later

summary(rowSums(clean_data@otu_table))
comparison[nrow(comparison) + 1, ] <- c("Median reads per exp_samp after clean", 19226, 47323)
comparison[nrow(comparison) + 1, ] <- c("Minimum reads per exp_samp after clean", 287, NA)

comparison[nrow(comparison) + 1, ] <- c("Total SVs before clean", 5144, NA)
comparison[nrow(comparison) + 1, ] <- c("Total SVs after clean", 4952, NA)
summary(make_SV_summary(as.data.frame(clean_data@otu_table))$SVs)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.0    42.0   113.5   150.7   233.0   489.0 


### Summary of decontamination ----
# FIXME
# CHANGE ME: Which types of negative controls did you use and why? 
# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?


## ------ Decontam method -------------
# Analysis branched off, see nasal_decontam_method.R

### Summary of decontamination ----
# CHANGE ME: Which types of negative controls did you use and why? 
# CHANGE ME: Did you notice clustering by batch before and after removing contaminants?








# Lab 6_A: DADA2 assigning taxonomy from a reference database file (Lab 6) ---------------------------------


# Be sure to have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 

# The Silva file shown below is used for 16S rRNA (prokaryotes) and nicely formatted versions can be downloaded from the DADA2 website, which also contains some (18S) eukaryotic versions.  Let me know if you need a custom one for your dataset.
#
# 1. take the SVs (otu table) from the phyloseq object and make it into a dataframe again


clean_data <- readRDS("clean_data.rds")

seqtab_nochim_decontam = as(otu_table(clean_data), "matrix")
seqtab_nochim_decontam = as.matrix(seqtab_nochim_decontam)

## 2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops.
clean_all_taxa <- assignTaxonomy(seqtab_nochim_decontam, 
                                 'silva_nr99_v138.1_train_set.fa.gz',        # CHANGE file path
                                 tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                 minBoot = 75, verbose = TRUE) # ; beep("treasure") 



saveRDS(clean_all_taxa, 'clean_all_taxa.rds')
write.csv(clean_all_taxa, 'clean_all_taxa.csv')
# clean_all.taxa

# OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.
# clean_all_taxa <- readRDS("clean_all_taxa.rds")

# clean_all_taxa <- readRDS("clean_all_taxa.rds")
clean_all_taxa_species <- addSpecies(clean_all_taxa, 'silva_species_assignment_v138.1.fa.gz', allowMultiple = FALSE, verbose = FALSE) # ; beep("treasure") 

saveRDS(clean_all_taxa_species, 'clean_all_taxa_species.rds') 
write.csv(clean_all_taxa_species, 'clean_all_taxa_species.csv') 

## Remake the phyloseq object with the new taxonomy file ----------------------

# reload metadata table as needed
meta <- readRDS("meta.rds")
clean_data <- readRDS("clean_data.rds")

#reload taxa table as needed
clean_all_taxa_species <- readRDS('clean_all_taxa_species.rds')

# troubleshooting
class(otu_table(clean_data))
class(meta)
class(tax_table(clean_all_taxa_species))
class(clean_all_taxa_species)

otu_t <- otu_table(clean_data)

## create a phyloseq object with all samples
phylo_clean_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                     sample_data(meta),
                                     tax_table(clean_all_taxa_species))


saveRDS(phylo_clean_with_species, "Ishaq_phylo_clean_with_species.rds")



# Lab 6_B: Clean out unwanted taxa ------------------------------------
# Regardless of whether you included sequenced negative controls, you can remove taxa which are of no interest to you

# clean out chloroplast and mitochondria. can also elect to remove other contaminating domains or kingdoms as needed. 

# phylo_clean_with_species <- readRDS("phylo_clean_with_species.rds")
# Optional: explore your taxonomy before filtering. Use the tax table you made

df <- as.data.frame(phylo_clean_with_species@tax_table)
table(df$Kingdom)
table(df$Phylum)
table(df$Class)
table(df$Order)
table(df$Family)

phylo_clean_with_species <- phylo_clean_with_species %>% # CHANGE ME to your phyloseq object name 
  subset_taxa(Family != "Mitochondria" & Order != "Chloroplast") # CHANGE ME to taxa you want to remove

saveRDS(phylo_clean_with_species, "phylo_clean_with_species.rds")

# Write out a description of experimental design (Homework)
# FIXME


# Rarefaction (Lab 7_A) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y



# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
#

tab <- otu_table(phylo_clean_with_species)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
# rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)


r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_clean_with_species)))
# smallest reads (sequences) in a sample ______ B271 - 122
# largest in a sample ______ B166- 186343
# number of samples with <5000 ______ 17
# If I cut at 1k I would only lose 2 samples
# Cut off at 2000 which sadly does mean 12 are lost
# In the paper they rarified to 12k, I don't understand how I lost so many reads
comparison[nrow(comparison) + 1, ] <- c("Median reads before rarify", 15255, NA)
comparison[nrow(comparison) + 1, ] <- c("Minimum reads before rarify", 122, 12000)


phylo_rar <- rarefy_even_depth(phylo_clean_with_species, 
                               sample.size = 2000, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                               replace = TRUE, #sampling with or without replacement
                               trimOTUs = TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                               rngseed = 711, 
                               verbose = TRUE)
# set.seed(711) was used
# 12 samples removed
# 358 OTUs removed


# this helps with plotting later
sample_data(phylo_rar)$Group <- factor(sample_data(phylo_rar)$Group, 
                                       levels = c("IPD_ATB", "IPD"), 
                                       labels = c("Antibiotics", "No antibiotics")) #CHANGE ME

saveRDS(phylo_rar, "clean_phylo_rarified.rds")

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(phylo_rar), 'clean_phylo_rarified.csv')

# Specifying taxa after rarification

phylo_clean_strep_rar  <- phylo_rar  %>% 
  subset_taxa(Family == "Streptococcaceae") 

phylo_clean_no_strep_rar   <- phylo_rar  %>% 
  subset_taxa(Family != "Streptococcaceae") 


## clean out taxa/SV columns that are no longer present
phylo_clean_strep_rar <- prune_taxa(taxa_sums(phylo_clean_strep_rar) > 0, phylo_clean_strep_rar)
phylo_clean_strep_rar <- prune_samples(sample_sums(phylo_clean_strep_rar) > 0, phylo_clean_strep_rar)
phylo_clean_strep_rar
# 51 samples, 54 taxa

phylo_clean_no_strep_rar <- prune_taxa(taxa_sums(phylo_clean_no_strep_rar) > 0, phylo_clean_no_strep_rar)
phylo_clean_no_strep_rar <- prune_samples(sample_sums(phylo_clean_no_strep_rar) > 0, phylo_clean_no_strep_rar)
phylo_clean_no_strep_rar
# 64 samples, 1360 taxa


saveRDS(phylo_clean_strep_rar, 'phylo_clean_strep_rar.rds') 
saveRDS(phylo_clean_no_strep_rar, 'phylo_clean_no_strep_rar.rds') 


# meta <- readRDS("meta.rds")
# phylo_rar <- readRDS("clean_phylo_rarified.rds")

# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

### Alpha diversity plotting ------------------

# This script was getting too long and busy so I have split the plotting
# See separate .R file named nasal_aspirate_plots



