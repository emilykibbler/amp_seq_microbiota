

## ------ Decontam method -------------
phylo <- readRDS("phylo.rds")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam) # This packages identify contaminants by frequency of SVs.

#### skipped, no quant data ----------
# The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?
# 
# # 1. Look for contaminants using DNA quantification data from your metadata file
# contamdf.freq <- isContaminant(EX_ps, method="frequency", conc="quant_reading") #CHANGE conc= "quant_reading" to the column heading that holds the concentration information in your metadata
# 
# # 2. See what it looks like 
# head(contamdf.freq) 
# 
# # 3. Make it into a table
# table(contamdf.freq$contaminant)
# 
# # 4. See which SVs are categorized as contaminants based on frequency
# head(which(contamdf.freq$contaminant)) 
# 
# # 5. Get rid of the contaminants 
# EX_ps.noncontamfreq <- prune_taxa(!contamdf.freq$contaminant, EX_ps)
# 
# # 6. How much is left
# EX_ps.noncontamfreg
# # CHANGE ME: number of samples and SVs left 



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

# clean_all_taxa <- readRDS("clean_all_taxa.rds")
decontam_all_taxa_species <- addSpecies(decontam_all_taxa, 'silva_species_assignment_v138.1.fa.gz', allowMultiple = FALSE, verbose = FALSE) # ; beep("treasure") 

saveRDS(decontam_all_taxa_species, 'decontam_all_taxa_species.rds') 
write.csv(decontam_all_taxa_species, 'decontam_all_taxa_species.csv') 

## Remake the phyloseq object with the new taxonomy file ----------------------

# reload metadata table as needed
meta <- readRDS("meta.rds")
phylo_clean_decontam_decontam <- readRDS("phylo_clean_decontam_decontam.rds")

#reload taxa table as needed
decontam_all_taxa_species <- readRDS('decontam_all_taxa_species.rds')

# troubleshooting
class(otu_table(decontam_data))
class(meta)
class(tax_table(decontam_all_taxa_species))
class(decontam_all_taxa_species)

otu_t <- otu_table(decontam_data)

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


# Rarefaction (Lab 7_A) --------------------------------------------------
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
# smallest reads (sequences) in a sample ______ #FIXME
# largest in a sample ______ 
# number of samples with <5000 ______ 
# If I cut at 1k I would only lose 2 samples
# Cut off at 2000 which sadly does mean 12 are lost
# In the paper they rarified to 12k, I don't understand how I lost so many reads
# comparison[nrow(comparison) + 1, ] <- c("Median reads before rarify", 15255, NA)
# comparison[nrow(comparison) + 1, ] <- c("Minimum reads before rarify", 122, 12000)


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






## Alpha diversity plotted against other metadata -----

# use phyloseq to measure alpha diversity
phylo_rar_rich <- estimate_richness(phylo_rar, measures = c("Observed", "Shannon")) #change to whatever measures you want

# # OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# install.packages("remotes")
# remotes::install_github("twbattaglia/btools")
# EX_faith <- estimate_pd(EX_ps_clean.rar)

# measure evenness for each sample
phylo_rar_even <- phylo_rar_rich$Shannon/log(phylo_rar_rich$Observed)

# Coerce to data.frame and add the metadata for these samples
phylo_rar_sd = as(sample_data(phylo_rar), "matrix")
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

# skipped, maybe come back later

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
