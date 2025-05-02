## Short intro -----

# This analysis is on the data from this paper:

# Crovetto, F., Selma-Royo, M., Crispi, F. et al. 
# Nasopharyngeal microbiota profiling of pregnant women with SARS-CoV-2 infection. 
# Sci Rep 12, 13404 (2022). https://doi.org/10.1038/s41598-022-17542-z

# 16S amplicon sequencing reads have generously been made available in the SRA
# These are used in this pipeline for analyzing microbial population composition

# This code is based on a previous semester project;
# However, the focus here is comparing this previously-developed pipeline to a new (to me) pipeline:
# Qiime2 pipeline, which I ran in a combination of interactive and SLURM scripts
# on the Ohio Supercomputer platform
# The command line code will also be available at github.com/emilykibbler/amp_seq_microbiota



## Start ----------
setwd("/Users/emilykibbler/Desktop/projects/bms690")

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/AVS554_packages.R")
# install_necessary()
# install_optional()

load_libraries()

## Preprocessing ----------
queso <- list.files("/Users/emilykibbler/Desktop/projects/bms690/raw", full.names = T)
read1fns <- queso[grep("_1.fastq", queso)]
read2fns <- queso[grep("_2.fastq", queso)]
# 
dir.create("/Users/emilykibbler/Desktop/projects/bms690/raw/fwd_raw")
dir.create("/Users/emilykibbler/Desktop/projects/bms690/raw/rev_raw")
file.copy(read1fns, "/Users/emilykibbler/Desktop/projects/bms690/raw/fwd_raw")
file.copy(read2fns, "/Users/emilykibbler/Desktop/projects/bms690/raw/rev_raw")

## Initial quality inspection, filter and trim ----------
file_info <- metadata_import(types = c("fwd", "rev")) # saves a copy as an rds
qualplots(file_info)

filtoutput <- filterAndTrim( 
  file_info$file_names_full[[1]], 
  file.path(file_info$filt_dir[1], paste0(file_info$sample_names[[1]], "_F_filt.fastq.gz")), 
  file_info$file_names_full[[2]], 
  file.path(file_info$filt_dir[2], paste0(file_info$sample_names[[1]], "_R_filt.fastq.gz")), 
  trimLeft = c(10, 10), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina
  truncLen = c(284, 224), 
  maxEE = c(5,5), # max errors tolerated
  verbose = FALSE)
saveRDS(filtoutput, "filtoutput.rds")

qualplots(file_info, type = "filtered")


## DADA2 learn error rates  ----------------------------


file_info <- readRDS("file_info.rds")
create_and_plot_error_profile(file_info, bases = 1e6) ; beep("treasure") # writes the error profile RDS files
# First time:
  # Fwd: 14274480 total bases in 59477 reads from 1 samples will be used for learning the error rates
  # Rev: 15115440 total bases in 62981 reads from 1 samples will be used for learning the error rates.
# Second time:
  # 17072666 total bases in 62309 reads from 1 samples will be used for learning the error rates.
  # 13792086 total bases in 64449 reads from 1 samples will be used for learning the error rates.

errF <- readRDS("error_profile_fwd.rds")
errR <- readRDS("error_profile_rev.rds")

dadaFs <- dada(list.files(file_info$filt_dir[1], full.names = TRUE),
               err = errF, 
               multithread = TRUE, verbose = FALSE) #; beep()
dadaRs <- dada(list.files(file_info$filt_dir[2], full.names = TRUE),
               err = errR, 
               multithread = TRUE, verbose = FALSE); beep()


## DADA2 pick sequence variants  ----------------------


mergers <- mergePairs(dadaF = dadaFs,
                      derepF = list.files(file_info$filt_dir[1], full.names = TRUE),
                      dadaR = dadaRs,
                      derepR = list.files(file_info$filt_dir[2], full.names = TRUE),
                      verbose = FALSE) 
seqtab <- makeSequenceTable(mergers)
rowSums(seqtab)
view(seqtab)
saveRDS(seqtab, "seqtab.rds"); beep("treasure")
# seqtab <- readRDS("seqtab.rds")

SVs_found_by_sample <- make_SV_summary(as.data.frame(seqtab)) # this function takes so long, look into that
summary(SVs_found_by_sample$SVs)
# min is 113, median is 365, max is 364


## Remove chimeras.-------------
# Leave the method as consensus. 
# multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab_nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", 
                                    multithread = TRUE,
                                    verbose = TRUE) 
# First time:
  # My samples have Identified 8341 bimeras out of 21279 input sequences.
# Second time:
 # 6131 bimeras out of 18202 input sequences.
# View(make_SV_summary(seqtab_nochim))
summary(make_SV_summary(as.data.frame(seqtab_nochim))$SVs)
# First time:
  # min 105, median 298, max 472
# Second time:
  # min 76, median 257, max 433
# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab_nochim) 
# 76 samples (row count) and 12938 SVs (column count)

# Calculate the percentage of chimeras identified out of the total
round(sum(seqtab_nochim)/sum(seqtab), digits = 3) # 0.886 first time, 0.844 second time
sort(rowSums(seqtab_nochim))
# first time:
  # Min reads is 29105, max is 83067
# second time:
  # min 28180, max 79055

saveRDS(seqtab_nochim, 'seqtab_nochim.rds')


# 2. Load the metadata file that goes with your sequencing data so you can match factors to seq data
meta <- read.table("/Users/emilykibbler/Documents/Classes/BMS_690/cov-project/metadata.tsv", sep = "\t", header = T)
head(meta)
row.names(meta) <- meta$sample.id
saveRDS(meta, "meta.rds")


# Lab 5 :Phyloseq First look (Lab 5_A) ------------------------------

# # load data if necessary
seqtab_nochim <- readRDS('seqtab_nochim.rds')

row.names(seqtab_nochim) <- row.names(meta)
# save it with the new row names
saveRDS(seqtab_nochim, 'seqtab_nochim.rds')

# Check the sample sames and see if they match
identical(row.names(meta), row.names(seqtab_nochim))


## create a phyloseq object with all samples (controls will be subset out later)
phylo <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = FALSE), # even though it's called an OTU table, it will use the SVs from my seqtab
                  sample_data(meta))

phylo # 12071 ASVs
saveRDS(phylo, "phylo.rds")
phylo <- readRDS("phylo.rds")



# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
plot(sort(taxa_sums(phylo), TRUE), 
     type = "h", 
     ylim = c(0, 20)) #limit the y-axis to better see the long tail of rare taxa


# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
phylo_ord <- ordinate(phylo, #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

all_taxa <- assignTaxonomy(seqtab_nochim, 
                                    '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_nr99_v138.1_train_set.fa.gz',
                                    tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                    minBoot = 75, verbose = TRUE); beep("treasure") 
saveRDS(all_taxa, 'all_taxa.rds')
write.csv(all_taxa, 'all_taxa.csv')

# all_taxa_species <- addSpecies(all_taxa, 
#                               '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_species_assignment_v138.1.fa.gz', 
#                               allowMultiple = FALSE, 
#                               verbose = FALSE) ; beep("treasure") 
# too big to do it all at once

all_taxa_species <- addSpecies(all_taxa[1:2,], 
                               '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_species_assignment_v138.1.fa.gz', 
                                     allowMultiple = FALSE, 
                                     verbose = FALSE)

all_taxa_species <- esk_add_species(all_taxa, 3:5000, all_taxa_species, fp = '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_species_assignment_v138.1.fa.gz')
all_taxa_species <- esk_add_species(all_taxa, 5001:10000, all_taxa_species, fp = '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_species_assignment_v138.1.fa.gz')#; beep("treasure")
all_taxa_species <- esk_add_species(all_taxa, 10001:nrow(all_taxa), all_taxa_species, fp = '/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/silva_species_assignment_v138.1.fa.gz'); beep("treasure")
write_rds(all_taxa_species, "all_taxa_species.rds")

#reload data as needed
all_taxa_species <- readRDS('all_taxa_species.rds')
meta <- readRDS("meta.rds")
phylo <- readRDS("phylo.rds")

otu_t <- otu_table(phylo, taxa_are_rows = FALSE)

## create a phyloseq object with all samples
phylo_with_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                        sample_data(meta),
                                        tax_table(all_taxa_species))


saveRDS(phylo_with_species, "phylo_with_species.rds")
summary(rowSums(phylo_with_species@otu_table))

df <- as.data.frame(phylo_with_species@tax_table)
table(df$Kingdom) # 5 eukaryotes ASVs
table(df$Phylum)
table(df$Class)
table(df$Order) # 33 chloroplast
table(df$Family) # 879 mitochondria 

phylo_relevant_species <- phylo_with_species %>% 
  subset_taxa(Kingdom != "Eukaryota" & Family != "Mitochondria" & Order != "Chloroplast")

saveRDS(phylo_relevant_species, "phylo_relevant_species.rds")


tab <- otu_table(phylo_relevant_species, taxa_are_rows = FALSE)
class(tab) <- "matrix"
## you get a warning here, ignore, this is what we need to have
tab <- t(tab) # transpose observations to rows

r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(phylo_relevant_species)))
summary((rowSums(otu_table(phylo_relevant_species))))
# at this point: min 7043
sort(rowSums(otu_table(phylo_with_species)))
summary((rowSums(otu_table(phylo_with_species))))
# woah min with all species is 28180, irrelevant is 3/4

# but that is closer to what they reported in the paper

phylo_rar <- rarefy_even_depth(phylo_relevant_species, 
                                        sample.size = 7043, # rarefaction depth
                                        replace = TRUE, #sampling with or without replacement
                                        trimOTUs = TRUE, #remove SVs left empty
                                        rngseed = 711, 
                                        verbose = TRUE)
# took out 310 OTUs and no samples, since I rarefied to min sample depth
saveRDS(phylo_rar, "phylo_rar.rds")



## Alpha diversity plotted against other metadata --------

# use phyloseq to measure alpha diversity
# Observed and Shannon were in the publication
phylo_rar_rich <- estimate_richness(phylo_rar, 
                                    measures = c("Observed", "Shannon")) #change to whatever measures you want

# # OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
# install.packages("remotes")
# remotes::install_github("twbattaglia/btools")
# install.packages("picante")
# library(picante)
r_faith <- pd(comm, phylo_rar)
remotes::install_github("RIVM-IIV-Microbiome/biomeUtils")
library(biomeUtils)
r_faith <- calculatePD(phylo_rar, include_root = FALSE)
remotes::install_github("leylabmpi/LeyLabRMisc")
install.packages("renv")
library(renv)
renv::install("leylabmpi/LeyLabRMisc")


# measure evenness for each sample
phylo_rar_even <- phylo_rar_rich$Shannon/log(phylo_rar_rich$Observed)

# Coerce to data.frame and add the metadata for these samples
phylo_rar_sd <- as(sample_data(phylo_rar), "matrix")
phylo_rar_sd <- as.data.frame(phylo_rar_sd)
phylo_rar_rich_df <- cbind(phylo_rar_rich, phylo_rar_even, phylo_rar_sd)

dim(phylo_rar_rich_df) # 76 x 5
head(phylo_rar_rich_df)

## Changing order of things ----------

summary((rowSums(otu_table(phylo))))
# Min 2810, median 43614, max 79055

dat <- as.data.frame(otu_table(phylo, taxa_are_rows = FALSE))
view(dat)

tab <- otu_table(phylo, taxa_are_rows = FALSE)
class(tab) <- "matrix" # you get a warning here, ignore, this is what we need to have
tab <- t(tab) # transpose observations to rows
r_curve <- rarecurve(tab, step = 10, cex = 0.5, label = FALSE) 
# I don't really get this, what is the Y axis? It's not SVs, I have way more of those

phylo_no_spec_rar <- rarefy_even_depth(phylo, 
                               sample.size = 28000, # rarefaction depth: just above the minimum sample depth (but a nice even number)
                               replace = TRUE, #sampling with or without replacement
                               trimOTUs = TRUE, #remove SVs left empty
                               rngseed = 711, 
                               verbose = TRUE)
phylo_no_spec_rar
# 46 OTUs removed

saveRDS(phylo_no_spec_rar, "phylo_no_spec_rar.rds")


# Add species in now
view(seqtab_nochim) # seqtab_nochim is what I used before as the data for this step
# sample names are row names, ASV sequences are column names, and reads are in the cells
class(seqtab_nochim) # matrix array
dim(seqtab_nochim) # 76 by 12071

view(phylo_no_spec_rar@otu_table)
# I think that is the equivalent table
class(phylo_no_spec_rar@otu_table) # phyloseq object
dim(phylo_no_spec_rar@otu_table) # 76 by 12025

all_taxa_rar <- assignTaxonomy(as.matrix(phylo_no_spec_rar@otu_table), 
                           'silva_nr99_v138.2_toGenus_trainset.fa.gz',
                           tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                           minBoot = 75, verbose = TRUE); beep("treasure") 
saveRDS(all_taxa_rar, 'all_taxa_rar.rds')
write.csv(all_taxa_rar, 'all_taxa_rar.csv')


dim(all_taxa_rar)

all_taxa_species_rar <- addSpecies(all_taxa_rar[1:2,], 
                               'silva_v138.2_assignSpecies.fa.gz', 
                               allowMultiple = FALSE, 
                               verbose = FALSE)

all_taxa_species_rar <- esk_add_species(all_taxa_rar, 3:5000, all_taxa_species_rar, fp = 'silva_v138.2_assignSpecies.fa.gz')
all_taxa_species_rar <- esk_add_species(all_taxa_rar, 5001:10000, all_taxa_species_rar, fp = 'silva_v138.2_assignSpecies.fa.gz')#; beep("treasure")
all_taxa_species_rar <- esk_add_species(all_taxa_rar, 10001:nrow(all_taxa_rar), all_taxa_species_rar, fp = 'silva_v138.2_assignSpecies.fa.gz'); beep("treasure")
write_rds(all_taxa_species_rar, "all_taxa_species_rars.rds")

# Put that together with the phyloseq object
otu_t <- otu_table(phylo_rar, taxa_are_rows = FALSE)
## create a phyloseq object with all samples
phylo_rar_then_species <- phyloseq(otu_table(otu_t, taxa_are_rows = FALSE), 
                                        sample_data(meta),
                                        tax_table(all_taxa_species_rar))


##Alpha diversity analysis, imported from Qiime analysis ----------
faith <- read.table("faith_group_signif.tsv", sep = "\t", header = TRUE)
faith <- rename(faith, "Value" = faith_pd)
# wilcox.test(subset(faith, SarsCov2 == "pos")$Value, subset(faith, SarsCov2 == "neg")$Value)
t.test(subset(faith, SarsCov2 == "pos")$Value, subset(faith, SarsCov2 == "neg")$Value)

evenness <- read.table("even_group_signif.tsv", sep = "\t", header = TRUE)
evenness <- rename(evenness, "Value" = pielou_evenness)
t.test(subset(evenness, SarsCov2 == "pos")$Value, subset(evenness, SarsCov2 == "neg")$Value)
# wilcox.test(subset(evenness, SarsCov2 == "pos")$Value, subset(evenness, SarsCov2 == "neg")$Value)

shannon <- read.table("shannon_group_signif.tsv", sep = "\t", header = TRUE)
shannon <- rename(shannon, "Value" = shannon_entropy)
t.test(subset(shannon, SarsCov2 == "pos")$Value, subset(shannon, SarsCov2 == "neg")$Value)
# wilcox.test(subset(shannon, SarsCov2 == "pos")$Value, subset(shannon, SarsCov2 == "neg")$Value)

obs <- read.table("observed_stats.tsv", sep = "\t", header = T)
obs <- rename(obs, "Value" = observed_features)
obs$Metric <- "Observed"
t.test(subset(obs, SarsCov2 == "pos")$Value, subset(obs, SarsCov2 == "neg")$Value)

## DADA analysis, alpha div------------------

# Rare, then species
diversity <- estimate_richness(phylo_no_spec_rar, measures = c("Shannon", "Observed"))
even <- evenness(phylo_no_spec_rar, index = "pielou")
diversity <- cbind(diversity, even)
diversity$sample.id <- row.names(diversity)
diversity <- merge(diversity, meta, by = "sample.id")
t.test(subset(diversity, SarsCov2 == "pos")$Observed, subset(diversity, SarsCov2 == "neg")$Observed)
t.test(subset(diversity, SarsCov2 == "pos")$Shannon, subset(diversity, SarsCov2 == "neg")$Shannon)
t.test(subset(diversity, SarsCov2 == "pos")$pielou, subset(diversity, SarsCov2 == "neg")$pielou)

# Species, filter, then rar
diversity <- estimate_richness(phylo_rar, measures = c("Shannon", "Observed"))
even <- evenness(phylo_rar, index = "pielou")
diversity <- cbind(diversity, even)
diversity$sample.id <- row.names(diversity)
diversity <- merge(diversity, meta, by = "sample.id")
t.test(subset(diversity, SarsCov2 == "pos")$Observed, subset(diversity, SarsCov2 == "neg")$Observed)
t.test(subset(diversity, SarsCov2 == "pos")$Shannon, subset(diversity, SarsCov2 == "neg")$Shannon)
t.test(subset(diversity, SarsCov2 == "pos")$pielou, subset(diversity, SarsCov2 == "neg")$pielou)




## DESeq-----------------------

# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(phylo_rar_then_species, ~ SarsCov2)

# calculate differential abundance
gm_mean =  function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = suppressMessages(estimateSizeFactors(diagdds, geoMeans = geoMeans))
diagdds = suppressMessages(DESeq(diagdds, fitType = "local"))

# calculate significance for those abundance calculations
res <- suppressMessages(results(diagdds))
res <- res[order(res$padj, na.last = NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(phylo_rar_then_species)[rownames(sigtab), ], 
                  "matrix")) #CHANGE ME if you didn't subset your data

head(sigtab)
dim(sigtab) 


# calculate log changes and set
# sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus")] #CHANGE ME add Order or Species - if you have it
sigtab <- subset(sigtab, select = -Species) # nothing in here got species-level assignment so just drop that
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels = names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels = names(x))


## if the Genus is empty, replace with the Family 
sigtab$Genus = ifelse(is.na(sigtab$Genus), 
                      paste(sigtab$Family), 
                      paste(sigtab$Genus)) 
view(sigtab)

DEseq_sig_SVs <- row.names(sigtab)

## Qiime DA ----------------------

da <- read.csv("qiime_da.csv")
head(da)
da <- subset(da, select = -X)
tax <- read.table("silva_tax_assignments.tsv", sep = "\t", header = TRUE)
head(tax)
tax <- rename(tax, "id" = Feature.ID)
da <- merge(da, tax, by = "id")

# da$kingdom <- str_remove_all(str_split_i(da$Taxon, ";", 1), "d__")
# da$Phylum <- str_remove_all(str_split_i(da$Taxon, ";", 2), "p__")
# da$Class <- str_remove_all(str_split_i(da$Taxon, ";", 3), "c__")
# da$Order <- str_remove_all(str_split_i(da$Taxon, ";", 4), "o__")
# da$Family <- str_remove_all(str_split_i(da$Taxon, ";", 5), "f__")
# da$Genus <- str_remove_all(str_split_i(da$Taxon, ";", 6), "g__")
# da$Species <- str_remove_all(str_split_i(da$Taxon, ";", 7), "s__")

da <- split_taxon_column(da)

head(da)
da <- subset(da, select = -Species)

da$Genus = ifelse(is.na(da$Genus), 
                      paste(da$Family), 
                      paste(da$Genus)) 

da$Genus = ifelse(da$Genus == "NA", 
                  paste(da$kingdom), 
                  paste(da$Genus)) 
da[5,"Genus"] <- "Mitochondria_A"
da[6,"Genus"] <- "Mitochondria_B"
da[9,"Genus"] <- "Mitochondria_C"

da[2,"Genus"] <- "Holdemanella_A"
da[8,"Genus"] <- "Holdemanella_B"


# write_rds(da, "da.rds")

# head(read.table("descriptive_stats.tsv", sep = "\t", header = TRUE))

seqs <- readDNAStringSet("sequences.fasta", format = "fasta")
seqs <- as.data.frame(seqs)
seqs$id <- row.names(seqs)
seqs <- rename(seqs, "Sequence" = x)

da <- merge(da, seqs, by = "id")

write_rds(da, "da.rds")

match(da$Sequence, DEseq_sig_SVs)
view(da[c(2,4),])

match(DEseq_sig_SVs, da$Sequence)

df <- as.data.frame(phylo_rar_then_species@tax_table)

view(subset(df, sequence %in% da$Sequence))
da <- subset(da, select = -Taxon)

colnames(df) <- paste(colnames(df), "phyloseq", sep = "_")
df$Sequence <- row.names(df)

df <- left_join(da, df, by = "Sequence")

new_col_order <- colnames(df)[1:10]
new_col_order <- c(new_col_order, "kingdom", "Kingdom_phyloseq",
                   "Phylum", "Phylum_phyloseq",
                   "Class", "Class_phyloseq",
                   "Order", "Order_phyloseq",
                   "Family", "Family_phyloseq",
                   "Genus", "Genus_phyloseq",
                   "Species_phyloseq",
                   "Sequence")

df <- df[,new_col_order]
view(df)

## Feature importance -----------------

feat_imp <- read.table("feature_importance.tsv", sep = "\t", header = TRUE)
head(feat_imp)
summary(feat_imp$importance)
dim(feat_imp) # 9388 features
feat_imp <- subset(feat_imp, importance != 0)
dim(feat_imp) # 704 with >0 importance

# Qiime generated a heat map with feature labels but no text files of them
# I used google docs to try to transcribe

heat_feats <- read.csv("heatmap_features.csv", header = FALSE)
# head(heat_feats)
colnames(heat_feats) <- ("id")

sum(heat_feats$id %in% feat_imp$id)
nrow(heat_feats)
# 49 out of 50, good job google

heat_feats <- left_join(heat_feats, feat_imp, by = "id")
head(heat_feats)
# features <- rename(features, "id" = "Feature.ID")
# head(features)
sum(heat_feats$id %in% features$id)
# heat_feats <- left_join(heat_feats, features)
view(heat_feats)


heat_feats <- read.csv("heatmap_features.csv", header = FALSE)
colnames(heat_feats) <- ("id")
heat_feats <- left_join(heat_feats, feat_imp, by = "id")
heat_feats <- left_join(heat_feats, tax, by = "id")
view(heat_feats)

heat_feats <- split_taxon_column(heat_feats)

write.csv(heat_feats, "heat_feats.csv", row.names = FALSE, na = "")
head(heat_feats)


## Random forest -------------------


# Make a dataframe of training data with OTUs as column and samples as rows, which is the phyloseq OTU table
phylo_rar_then_species <- readRDS("phylo_rar_then_species.rds")
predictors <- otu_table(phylo_rar_then_species, taxa_are_rows = FALSE)

dim(predictors)
# 76 samples, 4356 SVs

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(phylo_rar_then_species), stringsAsFactors = FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus = ifelse(is.na(tax.df$Genus), 
                      paste(tax.df$Family), 
                      paste(tax.df$Genus)) 

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)
tax.df$Genus.species <- str_remove_all(tax.df$Genus.species, " NA")

# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

# clean up some of the other taxa info
# I think I want to keep these
# colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
# colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))



### start here when choosing factors, can reuse the above lines as needed. one example for factorial data, and one for numeric data is provided. select as needed.

# Make one column for our outcome/response variable. Choose which one applies to the thing you want to test, and then follow factorial or numeric through the rest of the code section.
response <- as.factor(sample_data(phylo_rar_then_species)$SarsCov2) 

# Combine response and SVs into data frame
rf.data <- data.frame(response, predictors)


# set seed for random number generation reproducibility
set.seed(2)

# classify for factorial data
response.pf <- rfPermute(response ~. , data = rf.data, na.action = na.omit, ntree = 500, nrep = 100) #na.omit ignores NAs in data (not tolerated). ntrees is how many forests to build, nreps generates p-value

print(response.pf)

# paste the print out here, especially the OOB error. 1-(Out-of-the-box error) = accuracy of your model
  
  # neg pos pct.correct LCI_0.95 UCI_0.95
  # neg      36   2        94.7     82.3     99.4
  # pos       3  35        92.1     78.6     98.3
  # Overall  NA  NA        93.4     85.3     97.8


saveRDS(response.pf, "response_pf.rds")
response.pf <- readRDS("response_pf.rds")


# grab which features were labeled "important"
imp <- importance(response.pf, scale = TRUE)

# Make a data frame with predictor names and their importance
imp.df <- data.frame(predictors = rownames(imp), imp) 


# For factorial data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, MeanDecreaseAccuracy.pval <= 0.05) 
print(dim(imp.sig)) # 59 x 9

# or For factorial data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$MeanDecreaseAccuracy),]

#create levels to the factor based on SV table
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top so many predictors (more than 50 is a crowded graph)
# imp.top <- imp.sort[1:50, ] # 50 was too busy
imp.top <- imp.sort[1:30, ]

# figure out what they are and name them, otherwise will just be named the full SV
otunames <- imp.top$predictors
# grab the column names from the otu table that match those in the forest set
pred.abun.colnum <- which(colnames(rf.data) %in% otunames)

# when you find a match, grad the abudnance data
pred.abun <- rf.data[,sort(c(pred.abun.colnum))]

# make this into a dataframe for manipulation
pred.abun.df <- data.frame(pred.abun, stringsAsFactors = FALSE)

# use the row.names (sample names) from the phyloseq object to name the samples in your forest
row.names(pred.abun.df) <- row.names(sample_data(phylo_rar_then_species))

# add some factors that you can use to make your graph pretty, as many as you want
pred.abun.df$Sample <- row.names(sample_data(phylo_rar_then_species)) #always grab the sample names
pred.abun.df$SarsCov2 <- sample_data(phylo_rar_then_species)$SarsCov2

# head(pred.abun.df)
# melt and transform the data using ALL the factors you added
m <- melt(pred.abun.df, id.vars = c("Sample", "SarsCov2"))

head(m)


m$rescale <- log(1 + as.numeric(m$value))

rf_model <- m
rf_model$variable <- str_replace_all(rf_model$variable, "NA", "SV")
saveRDS(rf_model, "rf_model.rds")


# some alignment stuff -----------

head(heat_feats)

mito_seqs <- subset(heat_feats, Family == "Mitochondria")
mito_seqs <- subset(mito_seqs, select = -Taxon)

seqs <- readDNAStringSet("sequences.fasta", format = "fasta")
seqs <- as.data.frame(seqs)
seqs$id <- row.names(seqs)
seqs <- rename(seqs, "Sequence" = x)


mito_seqs <- left_join(mito_seqs, seqs, by = "id")
head(mito_seqs)
head(DNAStringSet(mito_seqs$Sequence))

write.fasta(sequences = DNAStringSet(mito_seqs$Sequence), 
            names = mito_seqs$id, 
            file.out = "mito_seqs.fasta")

# BiocManager::install("msa")
# library(msa)

# mySequences <- readDNAStringSet("mito_seqs.fasta")

mySequences <- DNAStringSet(mito_seqs$Sequence)
names(mySequences) <- mito_seqs$id
head(mySequences)
mySequences


myFirstAlignment <- msa(mySequences)

msaPrettyPrint(myFirstAlignment, 
               output = "pdf", 
               showNames = "left",
               showLogo = "none", 
               askForOverwrite = FALSE, 
               verbose = FALSE)





