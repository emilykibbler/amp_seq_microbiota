## Start -------
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
# 
# ggplot(phylo_rar_rich_df, aes(x = SarsCov2, y = Observed)) + 
#   theme_minimal() + 
#   geom_boxplot(aes(color = SarsCov2)) +
#   theme(legend.position = "none") +
#   xlab("") + 
#   ylab("Bacterial Richness (SVs)") + 
#   # theme(text = element_text(size = 20)) # increases font size
#   ggtitle("Observed bacterial richness")
