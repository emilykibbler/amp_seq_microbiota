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


# Lab 6_A: DADA2 assigning taxonomy from a reference database file (Lab 6) ---------------------------------


# Be sure to have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 

# The Silva file shown below is used for 16S rRNA (prokaryotes) and nicely formatted versions can be downloaded from the DADA2 website, which also contains some (18S) eukaryotic versions.  Let me know if you need a custom one for your dataset.
#
# 1. take the SVs (otu table) from the phyloseq object and make it into a dataframe again


phylo_decontam <- readRDS("phylo_clean_decontam_decontam.rds")
summary(rowSums(phylo_decontam@otu_table)[1:76]) # 1-76 are my experimental samples, after that are lab negative controls
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6112   43586   66853   70767   92416  200156 

seqtab_nochim_decontam = as(otu_table(phylo_decontam), "matrix")
seqtab_nochim_decontam = as.matrix(seqtab_nochim_decontam)

view(readDNAStringSet("trainset16_022016.rdp.tgz", format = "fasta"))
# Error in .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec,  : 
# reading FASTA file /Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/trainset16_022016.rdp.gz: ">" expected at beginning of line 1
view(readDNAStringSet("silva_nr99_v138.1_train_set.fa.gz", format = "fasta")) # shows a list of sequences with no names
# I can't unzip the above file but it just works


view(readDNAStringSet("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/trainset16_022016.rdp/trainset16_022016.rdp.fasta", format = "fasta"))

# works, name of first sequence, and they all look kind of like this:
# AJ000684|S000004347 Root;Bacteria;"Actinobacteria";Actinobacteria;Actinobacteridae;Actinomycetales;Corynebacterineae;Mycobacteriaceae;Mycobacterium
# the folder used above, the unzipped .tgz I downloaded, also has a file called .tax that I don't know how to open


RDP_seqs <- read.fasta("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/trainset16_022016.rdp/trainset16_022016.rdp.fasta")

RDP_tax_names <- read.delim("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/trainset16_022016.rdp/trainset16_022016.rdp.tax", sep = "\t", header = FALSE)

write.fasta(sequences = RDP_seqs, names = RDP_tax_names$V2, file.out = "rdp_train16.fas")

# view(read.fasta("rdp_train16.fas"))


RDP_decontam_all_taxa <- assignTaxonomy(seqtab_nochim_decontam, 
                                    "trainset16_022016.rdp.gz",       
                                    tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                    minBoot = 75, verbose = TRUE) # ; beep("treasure") 
# same error as above
# Error in .Call2("fasta_index", filexp_list, nrec, skip, seek.first.rec,  : 
#                   reading FASTA file trainset16_022016.rdp.gz: ">" expected at beginning of line 1

# This is from the assignTaxonomy help page, the snippet about the second argument:
  # refFasta	
  # (Required). The path to the reference fasta file, or an R connection Can be compressed. This reference fasta file should be formatted so that the id lines correspond to the taxonomy (or classification) of the associated sequence, and each taxonomic level is separated by a semicolon. Eg.
  # 
  # >Kingom;Phylum;Class;Order;Family;Genus; ACGAATGTGAAGTAA......

RDP_decontam_all_taxa <- assignTaxonomy(seqtab_nochim_decontam, 
                                        "/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/trainset16_022016.rdp/trainset16_022016.rdp.fasta",       
                                        tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                        minBoot = 75, verbose = TRUE) # ; beep("treasure") 

# Error: Invalid map from references to genus.
# In addition: Warning message:
#   In matrix(unlist(strsplit(genus.unq, ";")), ncol = td, byrow = TRUE) :
#   data length [108541] is not a sub-multiple or multiple of the number of rows [12061]
# 

# trying one I mashed together
RDP_decontam_all_taxa <- assignTaxonomy(seqtab_nochim_decontam, 
                                        "rdp_train16.fas",       
                                        tryRC = TRUE, # If you didn't get any IDs the first time, Use this to try the reverse complement of your sequences instead
                                        minBoot = 75, verbose = TRUE) # ; beep("treasure") 


saveRDS(RDP_decontam_all_taxa, 'RDP_decontam_all_taxa.rds')
write.csv(RDP_decontam_all_taxa, 'RDP_decontam_all_taxa.csv')

decontam_all_taxa <- readRDS("decontam_all_taxa.rds")

# how many have unassigned SVs?
# Silva:
decontam_all_taxa <- as.data.frame(decontam_all_taxa)
nrow(subset(decontam_all_taxa, is.na(Kingdom))) # 227
nrow(decontam_all_taxa) # 5074
nrow(subset(decontam_all_taxa, is.na(Kingdom)))/nrow(decontam_all_taxa) # 0.045
# RDP:
RDP_decontam_all_taxa <- as.data.frame(RDP_decontam_all_taxa)
nrow(subset(RDP_decontam_all_taxa, is.na(Kingdom))) # 927
nrow(RDP_decontam_all_taxa) # 5074
nrow(subset(RDP_decontam_all_taxa, is.na(Kingdom)))/nrow(RDP_decontam_all_taxa) # 0.183

# OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.





RDP_decontam_all_taxa_species <- addSpecies(RDP_decontam_all_taxa, 
                                        'silva_species_assignment_v138.1.fa.gz', 
                                        allowMultiple = FALSE, 
                                        verbose = FALSE) ; beep("treasure") 

saveRDS(RDP_decontam_all_taxa_species, 'RDP_decontam_all_taxa_species.rds') 
write.csv(RDP_decontam_all_taxa_species, 'RDP_decontam_all_taxa_species.csv') 