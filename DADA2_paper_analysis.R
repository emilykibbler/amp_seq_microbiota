setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")


library(tidyverse)
library(devtools)
library(Rqc)
library(BiocParallel)
library(dada2); packageVersion('dada2')
library(beepr)

## Preliminary data processing---------------------------

# prelim data processing, splitting by read 1 and read 2
# queso <- list.files(getwd())
# read1fns <- queso[grep("_1.fastq", queso)]
# read2fns <- queso[grep("_2.fastq", queso)]
# 
# dir.create("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_raw")
# dir.create("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_raw")
# 
length(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_raw")))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_raw"))))
length(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_raw")))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_raw"))))

# missing files for some reason, investigating
all_files <- read.csv("filenames.csv", header = FALSE)
downloaded_files <- list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_raw")
downloaded_files <- c(downloaded_files, list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_raw"))
all_files <- all_files[,1]

downloaded_files <- str_remove_all(downloaded_files, ".fastq.gz")


#  missing_read1 <- list()


salsa <- sapply(downloaded_files, FUN = grep, all_files) # the indexes of all files that we have

queso <- all_files[-salsa]
# This generated a new wget list to download the ones I missed somehow
write.csv(queso, "missing_files.csv", row.names = FALSE)

## Lab 2-------------------------------

# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper")
# source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
# source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")

# "Types" of data could be forward and reverse, batches, etc.
# Code expects each "type" to have a corresponding folder with format type_raw in the wd

file_info <- metadata_import(types = c("read-1", "read-2"))
# saveRDS(file_info, "file_info.RDS")
file_info <- readRDS("file_info.rds")

# One plot for each data type will be in the wd in .png format
# Default data type is raw; can change to "filtered" and run again after filtering
qualplots(file_info); beep()

# run this command to create a quality control report

for (i in 1:nrow(file_info)) {
  rqc(path = file_info$raw_data_dir[i][1], 
      sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
      file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())
} # ; beep()


    # # 5. Optional: search for strings in your data
    # search_strings("AAAAAAAAAAA", file_info$file_names_full[[1]][[9]]) # pick a random one to search, as an example

# Filter separately -- do not use
    # subset(file_info, data_subset == "read-1") %>%
    #   filtering(trunclen = 150)#; beep() # trimleft default is 10, trimright default is 0, trunclen default is 0
    # 
    # subset(file_info, data_subset == "read-2") %>%
    #   filtering(trunclen = 150) #; beep()

# Filter together: this needs to happen for merging seqtab later

filtoutput <- filterAndTrim( # filter reads, takes about 25 minutes
  fwd = file_info$file_names_full[[1]], 
  filt = file_info$filt_dir[1], 
  rev = file_info$file_names_full[[2]], 
  filt.rev = file_info$filt_dir[2],
  trimLeft = c(10, 10), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina
  # trimRight = c(50,50), # cuts of last XX bases from F,R reads, or use truncLen
  truncLen = c(150, 150), # cuts off end of F reads by trimming all to same length, use instead of trimRight
  maxEE = c(2,3), ## The maxEE parameter sets the maximum number of expected errors allowed in a read. Always >1
  verbose = TRUE); beep() # verbose = print to screen
saveRDS(filtoutput, "filtoutput.rds")

# FIXME qualplots isn't working because filtering file names didn't come out correctly.
# Manual fix:
# old_names <- as.character(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_filtered", full.names = TRUE))
# new_names <- str_replace_all(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-1_filtered", full.names = TRUE), ".fastq.gz", "_filt.fastq.gz")
# file.rename(from = old_names, to = new_names)
# old_names <- as.character(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_filtered", full.names = TRUE))
# new_names <- str_replace_all(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper/read-2_filtered", full.names = TRUE), ".fastq.gz", "_filt.fastq.gz")
# file.rename(from = old_names, to = new_names)


qualplots(file_info, type = "filtered") #; beep()
# Reads out: 22086727


ggplot(filtoutput) +
  geom_point(aes(row.names(filtoutput), reads.in, color = "Raw reads")) +
  geom_point(aes(row.names(filtoutput), reads.out, color = "Filtered reads")) +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1)) +
  theme(legend.title = element_blank()) +
  ylab("Reads") +
  xlab("Sample") +
  ggtitle("DADA2_paper_data, raw reads and filtered reads")
ggsave("reads_in_reads_out_DADA2_combined_filtering.png")



## Lab 3  ---------------------------

# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/DADA2_paper")
# source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
# source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")

# Re load if necessary
file_info <- readRDS("file_info.rds")

create_and_plot_error_profile(file_info, bases = 1e6) #; beep()
# 7753760 total bases in 55384 reads from 2 samples will be used for learning the error rates.
# 33610500 total bases in 240075 reads from 1 samples will be used for learning the error rates.

errF <- readRDS("error_profile_read-1.rds")
errR <- readRDS("error_profile_read-2.rds")


# Making the dadas will take about 1.5 hours to run
# Sys.time()
dadaFs <- dada(list.files(subset(file_info, data_subset == "read-1")$filt_dir[1], full.names = TRUE),
               err = errF, 
               multithread = TRUE) #; beep()
saveRDS(dadaFs, "dadaFs.rds")

dadaRs <- dada(list.files(subset(file_info, data_subset == "read-2")$filt_dir[1], full.names = TRUE),
               err = errR, 
               multithread = TRUE) ; beep()
saveRDS(dadaRs, "dadaRs.rds")
# Sys.time()

# This takes about five minutes
mergers <- mergePairs(dadaF = dadaFs,
                      derepF = list.files(file_info$filt_dir[1], full.names = TRUE),
                      dadaR = dadaRs,
                      derepR = list.files(file_info$filt_dir[2], full.names = TRUE),
                      verbose = FALSE) ; beep()
seqtab <- makeSequenceTable(mergers)
# Sys.time()
saveRDS(seqtab, "seqtab.rds")
dim(seqtab) # 64 74948

SVs_found_by_sample <- make_SV_summary(seqtab)
View(SVs_found_by_sample)
# Most SVs: ERR777710_1_filt.fastq.gz -- 15107
# Most reads: ERR777732_1_filt.fastq.gz -- 1541019
# Fewest SVs: there are a few with 0
# Fewest reads: there are a few with 0

# DADA2 Remove chimeras from seqtab (Lab 4) ----------------------------	

# Chimeras are accidentally created sequences during lab protocols. Remove them.	

# 1. Remove chimeras. Leave the method as consensus. multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE) # takes maybe ten minutes
# Identified 14619 bimeras out of 74948 input sequences.
saveRDS(seqtab.nochim, 'seqtab.nochim.rds')

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab.nochim) 
# 64 samples and 60329 SVs (down from 74948)

# Calculate the percentage of chimeras identified out of the total
round(sum(seqtab.nochim)/sum(seqtab), digits = 3) # 0.968 retained
round(1 - sum(seqtab.nochim)/sum(seqtab), digits = 3) # meaning 0.032 removed


