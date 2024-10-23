## Intro and package loading ---------------------------
setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/cow_acidosis_rumen_data") 

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab2_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab3_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")
source("/Users/emilykibbler/Desktop/projects/R/AVS_554/load_data.R")
# source("AVS554_packages.R")
# install_necessary()
# install_optional()

library(tidyverse)
library(devtools)
library(Rqc)
library(BiocParallel)
library(Rqc)
library(Biostrings)
library(dada2); packageVersion('dada2')
library(beepr)
library(phyloseq)
library(decontam)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)
library(ggtext)
library(lme4)
library(lmerTest)
library(emmeans)
library(remotes)
library(PerformanceAnalytics)

#this wd folder name will be recorded in the metadata data so try to make it something descriptive
setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554/cow_acidosis_rumen_data") 

# "Types" of data could be forward and reverse, batches, etc.
# Code expects each "type" to have a corresponding folder with format type_raw in the wd

file_info <- metadata_import(types = c("forward", "reverse"))

# One plot for each data type will be in the wd in .png format
# Default data type is raw
qualplots(file_info); beep()


# still needs troubleshooting
for (i in 1:nrow(file_info)) {
    rqc(path = file_info$raw_data_dir[i][1], 
      sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
       file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())
}; beep()


# 5. Optional: search for strings in your data
search_strings("AAAAAAAAAAA", file_info$file_names_full[[1]][[9]]) # pick a random one to search, as an example


# Filter based on the plots
  # For the cow acidosis data, there is nothing usable for reverse

# See the functions R file for other arguments I have left as defaults, such as remove phi x = TRUE


subset(file_info, data_subset == "forward") %>%
  filtering(trunclen = 150); beep() # trimleft default is 10, trimright default is 0, trunclen default is 0

filtoutput <- readRDS("cow_acidosis_rumen_data_forward_filtered_output.rds")
file_info <- readRDS("file_info.rds")
qualplots(file_info, type = "filtered"); beep()
# qualplots(file_info, type = "raw"); beep()



ggplot(as.data.frame(filtoutput)) + 
  geom_point(aes(row.names(filtoutput), reads.in),color = "blue") + geom_point(aes(row.names(filtoutput), reads.out), color = "orange")



full_filt_fns_rev <- list.files(path_filt_rev, full.names = TRUE)

plotQualityProfile(full_filt_fns_rev, aggregate = T); beep() # aggregates all forward reads/read1
ggsave("filt_rev_qplot.png")

rqc(path = path_filt_rev, 
    sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
    outdir = tempdir(), file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())

# You want to end up with

# check the dimensions of the variable you created, outputs two numbers: rows (# samples), columns (# of info it added)
dim(filtoutput) # 34 x 2

# take a look at the counts
head(filtoutput)
filtoutput <- as.data.frame(filtoutput)

# order the info by the first column (reads.in)
View(filtoutput[order(filtoutput[,1], decreasing = FALSE),])

# order the info by the second column (reads.out) which are filtered reads
View(filtoutput[order(filtoutput[,2], decreasing = FALSE),])

# get a sum total for raw reads in and filtered reads out
colSums(filtoutput)

# look at trends

ggplot(filtoutput) + 
  geom_point(aes(row.names(filtoutput), reads.in),color = "blue") + 
  geom_point(aes(row.names(filtoutput), reads.out), color = "orange")


### Summary of filtering ----
# Which parameters did you use and why? 
# I threw out reverse reads due to poor quality. I trimmed at 150, where the quality was starting to dip.
# What was the largest and small number of raw reads?
  # Largest: E53_ITS_for.fastq.gz -- 99336
  # Smallest: F29_ITS_for.fastq.gz -- 10545
# What was the largest and small number of filtered reads?
  # Largest: E53_ITS_for.fastq.gz -- 93899
  # Smallest: F29_ITS_for.fastq.gz -- 9756

# How many total raw reads?
# FIXME
# How many total filtered reads?
# FIXME

