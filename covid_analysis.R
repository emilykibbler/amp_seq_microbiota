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

## Initial quality inspection, filter and trim
file_info <- metadata_import(types = c("fwd", "rev")) # saves a copy as an rds
qualplots(file_info)

filtoutput <- filterAndTrim( 
  file_info$file_names_full[[1]], 
  file.path(file_info$filt_dir[1], paste0(file_info$sample_names[[1]], "_F_filt.fastq.gz")), 
  file_info$file_names_full[[2]], 
  file.path(file_info$filt_dir[2], paste0(file_info$sample_names[[1]], "_R_filt.fastq.gz")), 
  trimLeft = c(10, 10), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina
  truncLen = c(250, 250), 
  maxEE = c(5,5), # max errors tolerated
  verbose = FALSE)
saveRDS(filtoutput, "filtoutput.rds")

qualplots(file_info, type = "filtered")


