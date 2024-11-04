# 5. Optional: search for strings in your data

search_strings <- function(bases, fastq_dot_gz) {
  ## specify one file name to look through
  DNA_to_read <- readDNAStringSet(fastq_dot_gz, format = "fastq", with.qualities = TRUE) # FIXME
  # count the number of different nucleotides
  alphabetFrequency(DNA_to_read)
  # look for ambiguous bases
  vmatchPattern("N", DNA_to_read)
  # look for a homopolymer. You can change the AAAs to other letters or number of letters. you can also use this to search for indexes, primers, adaptors, or other strings
  vmatchPattern(bases, DNA_to_read, fixed = FALSE) 
}

# Importing data, preliminary formatting
# Default expectation is fwd and rev reads but can be changed to batch 1 and batch 2, for example
# "Types" must match folder names. If types are "batch1, batch2" the raw data must be found in folders of those names followed by "_raw"
# "Parent directory" needs to be defined previously (unless everything is in the wd, then p_dir is defined as "")

metadata_import <- function(types = c("fwd", "rev")) {
  info <- as.data.frame(matrix(nrow = 0, ncol = 6))
  colnames(info) <- c("data_subset", "raw_data_dir", "file_names", "file_names_full", "sample_names", "filt_dir")
  for (i in 1:length(types)) {
    temp <- as.data.frame(matrix(nrow = 1, ncol = 6))
    colnames(temp) <- c("data_subset", "raw_data_dir", "file_names", "file_names_full", "sample_names", "filt_dir")
    temp$data_subset[1] <- types[i] #name of this sub-data set
    temp$raw_data_dir[1] <- paste0(getwd(), "/", types[i], "_raw") # folder for data files
    if (!(temp$raw_data_dir[1] %in% list.dirs(getwd(), recursive = FALSE))) {
      stop("Folder names not as expected, check that types have [type]_raw folders in the wd")
    }
    temp$file_names[1] <- list(list.files(paste0(getwd(), "/", types[i], "_raw"))) # files, names only
    temp$file_names_full[1] <- list(list.files(paste0(getwd(),"/", types[i], "_raw"), full.names = TRUE))    # files with full paths
    temp$sample_names[1] <- list(str_remove_all(temp$file_names[[1]], ".fastq.gz"))
    # If this code has already been run, the _filtered folder will already exist
    if (!(paste0(getwd(), "/", types[i], "_filtered") %in% list.dirs(paste0(getwd()), recursive = FALSE))) {
      # If it doesn't already exist, create now
      dir.create(paste0(getwd(), "/", types[i], "_filtered"))
    }
    temp$filt_dir <- paste0(getwd(), "/", types[i], "_filtered")
    info <- rbind(info, temp)
  }
  info$project <- gsub("^.*/", "", getwd()) # name the project by whatever the wd is
  saveRDS(info, file = "file_info.rds")
  return(info)
}

# This function works but it definitely needs to be simplified

qualplots <- function(info, type = "raw") {
  if (type == "raw" | type == "Raw" | type == "RAW") {
    for (i in 1:nrow(info)) {
      plot_name <- paste(info[i,]$project, info[i,]$data_subset, "quality plot")
      temp <- info[i,]$file_names_full[[1]] # the files we want to analyze for this step
      p <- plotQualityProfile(temp, aggregate = TRUE) +
          ggtitle(str_to_title(str_replace_all(plot_name, "_", " ")))
      print(p)
      ggsave(plot = p, filename = paste0(str_replace_all(plot_name, " ", "_"), ".png"))
    }
  }
  else if (type == "filtered" | type == "Filtered" | type == "FILTERED") {
    if (sum(grep("filt.fastq.gz", list.files(getwd(), recursive = TRUE))) == 0) {
      stop("No filtered files found")
    }
    for (i in 1:nrow(info)) {
      if (sum(grep("filt.fastq.gz", list.files(paste0(getwd(), "/", info$data_subset[i], "_filtered")))) != 0) {
        plot_name <- paste(info[i,]$project, info[i,]$data_subset, "filtered quality plot")
        temp <- list.files(paste0(getwd(), "/", info$data_subset[i], "_filtered"), full.names = TRUE)
        temp <- temp[grep("fastq.gz", temp)] # take out non-fastq files
        p <- plotQualityProfile(temp, aggregate = TRUE) +
          ggtitle(str_to_title(str_replace_all(plot_name, "_", " ")))
        print(p)
        ggsave(plot = p, filename = paste0(str_replace_all(plot_name, " ", "_"), ".png"))
      }
    }
  } else {print("Type of data should be raw or filtered")}
}

filtering <- function(dat, trimleft = 10, trimright = 0, trunclen = 0) {
  if (sum(grep("filt.fastq", list.files(paste0(getwd(), "/", dat$data_subset[1], "_filtered")))) != 0) {
    stop("There are already filtered files. Delete or move them.")
  }
  filtoutput <- filterAndTrim(
    file.path(dat$raw_data_dir[1], dat$file_names[[1]]), #, fastqs_fwd), # path to and file names of the raw reads
    file.path(paste0(getwd(), "/", dat$data_subset[1], "_filtered"), paste0(dat$sample_names[[1]], "_filt.fastq.gz")), # set filtered file paths and filtered file names to create
    trimLeft = trimleft, # cuts off the first XX bases from the F reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
    trimRight = trimright, # cuts of last XX bases, or hash this out and use truncLen
    truncLen = trunclen, # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
    #minLen=10, # optional: remove any sequences longer than max length, for use with 454 if don't want to truncate
    maxEE = 2, # the maximum number of expected errors allowed in a read, always >1
    maxN = 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any!!
    rm.phix = TRUE, # remove any PhiX DNA (used as positive control),
    verbose = TRUE, matchIDs = TRUE) # verbose = print to screen
  saveRDS(filtoutput, paste0(dat$project, "_", dat$data_subset, "_filtered_output.rds")) 
  
  plot_name <- paste(dat$project, dat$data_subset, "reads in and out")
  df <- as.data.frame(filtoutput)
  df$sample_name <- row.names(df)
  p <- ggplot(df) +
        geom_point(aes(sample_name, reads.in, color = "Reads in")) + 
        geom_point(aes(sample_name, reads.out, color = "Reads out")) +
        theme(axis.text.x = element_text(angle = 45, size = 4)) +
        theme(legend.title = element_blank()) +
        ylab("Reads") +
        xlab("Sample") +
        ggtitle(str_to_title(str_replace_all(plot_name, "_", " ")))
  print(p)
  ggsave(plot = p, filename = paste0(str_replace_all(plot_name, " ", "_"), ".png"))
}



# 
# ggplot(as.data.frame(filtoutput)) + 
#   geom_point(aes(row.names(filtoutput), reads.in),color = "blue") + geom_point(aes(row.names(filtoutput), reads.out), color = "orange")



# for (i in 1:nrow(df)) {
#   temp <-  df[i,]$y_vals[[1]]
#   p <- ggplot() +
#     geom_boxplot(aes(y = temp))
#   ggsave(filename = paste0(i, "_randplot.png"), plot = p)
#   print(p)
# }


 

# p <- plotQualityProfile(file_info[1,]$file_names_full[[1]], aggregate = TRUE) +
#   ggtitle(paste(file_info[1,]$project, file_info[1,]$data_subset, "quality plot"))
# ggsave(plot = p, filename = paste(file_info[1,]$project, file_info[1,]$data_subset, "quality_plot.png", sep = "_")) 
  
#   
#   plotQualityProfile(full_fns_raw_fwd, aggregate = T); beep() # aggregates all forward reads/read1
# ggsave("raw_lamb_f_qplot.png")
# 
# plotQualityProfile(full_fns_raw_rev, aggregate = T); beep() # aggregates all reverse reads/read2
# ggsave(paste0(parent_dir, "/raw_lamb_r_qplot.png"))
# 
# 
# # dir.create(paste0(parent_dir, "/filtered_b1")) 
# 
# 
# # parent_dir <- "/Users/emilykibbler/Desktop/projects/R/AVS_554/scallop_data"
# 
# # "Batch 1" will be abbreviated as b1 and "Batch 2" as b2
# # Sequencing data is all forward reads but was done twice in technical repeats
# # Rev data was so poor it was not made available by Dr. Isqhaq
# 
# path_raw_b1 <- paste0(parent_dir, "/batch1_raw") # folder for data files
# path_raw_b2 <- paste0(parent_dir, "/batch2_raw")
# fns_raw_b1 <- list.files(path_raw_b1) # file names
# fns_raw_b2 <- list.files(path_raw_b2) # Fully specify the path, and file names (fns) for the reads.This will come in handy later
# 
# full_fns_raw_b1 <- list.files(path_raw_b1, full.names = TRUE) # file names with the path
# full_fns_raw_b2 <- list.files(path_raw_b2, full.names = TRUE)
# 
# 
# # 2. Specify what file names (fastq files) to look for in that folder based on the type of file extension. As written, this looks for only zipped (gz) files.
# 
# fastqs_b1 <- fns_raw_b1[grepl('.gz', fns_raw_b1)]
# fastqs_b2 <- fns_raw_b2[grepl('.gz', fns_raw_b2)] 
# # head(fastqs_b1)
# 
# # 4. Pull the sample names by reading the forward read file names, cutting it, and taking the sample name chunk.
# sample_names_b1 <- sapply(strsplit(basename(fastqs_b1), '.fas'), `[`, 1) 
# sample_names_b2 <- sapply(strsplit(basename(fastqs_b2), '.fas'), `[`, 1) 
# 
# 
# # 5. Assign the sample names to the list of Read 1 forward fastq files
# names(fastqs_b1) <- sample_names_b1
# names(fastqs_b2) <- sample_names_b2