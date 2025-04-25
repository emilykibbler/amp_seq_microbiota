## Lab 2 ---------

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


## Lab 3--------

create_and_plot_error_profile <- function(info, bases = 1e6) {
  set.seed(1234) 
  for (i in 1:nrow(info)) {
    # 4. Learn the error rates from your sequencing run. Some error is inherent to sequencing run, so this should be done on each run separately, and you can combine multiple datasets later at the sequence table step. Be sure to include controls or low quality samples for learning error. When using smaller datasets, DADA2 recommended nbases=1e6 reads, but for the big data workflow, it recommended nbases=2e6. Multithread = TRUE if you are using a mac that can handle it.
    filtered_file_names <- list.files(info$filt_dir[i], full.names = TRUE)
    # if there are any other types of files, take them out
    filtered_file_names <- filtered_file_names[grep("_filt.fastq", filtered_file_names)]
    if (length(filtered_file_names) == 0) {
      print("No _filt data files found in", gsub("^.*/", "", info$filt_dir[i]))
      if (i != nrow(info)) {
        print("Proceeding, looking for other filtered data subsets")
      }
    }
    err <- learnErrors(filtered_file_names, nbases = bases, multithread = TRUE, randomize = TRUE) #I can use multithread=TRUE because I'm on a mac
    # 5. Save the error profiles to output files
    saveRDS(err, paste0("error_profile_", info$data_subset[i], ".rds")) 
    plot_name <- paste(info$project[i], info$data_subset[i], "filtered_error_profile", sep = "_")
    p <- plotErrors(err, nominalQ = TRUE) +
      ggtitle(str_to_title(str_replace_all(plot_name, "_", " ")))
    ggsave(paste0(plot_name, ".png"), p)
    print(p)
  }
}



make_SV_summary <- function(seqtab_input) {
  if (class(seqtab_input) != "data.frame") {
    # print("Data frame is the expected input")
    # print(paste("Your input is:", class(seqtab_input)))
    stop(paste("Data frame is the expected input.", "\n", "Your input is:", class(seqtab_input)))
  }
  summary <- data.frame(matrix(ncol = 2, nrow = 0)) # empty df to dump row-wise data into
  colnames(summary) <- c("sample", "SVs")
  for (i in 1:nrow(seqtab_input)) {
    temp <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(temp) <- c("sample", "SVs")
    # For each row in seqtab, report how many SVs (columns of seqtab) with reads there were
    temp[1,] <- c(rownames(seqtab_input)[i], length(seqtab_input[i,][seqtab_input[i,] > 0]))
    summary <- rbind(temp, summary)
  }
  read_totals <- as.data.frame(rowSums(seqtab_input))
  read_totals$sample <- row.names(read_totals)
  colnames(read_totals) <- c("filt_reads", "sample")
  summary <- merge(summary, read_totals, by = "sample")
  summary$SVs <- as.numeric(summary$SVs)
  return(summary)
}
