# setwd("/Users/emilykibbler/Desktop/projects/R/AVS_554")

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
