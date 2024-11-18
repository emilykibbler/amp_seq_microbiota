# For some reason the subset_samples function, which is part of the phyloseq packages, does not work in a loop or function
# This is a workaround I found on github
# Note to self, rename function to something more profesh later

subset_me_bro <- function(starting_object, category, subset_on) {
  search_for <- as.character(get_variable(starting_object, category)) %in% subset_on
  the_subset <- prune_samples(search_for, starting_object)
  return(the_subset)
}

subset_and_trim <- function(starting_object, category, subset_on) {
  search_for <- as.character(get_variable(starting_object, category)) %in% subset_on
  the_subset <- prune_samples(search_for, starting_object)
  the_subset <- prune_taxa(taxa_sums(the_subset) > 0, the_subset)
  return(the_subset)
}


# This function takes a phyloseq object and returns clean data
# At least type of negative control is required
# neg_con is for experimental controls to be applied to all data
# batch_con is for extraction batch controls
# Controls must be in Sample_type column
# Note to self, update later to accept any column


clean_phylo_data <- function(phyloseq_object, neg_con = list(), batch_con = list()) {
  control_names <- c(neg_con, batch_con)
  summary <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(summary) <- c("Description", "taxa_found", "samples_used", "control_name")
  summary[1,] <- c("Start", ncol(phyloseq_object@otu_table), nrow(phyloseq_object@otu_table), NA) # initiate a summary table to tell user what happened at each step
  # Error checks
  if (length(control_names) == 0) {
    stop("No negative controls provided to the function")
  }
  if (length(as.character(get_variable(phyloseq_object, "Sample_type")) %in% control_names) == 0) {
    stop("Negative control names not found in the data")
  }
  # Experimental clean
  experimental_clean <- phyloseq_object
  if (length(neg_con) > 0) {
    temp <- subset_me_bro(phyloseq_object, "Sample_type", neg_con)
    temp <- prune_taxa(taxa_sums(temp) > 0, temp) # Take out taxa aka SVs with no reads
    control_vec <- as.vector(taxa_names(temp)) #Make the taxa names into a vector so you can remove them. 
    vec <- as.vector(taxa_names(phyloseq_object))
    keep <- setdiff(vec, control_vec)
    experimental_clean <- prune_taxa(keep, experimental_clean) # Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
    experimental_clean <- prune_samples(sample_sums(experimental_clean) > 0, experimental_clean) # Then remove the samples which are now empty (neg cons)
    summary[nrow(summary) + 1, ] <- c("Experimental negative", ncol(temp@otu_table), nrow(temp@otu_table), neg_con)
    summary[nrow(summary) + 1, ] <- c("Data after experimental control clean", ncol(experimental_clean@otu_table), nrow(experimental_clean@otu_table), NA)
  } # end of experimental clean
  # At this point: experimental_cleaned either has the first clean done, if relevant, and if not, is just the phyloseq input object
  # Batch cleanup
  clean <- experimental_clean
  if (length(batch_con) > 0) {
    batches_list <- unique(as.character(get_variable(experimental_clean, "DNA_extraction_batch")))
    if (length(batches_list) == 0) {
      print("Column called DNA_extraction_batch must be in metadata; go back and rename if necessary")
    } else {
      for (i in 1:length(batches_list)) {
        batch <- subset_me_bro(experimental_clean, "DNA_extraction_batch", batches_list[i])
        summary[nrow(summary) + 1, ] <- c(paste("Batch", batches_list[i], "samples"), ncol(batch@otu_table), nrow(batch@otu_table), NA)
        if (batch_con %in% get_variable(batch, "Sample_type")) { # if the control is not in the batch, skip 
          batch_controls <- subset_me_bro(batch, "Sample_type", batch_con)
          batch_controls <- prune_taxa(taxa_sums(batch_controls) > 0, batch_controls) # Take out taxa aka SVs with no reads
          summary[nrow(summary) + 1, ] <- c(paste("Batch", batches_list[i], "control"), ncol(batch_controls@otu_table), nrow(batch_controls@otu_table), batch_con)
          control_vec <- as.vector(taxa_names(batch_controls)) # make the taxa names into a vector so you can remove them 
          vec <- as.vector(taxa_names(batch))
          keep <- setdiff(vec, control_vec) ## then, use the keep vector for the prune taxa argument
          batch_clean <- prune_taxa(keep, batch)
          
          if (identical(clean, experimental_clean)) { # this would only be true the first time this loop executes so it's OK to overwrite
            clean <- batch_clean
          } else {
            clean <- merge_phyloseq(clean, batch_clean) # every subsequent loop adds in
          }
        }
      } # end of batches for-loop
      clean <- prune_taxa(taxa_sums(clean) > 0, clean)
      clean <- prune_samples(sample_sums(clean) > 0, clean)
    } 
  } # end of if statement that executes if batches are included
  # Put it all together
  summary[nrow(summary) + 1, ] <- c("Total removed", ncol(phyloseq_object@otu_table) - ncol(clean@otu_table), nrow(phyloseq_object@otu_table) - nrow(clean@otu_table), NA)
  summary[nrow(summary) + 1, ] <- c("Final data after clean", ncol(clean@otu_table), nrow(clean@otu_table), NA)
  print(summary)
  saveRDS(summary, "cleaning_summary.rds")
  return(clean)
}

# Making parameters a little less strict  
# Threshold is where you want to cut off how many reads in the negative is cleared out
# That is, if it shows up in the negative with more than the threshold number of reads,
# the taxa is removed from the data
# Default value is 5
relaxed_clean_phylo_data <- function(phyloseq_object, neg_con = list(), batch_con = list(), threshold = 5) {
  control_names <- c(neg_con, batch_con)
  summary <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(summary) <- c("Description", "taxa_found", "samples_used", "control_name")
  summary[1,] <- c("Start", ncol(phyloseq_object@otu_table), nrow(phyloseq_object@otu_table), NA) # initiate a summary table to tell user what happened at each step
  # Error checks
  if (length(control_names) == 0) {
    stop("No negative controls provided to the function")
  }
  if (length(as.character(get_variable(phyloseq_object, "Sample_type")) %in% control_names) == 0) {
    stop("Negative control names not found in the data")
  }
  # Experimental clean
  experimental_clean <- phyloseq_object
  if (length(neg_con) > 0) {
    temp <- subset_me_bro(phyloseq_object, "Sample_type", neg_con)
    temp <- prune_taxa(taxa_sums(temp) > threshold, temp) # Take out taxa aka SVs with almost no reads
    control_vec <- as.vector(taxa_names(temp)) #Make the taxa names into a vector so you can remove them. 
    vec <- as.vector(taxa_names(phyloseq_object))
    keep <- setdiff(vec, control_vec)
    experimental_clean <- prune_taxa(keep, experimental_clean) # Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
    experimental_clean <- prune_samples(sample_sums(experimental_clean) > 0, experimental_clean) # Then remove the samples which are now empty (neg cons)
    summary[nrow(summary) + 1, ] <- c("Experimental negative", ncol(temp@otu_table), nrow(temp@otu_table), neg_con)
    summary[nrow(summary) + 1, ] <- c("Data after experimental control clean", ncol(experimental_clean@otu_table), nrow(experimental_clean@otu_table), NA)
  } # end of experimental clean
  # At this point: experimental_cleaned either has the first clean done, if relevant, and if not, is just the phyloseq input object
  # Batch cleanup
  clean <- experimental_clean
  if (length(batch_con) > 0) {
    batches_list <- unique(as.character(get_variable(experimental_clean, "DNA_extraction_batch")))
    if (length(batches_list) == 0) {
      print("Column called DNA_extraction_batch must be in metadata; go back and rename if necessary")
    } else {
      for (i in 1:length(batches_list)) {
        batch <- subset_me_bro(experimental_clean, "DNA_extraction_batch", batches_list[i])
        summary[nrow(summary) + 1, ] <- c(paste("Batch", batches_list[i], "samples"), ncol(batch@otu_table), nrow(batch@otu_table), NA)
        if (batch_con %in% get_variable(batch, "Sample_type")) { # if the control is not in the batch, skip 
          batch_controls <- subset_me_bro(batch, "Sample_type", batch_con)
          batch_controls <- prune_taxa(taxa_sums(batch_controls) > threshold, batch_controls) # Take out taxa aka SVs with no reads
          summary[nrow(summary) + 1, ] <- c(paste("Batch", batches_list[i], "control"), ncol(batch_controls@otu_table), nrow(batch_controls@otu_table), batch_con)
          control_vec <- as.vector(taxa_names(batch_controls)) # make the taxa names into a vector so you can remove them 
          vec <- as.vector(taxa_names(batch))
          keep <- setdiff(vec, control_vec) ## then, use the keep vector for the prune taxa argument
          batch_clean <- prune_taxa(keep, batch)
          
          if (identical(clean, experimental_clean)) { # this would only be true the first time this loop executes so it's OK to overwrite
            clean <- batch_clean
          } else {
            clean <- merge_phyloseq(clean, batch_clean) # every subsequent loop adds in
          }
        }
      } # end of batches for-loop
      clean <- prune_taxa(taxa_sums(clean) > 0, clean)
      clean <- prune_samples(sample_sums(clean) > 0, clean)
    } 
  } # end of if statement that executes if batches are included
  # Put it all together
  summary[nrow(summary) + 1, ] <- c("Total removed", ncol(phyloseq_object@otu_table) - ncol(clean@otu_table), nrow(phyloseq_object@otu_table) - nrow(clean@otu_table), NA)
  summary[nrow(summary) + 1, ] <- c("Final data after clean", ncol(clean@otu_table), nrow(clean@otu_table), NA)
  print(summary)
  saveRDS(summary, "cleaning_summary.rds")
  return(clean)
}


# Version 1 of this code, which combined all the controls
# Prints where the taxa comes from but does not consider them separately
# I don't know why that makes a difference but it does
# 
# clean_phylo_data <- function(phyloseq_object, neg_con = list(), batch_con = list()) {
#   control_names <- c(neg_con, batch_con)
#   summary <- data.frame(matrix(nrow = 0, ncol = 4))
#   colnames(summary) <- c("Type", "taxa_found", "samples_used", "control_name")
#   summary[1,] <- c("Start", ncol(phyloseq_object@otu_table), nrow(phyloseq_object@otu_table), NA)
#   # Error checks
#   if (length(control_names) == 0) {
#     stop("No negative controls provided to the function")
#   }
#   if (length(as.character(get_variable(phyloseq_object, "Sample_type")) %in% control_names) == 0) {
#     stop("Negative control names not found in the data")
#   }
#   if (length(neg_con) > 0) {
#     search_for <- as.character(get_variable(phyloseq_object, "Sample_type")) %in% neg_con
#     temp <- prune_samples(search_for, phyloseq_object)
#     temp <- prune_taxa(taxa_sums(temp) > 0, temp) # Take out taxa aka SVs with no reads
#     summary[nrow(summary) + 1, ] <- c("Experimental negative", ncol(temp@otu_table), nrow(temp@otu_table), neg_con)
#     summary[nrow(summary) + 1, ] <- c("Data after experimental control clean", ncol(phyloseq_object@otu_table) - ncol(temp@otu_table), nrow(phyloseq_object@otu_table) - nrow(temp@otu_table), NA)
#   }
#   if (length(batch_con) > 0) {
#     search_for <- as.character(get_variable(phyloseq_object, "Sample_type")) %in% batch_con
#     # print(search_for)
#     temp <- prune_samples(search_for, phyloseq_object)
#     batches <- unique(as.character(get_variable(temp, "DNA_extraction_batch")))
#     # print(batches)
#         if (length(batches) == 0) {
#           print("Column called DNA_extraction_batch must be in metadata, go back and rename if necessary")
#         } else {
#             for (i in 1:length(batches)) {
#               search_for <- as.character(get_variable(temp, "DNA_extraction_batch")) == batches[i]
#               temp2 <- prune_samples(search_for, temp)
#               temp2 <- prune_taxa(taxa_sums(temp2) > 0, temp2) # Take out taxa aka SVs with no reads
#               # print(temp2)
#               summary[nrow(summary) + 1, ] <- c(paste("Batch negative", batches[i]), ncol(temp2@otu_table), nrow(temp2@otu_table), batch_con)
#               
#             }
#         }
#   }
# 
#   # controls = subset_samples(phylo, Sample_type == "NegCon_PCR") # can't do this in the function, do the prune workaround
#   controls <- phyloseq_object
#   search_for <- as.character(get_variable(phyloseq_object, "Sample_type")) %in% control_names
#   controls <- prune_samples(search_for, controls)
#   controls <- prune_taxa(taxa_sums(controls) > 0, controls) # Take out taxa aka SVs with no reads
#   # Make the taxa names into a vector so you can remove them. 
#   control_vec <- as.vector(taxa_names(controls))
#   vec <- as.vector(taxa_names(phyloseq_object))
#   keep <- setdiff(vec, control_vec)
#   # Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
#   clean <- prune_taxa(keep, phyloseq_object) # the "keep"s are taxa that are in the data set but not in the negative controls
#   # Then remove the samples which are now empty, namely the NegCon_ samples
#   clean <- prune_samples(sample_sums(clean) > 0, clean)
#   summary[nrow(summary) + 1, ] <- c("Total removed", ncol(phyloseq_object@otu_table) - ncol(clean@otu_table), nrow(phyloseq_object@otu_table) - nrow(clean@otu_table), NA)
#   summary[nrow(summary) + 1, ] <- c("Final data after clean", ncol(clean@otu_table), nrow(clean@otu_table), NA)
#   
#   print(summary)
#   return(clean)
# }
# 
# 
# 



