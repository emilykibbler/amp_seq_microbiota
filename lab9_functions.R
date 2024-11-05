# Straight from the example script, just separating it out for simplicity

### test significance of correlations, copy and paste this whole chunk together
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}


core_taxa_finder <- function(phylo_input, params) {
  my_file_name <- paste(deparse(substitute(phylo_input)), "core_taxa.csv", sep = "_")
  if (length(params) != 2) {
    stop("Supply two parameters for the core_members function, detection and prevalence")
  }
  df <- data.frame(tax_table(phylo_input), stringsAsFactors = FALSE) # grab the taxonomy to clean up first
  ## if the Genus is empty, replace with the Family
  df$Genus[is.na(df$Genus)] <- df$Family[is.na(df$Genus)]
  # If the Genus field is over 40 characters, it's a combination of so many things, it's not meaningfully defined
  df$Genus[which(nchar(df$Genus) > 40)] <- "Undefined"
  tax_table(phylo_input) <- as.matrix(df)
  #calculate compositional version of the data (relative abundances)
  core_biomeW_rel <- microbiome::transform(phylo_input, "compositional") 
    # #This returns the taxa that exceed the given prevalence and minimum abundance detection thresholds.
    # Defined but never used? Take out for now. Add to return if I figure out the need
    # core_taxa_standardW <- core_members(core_biomeW.rel, detection = params[1], prevalence = params[2])
  #A full phyloseq object of the core microbiota at those limits is obtained as follows
  phylo_coreW <- core(core_biomeW_rel, detection = params[1], prevalence = params[2])
  
  ###retrieving the associated taxa names from the phyloseq object and add it to what you just made.
  core_taxaW <- taxa(phylo_coreW)
  
  # get the taxonomy data and assign it to what you just made.
  tax_matW <- tax_table(phylo_coreW)
  tax_dfW <- as.data.frame(tax_matW)
  # add the SVs to last column of what you just made.
  tax_dfW$SV <- rownames(tax_dfW)
  # select taxonomy of only those OTUs that are core members based on the thresholds that were used.
  core_taxa_classW <- dplyr::filter(tax_dfW, rownames(tax_dfW) %in% core_taxaW)
  knitr::kable(head(core_taxa_classW))

  # save the list so you can access later, can report just the list
  write.csv(core_taxa_classW, file = my_file_name)

  return(phylo_coreW)
}

# paste(deparse(substitute(phylo)), "core_taxa.csv", sep = "_")


create_top_N_corr_df <- function(phylo_object, number, rich_df) {
  topOTUs <- names(sort(taxa_sums(phylo_object), TRUE)[1:number])
  top <- prune_taxa(topOTUs, phylo_object)
  
  top_sd = as(otu_table(top, taxa_are_rows = FALSE), "matrix")
  top_sd = as.data.frame(top_sd)
  
  # add Genus names in place of the full sequence name that is in the SV columns
  top_tax <- as.data.frame(tax_table(top))
  ## if the Genus is empty, replace with the Family
  top_tax$Genus[is.na(top_tax$Genus)] <- top_tax$Family[is.na(top_tax$Genus)]
  colnames(top_sd) = top_tax$Genus
  
  top_corr_df <- cbind(top_sd, rich_df)
  return(top_corr_df)
}
# 
#   # corrplot requires one dataframe of metadata and/or seqtab data.  Any columns in the dataframe will be correlated against all others.  Too many columns makes the graph unreadable, try to keep it to <50.
#   
#   # select to top 15 most abundant SVs from your phyloseq object and extract to a dataframe
#   #take out top N taxa based on abundance
#   top100OTUs = names(sort(taxa_sums(phylo_decontam_rar), TRUE)[1:100])
#   top100 = prune_taxa(top100OTUs, phylo_decontam_rar)
#   
#   # combine with your metadata and create one dataframe object. you can include other info that you created for a previous dataframe, as long as those objects are still in your R environment. Reminder, you can only use numeric data in a correlation matrix, so you will have to drop certain columns or make them numeric instead.
#   # Coerce to data.frame and add the metadata for these samples
#   top100_sd = as(otu_table(top100, taxa_are_rows = FALSE), "matrix")
#   top100_sd = as.data.frame(top100_sd)
#   
#   # add Genus names in place of the full sequence name that is in the SV columns
#   top100_tax <- as.data.frame(tax_table(top100))
#   ## if the Genus is empty, replace with the Family
#   top100_tax$Genus[is.na(top100_tax$Genus)] <- top100_tax$Family[is.na(top100_tax$Genus)]
#   colnames(top100_sd) = top100_tax$Genus
# 
# # paste all the components together
# # top100_decontam_rar_corr_df <- cbind(top100_sd, phylo_decontam_rar_rich, phylo_decontam_rar_rich, phylo_decontam_rar_sd)
# 
# top100_decontam_rar_corr_df <- cbind(top100_sd, phylo_decontam_rar_rich_df)
# # check header to make sure it looks good
# head(top100_decontam_rar_corr_df)

# # change any column factor names to make them prettier
# names(EX_ps_clean.rar.corr.df)[names(EX_ps_clean.rar.corr.df) == "EX_ps_clean.rar.even"] <- "Evenness"
# 
# # check header to make sure it looks good
# head(EX_ps_clean.rar.corr.df)
# 
# # for example, drop a column with factorial data
# EX_ps_clean.rar.corr.df <- subset(EX_ps_clean.rar.corr.df, select = -c(Treatment, Diet, Sheep_ID, ICTERIC_INDEX, LIPEMIC_INDEX))
# 
# # check header to make sure it looks good
# head(EX_ps_clean.rar.corr.df)
# 
# # check that all remaining columns are numeric instead of factor or character
# str(EX_ps_clean.rar.corr.df)

# clean up any columns which are not registering as numeric
# EX_ps_clean.rar.corr.df <- sapply(EX_ps_clean.rar.corr.df, as.numeric)
# 
# # This is from browsing the metadata
# numeric_columns <- c("Age", "Birth_weight", "Gestational_age", "House_surface", "Breastfeeding_time", "Pneumococcal_load", "Fever_time_before_sampling", "C_reactive_protein" , "Hemoglobin", "Leukocytes", "Hospitalization_days" )
# colnames(top100_decontam_rar_corr_df)
# sapply(top100_decontam_rar_corr_df, class)
# top100_decontam_rar_corr_df[,numeric_columns] <- sapply(top100_decontam_rar_corr_df[numeric_columns],as.numeric)
# sapply(top100_decontam_rar_corr_df, class)
# numeric_columns <- colnames(top100_decontam_rar_corr_df)[which(as.data.frame(sapply(top100_decontam_rar_corr_df, class))[,1] == "numeric")]
# numeric_columns