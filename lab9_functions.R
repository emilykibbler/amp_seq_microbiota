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

# This is only used in the plotting script but could be useful for something I guess
# actually can this be replaced by pivot_longer?

amalgamate_SV_data <- function(input_otu_table, group_name){
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df) <- c("SV", "Group", "Abundance")
  for (i in 1:ncol(input_otu_table)) {  
    temp <- data.frame(matrix(nrow = nrow(input_otu_table), ncol = 3))
    colnames(temp) <- c("SV", "Group", "Abundance")
    temp$SV <- colnames(input_otu_table)[i]
    temp$Group <- group_name
    temp$Abundance <- as.data.frame(input_otu_table)[,i]
    df <- rbind(df, temp)
    
    
  }
  return(df)
}

# the input_df is a df of SVs considered significant by statistical tests already done
# input_phylo is what the input_df was made from

SV_abundance_df_creator <- function(input_df, input_phylo) {
  # the taxa I'll merge back into the data later
  tax_key <- as.data.frame(input_phylo@tax_table)
  tax_key$SV <- row.names(tax_key)
  df <- left_join(input_df, tax_key, by = "SV")
  
  phylo_abundance <- microbiome::transform(input_phylo, "compositional")
  abundance_df <- as.data.frame(phylo_abundance@otu_table)
  # mean abundance per SV will make it easier for plotting to sort out the more abundant ones
  abundance_means <- as.data.frame(cbind(colnames(abundance_df), 
                                         colMeans(abundance_df)))
  colnames(abundance_means) <- c("SV", "mean_abundance")
  
  abundance_df$SampleID <- row.names(abundance_df)
  abundance_df <- pivot_longer(abundance_df, 
                              cols = colnames(abundance_df)[1:length(abundance_df) - 1], # pivot everything except the last column, which we just added in: the SampleID
                              names_to = "SV",
                              values_to = "Abundance")
  abundance_df <- left_join(abundance_df, 
                             as.data.frame(input_phylo@sam_data)[,1:2],
                             by = "SampleID")
  
  df <- left_join(df, abundance_df, by = "SV")
  df <- left_join(df, abundance_means, by = "SV")
  df$mean_abundance <- as.numeric(df$mean_abundance)
  return(df)
}


