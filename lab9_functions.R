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



# 
# 
# 
# 
# 
# 
# 
# # graph the abundance of those shared taxa, here are some example: https://microbiome.github.io/tutorials/Core.html
# 
# #aggregate the genera so we don't get a lot of lines separating all the SVs
# plot.gen <- aggregate_taxa(phylo.coreW, "Genus")
# 
# prevalences <- seq(.05, 1, .05)
# detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)
# 
# plot_core(plot.gen, 
#           plot.type = "heatmap", 
#           prevalences = prevalences, 
#           detections = detections, min.prevalence = .5) + #CHANGE min prevalence
#   xlab("Detection Threshold (Relative Abundance (%))") + ylab("Bacterial SVs") +
#   theme_minimal() + scale_fill_viridis()