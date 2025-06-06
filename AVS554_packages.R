# FIXME: This whole repo should be a proper library with dependencies
# load_libraries() and install_x() functions are workarounds
# "Necessary" versus "Optional" are defined per Dr. Ishaq; nasal_analysis.R requires all

install_necessary <- function() {
  # Lab 1
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2")

  install.packages("tidyverse") 
  # install devtools which will allow you to install the package from GitHub
  install.packages("devtools")
  install.packages("remotes")
  # install_github("labbcb/Rqc")
  remotes::install_github("labbcb/Rqc")
  
  ## install these packages at your convenience ----
  BiocManager::install('phyloseq') # the data interpretation program
  # install.packages("microbiome")
  BiocManager::install("microbiome") 
  install.packages("vegan") # statistical package
  install.packages("ape") # reading, writing, plotting, and manipulating phylogenetic trees
  install.packages("lme4") # linear mixed models
  install.packages("lmerTest") # and permutational capacity to go along with it
  install.packages("emmeans")
  install.packages("CCA")
  BiocManager::install('DESeq2') # differential expression/abundance
  install.packages("permute") # permutational capacity
  install.packages("randomForest") # random forest feature selection
  install.packages("rfPermute") # random forest with permutational capacity
  install.packages("RColorBrewer") # graphics color customization
  install.packages("corrplot")
  BiocManager::install("decontam")
  
  install.packages("readxl")
  install.packages("conover.test")
  
  install.packages("asbio")
  install.packages('vegan3d')
  install.packages('scatterplot3d')
  
  
  BiocManager::install("msa")

}

install_optional <- function() {

  BiocManager::install("GenomeInfoDBData") # used by DADA 2 for adding species level taxonomy
  BiocManager::install("DECIPHER") # for building phylogenetic trees
  install.packages('phangorn') # for building phylogenetic trees
  install.packages("beepr") # R will beep when it is done running a command
  BiocManager::install("GenomeInfoDbData", force = TRUE)
  install.packages("ggpubr")
  install.packages("ggtext")

  remotes::install_github("twbattaglia/btools")
  # if (!require(devtools)) install.packages("devtools") # devtools is included in install_necessary
  devtools::install_github("gaospecial/ggVennDiagram")
  install.packages("seqinr")
  remotes::install_github("vqv/ggbiplot")
  
  install.packages("VennDiagram")
  install.packages("paletteer")
  install.packages("PerformanceAnalytics")
  install.packages("asbio")
  install.packages("gridExtra")
  install.packages("gtable")
  # install.packages("grid")
}


load_libraries <- function() {
  library(devtools)
  library(Rqc)
  library(BiocParallel)
  library(Biostrings)
  library(dada2); packageVersion('dada2')
  library(beepr)
  library(phyloseq)
  library(decontam)
  library(vegan)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(remotes)
  library(PerformanceAnalytics)
  library(readxl)
  library(conover.test)
  library(asbio)
  library(microbiome)
  library(DESeq2)
  library(seqinr)
  library(plyr)

  library(magrittr)

  library(knitr)
  library(randomForest)
  library(rfPermute)

  library(msa)
  
  #plotting
  library(RColorBrewer)
  library(ggpubr)
  library(ggsignif)
  library(ggtext)
  library(ggVennDiagram)
  library(RColorBrewer)
  library(corrplot)
  library(viridis)
  require(scales)
  require(grid)
  require(reshape2)
  library(ggbiplot)
  library(vegan3d)
  library(scatterplot3d)
  library(paletteer)
  library(VennDiagram)
  
  library(gridExtra)
  library(gtable)
  # library(grid)
  
  # Load tidyverse last so it has final say in masking functions
  library(tidyverse) 
}


