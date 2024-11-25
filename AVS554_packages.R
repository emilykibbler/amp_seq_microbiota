# This whole repo should be a proper library with dependencies
# I am aware that load_libraries() is not an ideal way to do it
# don't @ me in the github comments
# I'm doing my best and learning as I go

install_necessary <- function() {
  # Lab 1
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2")

  install.packages("tidyverse") 
  # install devtools which will allow you to install the package from GitHub
  install.packages("devtools")
  install_github("labbcb/Rqc")
  
  
  ## install these packages at your convenience ----
  BiocManager::install('phyloseq') # the data interpretation program
  BiocManager::install("microbiome") 
  install.packages("vegan") # statistical package
  install.packages("ape") # reading, writing, plotting, and manipulating phylogenetic trees
  install.packages("plyr")  # data manipulation and a dependancy of dplyr
  install.packages("dplyr") # data manipulation 
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
  
  
  # Sequence quality, lab 2
  # install devtools which will allow you to install the package from GitHub
  # install.packages("devtools") # in the install_necessary function
  install.packages("devtools")
  # more info on this package: https://www.rdocumentation.org/packages/Rqc/versions/1.6.2
  install_github("labbcb/Rqc")
  
  install.packages("readxl")
  install.packages("conover.test")
  
  install.packages("asbio")
  install.packages("microbiome")
  
  install.packages('vegan3d')
  install.packages('scatterplot3d')

  # library(vegan)

  
}

install_optional <- function() {
  
  # These packages are optional to install
  BiocManager::install("GenomeInfoDBData") # used by DADA 2 for adding species level taxonomy
  BiocManager::install("DECIPHER") # for building phylogenetic trees
  install.packages('phangorn') # for building phylogenetic trees
  install.packages("beepr") # R will beep when it is done running a command
  BiocManager::install("GenomeInfoDbData", force = TRUE)
  install.packages("ggpubr")
  install.packages("ggtext")
  install.packages("remotes")
  remotes::install_github("twbattaglia/btools")
  # if (!require(devtools)) install.packages("devtools")
  devtools::install_github("gaospecial/ggVennDiagram")
  install.packages("seqinr")
  install_github("vqv/ggbiplot")
  
  # install.packages("paletteer")
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
  
  # Load tidyverse last so it has final say in masking functions
  # Too gay to function without the tidyverse, as they say
  library(tidyverse) 
}


