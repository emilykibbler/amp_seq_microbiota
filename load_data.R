load_lab_three <- function() {
  file_info <<- readRDS("file_info.rds") # Created in lab 2
  # errF <<- readRDS("error_profile_Forward.RDS") # created in lab 3
}


load_lab_four_partA <- function() {
  seqtab <<- readRDS("seqtab.rds") # Created in lab 3
  # FIXME make this generalizable
  filtoutput <<- readRDS("lamb_data_Forward_filtered_output.rds") # Created in lab 2
}

load_lab_four_partB <- function() {
  seqtab.nochim <<- readRDS('seqtab.nochim.rds') # Created in lab 4A
  # FIXME make this generalizable
  filtoutput <- readRDS("lamb_data_Forward_filtered_output.rds")
}


load_lab_five_partA <- function() {
  seqtab.nochim <<- readRDS('seqtab_nochim.RDS') 
  meta <<- readRDS("meta.rds")
}

load_lab_five_partB <- function() {
  meta <<- readRDS("meta.rds")
  phylo <<- readRDS("phylo.rds") # Created in 5A
}

load_lab_seven_partB <- function() {
  phylo_clean_wanted_species <<- readRDS('phylo_clean_wanted_species.rds')
  meta <<- readRDS("meta.rds")
  
  phylo_no_replace_rar <<- readRDS('clean_phylo_no_replace_rarified.rds')
  phylo_replace_rar <<- readRDS("clean_phylo_replace_rarified.rds")
}
