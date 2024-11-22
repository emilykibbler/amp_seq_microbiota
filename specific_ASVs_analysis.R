library(Biostrings)
asv3_and_asv12 <- suppressWarnings(readDNAStringSet("HG315101_strep_ASV3_ASV12.fasta", format = "fasta"))
asv2_and_asv11 <- suppressWarnings(readDNAStringSet("AY612844_ASV2_strep_ASV11_strep.fasta", format = "fasta"))
paper_sig_asvs <- suppressWarnings(c(asv3_and_asv12, asv2_and_asv11))
paper_sig_asvs

# phylo.coreW <- readRDS("phylo.coreW.rds")
atb_phylo.coreW <- readRDS("atb_phylo.coreW.rds")
no_atb_phylo.coreW <- readRDS("no_atb_phylo.coreW.rds")

atb_core_svs <- colnames(atb_phylo.coreW@otu_table)
no_atb_core_svs <- colnames(no_atb_phylo.coreW@otu_table)


# vmatchPattern(atb_core_svs[2], asv2_and_asv11, fixed = FALSE, max.mismatch = 1)
# 
# vmatchPattern(atb_core_svs[2], paper_sig_asvs[1], fixed = FALSE)
# vmatchPattern(atb_core_svs[2], paper_sig_asvs[2], fixed = FALSE)
# length(vmatchPattern(atb_core_svs[2], paper_sig_asvs[1], fixed = FALSE, max.mismatch = 1))
# length(unlist(vmatchPattern("AA", asv2_and_asv11, fixed = FALSE, max.mismatch = 2)))

for (i in 1:length(atb_core_svs)) {
  for (j in 1:length(paper_sig_asvs)) {
    if (length(unlist(vmatchPattern(atb_core_svs[i], paper_sig_asvs[j], fixed = FALSE))) > 0) {
      print("SV sequence:")
      print(atb_core_svs[i])
      print(paste("Which is column number", i))
      print(vmatchPattern(atb_core_svs[i], paper_sig_asvs[j], fixed = FALSE))
    }
  }
}
# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "Which is column number 2"
# This is ASV2 or ASV11
# $`AY612844.1 Streptococcus pseudopneumoniae strain ATCC BAA-960 16S ribosomal RNA gene, partial sequence`
# IRanges object with 1 range and 0 metadata columns:
#   start       end     width
# <integer> <integer> <integer>
#   [1]       338       762       425
dat <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(dat) <- c("sample", "group", "abundance", "asv")

temp <- data.frame(matrix(nrow = nrow(atb_phylo.coreW@otu_table), ncol = 4))
colnames(temp) <- c("sample", "group", "abundance", "asv")
temp$sample <- row.names(atb_phylo.coreW@otu_table)
temp$group <- "Antibiotics"
temp$abundance <- atb_phylo.coreW@otu_table[,2]
temp$asv <- "ASV2 or ASV11"
dat <- rbind(dat, temp)


summary(atb_phylo.coreW@otu_table[,2])
# Min.   :0.0000                                                                                                                                                                                                                                                                                                                                                                                                                           
# 1st Qu.:0.0036                                                                                                                                                                                                                                                                                                                                                                                                                           
# Median :0.0196                                                                                                                                                                                                                                                                                                                                                                                                                           
# Mean   :0.1216                                                                                                                                                                                                                                                                                                                                                                                                                           
# 3rd Qu.:0.1032                                                                                                                                                                                                                                                                                                                                                                                                                           
# Max.   :0.9682   

# 
# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "Which is column number 3"
# ASV3 or ASV12
# $`HG315101.1 Streptococcus dentisani partial 16S rRNA gene, strain DSM 27088, isolate 7747`
# IRanges object with 1 range and 0 metadata columns:
#   start       end     width
# <integer> <integer> <integer>
#   [1]       358       782       425
temp <- data.frame(matrix(nrow = nrow(atb_phylo.coreW@otu_table), ncol = 4))
colnames(temp) <- c("sample", "group", "abundance", "asv")
temp$sample <- row.names(atb_phylo.coreW@otu_table)
temp$group <- "Antibiotics"
temp$abundance <- atb_phylo.coreW@otu_table[,3]
temp$asv <- "ASV3 or ASV12"
dat <- rbind(dat, temp)

summary(atb_phylo.coreW@otu_table[,3])
# Min.   :0.0000                                                                                                                                                                                                                                                                                                                                                                                                                           
# 1st Qu.:0.0036                                                                                                                                                                                                                                                                                                                                                                                                                           
# Median :0.0196                                                                                                                                                                                                                                                                                                                                                                                                                           
# Mean   :0.1216                                                                                                                                                                                                                                                                                                                                                                                                                           
# 3rd Qu.:0.1032                                                                                                                                                                                                                                                                                                                                                                                                                           
# Max.   :0.9682   

for (i in 1:length(no_atb_core_svs)) {
  for (j in 1:length(paper_sig_asvs)) {
    if (length(unlist(vmatchPattern(atb_core_svs[i], paper_sig_asvs[j], fixed = FALSE))) > 0) {
      print("SV sequence:")
      print(atb_core_svs[i])
      print(paste("Which is column number:", i))
      print(vmatchPattern(atb_core_svs[i], paper_sig_asvs[j], fixed = FALSE))
    }
  }
}

# 

# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "Which is column number: 2"
# ASV2 or ASV11
# $`AY612844.1 Streptococcus pseudopneumoniae strain ATCC BAA-960 16S ribosomal RNA gene, partial sequence`
# IRanges object with 1 range and 0 metadata columns:
#   start       end     width
# <integer> <integer> <integer>
#   [1]       338       762       425
temp <- data.frame(matrix(nrow = nrow(no_atb_phylo.coreW@otu_table), ncol = 4))
colnames(temp) <- c("sample", "group", "abundance", "asv")
temp$sample <- row.names(no_atb_phylo.coreW@otu_table)
temp$group <- "No antibiotics"
temp$abundance <- no_atb_phylo.coreW@otu_table[,2]
temp$asv <- "ASV2 or ASV11"
dat <- rbind(dat, temp)

summary(no_atb_phylo.coreW@otu_table[,2])
# Min.   :0.0000                                                                                                                                                                                                                                                                                                                                                                                                                           
# 1st Qu.:0.0025                                                                                                                                                                                                                                                                                                                                                                                                                           
# Median :0.1123                                                                                                                                                                                                                                                                                                                                                                                                                           
# Mean   :0.2587                                                                                                                                                                                                                                                                                                                                                                                                                           
# 3rd Qu.:0.3931                                                                                                                                                                                                                                                                                                                                                                                                                           
# Max.   :0.9766  

# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "Which is column number: 3"
# ASV3 or ASV12
# 
# $`HG315101.1 Streptococcus dentisani partial 16S rRNA gene, strain DSM 27088, isolate 7747`
# IRanges object with 1 range and 0 metadata columns:
#   start       end     width
# <integer> <integer> <integer>
#   [1]       358       782       425

temp <- data.frame(matrix(nrow = nrow(no_atb_phylo.coreW@otu_table), ncol = 4))
colnames(temp) <- c("sample", "group", "abundance", "asv")
temp$sample <- row.names(no_atb_phylo.coreW@otu_table)
temp$group <- "No antibiotics"
temp$abundance <- no_atb_phylo.coreW@otu_table[,3]
temp$asv <- "ASV3 or ASV12"
dat <- rbind(dat, temp)

summary(no_atb_phylo.coreW@otu_table[,3])
# Min.   :0.00000                                                                                                                                                                                                                                                                                                                                                                                                                          
# 1st Qu.:0.00000                                                                                                                                                                                                                                                                                                                                                                                                                          
# Median :0.00000                                                                                                                                                                                                                                                                                                                                                                                                                          
# Mean   :0.02503                                                                                                                                                                                                                                                                                                                                                                                                                          
# 3rd Qu.:0.00725                                                                                                                                                                                                                                                                                                                                                                                                                          
# Max.   :0.29200     

dat %>% ggplot(aes(x = group, y = abundance, color = asv)) +
  geom_boxplot()

SVs_sig_diff_on_t_test <- c("GGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGCTTATGGTTAATACCCATAAGCCCTGACGTTACCCACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTATTTAAGTCAGATGTGAAAGCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGGATAACTAGAGTAGGTGAGAGGGGAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTCCCTGGCATCATACTGACACTGAGGTGCGAAAGCGTGGGTAGCAAACAG",
                            "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                            "GGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACAAATGTGTAAGTAACTATGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG")


for (i in 1:length(SVs_sig_diff_on_t_test)) {
  for (j in 1:length(paper_sig_asvs)) {
    if (length(unlist(vmatchPattern(SVs_sig_diff_on_t_test[i], paper_sig_asvs[j], fixed = FALSE))) > 0) { # only prints matches
      print(i)
      print("SV sequence:")
      print(SVs_sig_diff_on_t_test[i])
      print(vmatchPattern(SVs_sig_diff_on_t_test[i], paper_sig_asvs[j], fixed = FALSE))
    }
  }
}

my_top_SVs <- colnames(atb_phylo.coreW_35@otu_table)
         
for (i in 1:length(my_top_SVs)) {
  for (j in 1:length(paper_sig_asvs)) {
    if (length(unlist(vmatchPattern(my_top_SVs[i], paper_sig_asvs[j], fixed = FALSE))) > 0) { # only prints matches
      print(i)
      print("SV sequence:")
      print(my_top_SVs[i])
      print(vmatchPattern(my_top_SVs[i], paper_sig_asvs[j], fixed = FALSE))
    }
  }
}         
