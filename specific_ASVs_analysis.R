## Alignments based on core SV numbers --------

library(Biostrings)
# in their supplemental material, they say ASV3 and ASV12 are in HG315101. Not really sure what that means but ok
asv3_and_asv12 <- suppressWarnings(readDNAStringSet("HG315101_strep_ASV3_ASV12.fasta", format = "fasta"))

# in their supplemental material, they say ASV2 and ASV11 are in AY612844
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

### Trying to decode paper key -------------


phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")

# this is the one with SVs as colnames, sample as row names, and read numbers as values
view(phylo_decontam_rar@otu_table)
# this is the one with SVs as row names, classifications as colnames, class names as values
view(phylo_decontam_rar@tax_table)

strep_SVs <- row.names(subset(as.data.frame(phylo_decontam_rar@tax_table), Genus == "Streptococcus"))

staph_SVs <- row.names(subset(as.data.frame(phylo_decontam_rar@tax_table), Genus == "Staphylococcus"))

strep_reads <- as.data.frame(phylo_decontam_rar@otu_table)[,strep_SVs]
strep_reads_sum <- sum(rowSums(strep_reads))
strep_reads_df <- data.frame(matrix(nrow = ncol(strep_reads), ncol = 2))
colnames(strep_reads_df) <- c("SV", "Reads")
strep_reads_df$SV <- colnames(strep_reads)
strep_reads_df$Reads <- colSums(strep_reads)
strep_reads_df$Fraction <- strep_reads_df$Reads/strep_reads_sum
strep_reads_df <- strep_reads_df[order(strep_reads_df$Fraction, decreasing = TRUE),]

subset(strep_reads_df, Fraction > 0.001) %>% 
  ggplot(aes(x = reorder(SV, Fraction, decreasing = TRUE), y = Fraction)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_blank()) +
  xlab("SV") +
  ylab("Fraction of Streptococcus Reads") +
  ggtitle("Streptococcus Reads Attributed to Individual SV",
          subtitle = "Frequency > 1 in 1000")

strep_reads_df$esk_SV <- NA
strep_reads_df$esk_SV[1] <- "strepASV1"
strep_reads_df$esk_SV[2] <- "strepASV2"
strep_reads_df$esk_SV[3] <- "strepASV3"
strep_reads_df$esk_SV[4] <- "strepASV4"
strep_reads_df$esk_SV[5] <- "strepASV5"

strep_reads_df$paper_SV <- NA
strep_reads_df$paper_SV[1] <- "ASV2"
strep_reads_df$paper_SV[2] <- "ASV3"
strep_reads_df$paper_SV[3] <- "ASV11"
strep_reads_df$paper_SV[4] <- "ASV22"
strep_reads_df$paper_SV[5] <- "ASV12"

saveRDS(strep_reads_df, "strep_reads_df.rds")

SVs_sig_diff_on_t_test <- c("GGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGCTTATGGTTAATACCCATAAGCCCTGACGTTACCCACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTATTTAAGTCAGATGTGAAAGCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGGATAACTAGAGTAGGTGAGAGGGGAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTCCCTGGCATCATACTGACACTGAGGTGCGAAAGCGTGGGTAGCAAACAG",
                            "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                            "GGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACAAATGTGTAAGTAACTATGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG")
for (i in 1:5) {
  if (strep_reads_df$SV[i] %in% SVs_sig_diff_on_t_test) {
    print(strep_reads_df$esk_SV[i])
    print(strep_reads_df$paper_SV[i])
  }
}

# 
# strep_reads$SampleID <- row.names(strep_reads)
# strep_reads <- left_join(strep_reads, as.data.frame(phylo_decontam_rar@sam_data)[,1:2], by = "SampleID")
# view(strep_reads)
# 
# # strep_reads <- t(strep_reads)
# strep_reads <- pivot_longer(strep_reads, cols = colnames(strep_reads)[1:61], names_to = "SV", values_to = "reads")

ggplot(strep_reads, aes(x = SV, y = reads)) +
  geom_bar(stat = "sum") +
  facet_wrap(~Group)

staph_reads <- as.data.frame(phylo_decontam_rar@otu_table)[,staph_SVs]


strep_reads_df <- readRDS("strep_reads_df.rds")

for (i in 1:5) {
  for (j in 1:length(paper_sig_asvs)) {
    if (length(unlist(vmatchPattern(strep_reads_df$SV[i], paper_sig_asvs[j], fixed = FALSE))) > 0) {
      print("SV sequence:")
      print(strep_reads_df$SV[i])
      print(strep_reads_df$esk_SV[i])
      print(strep_reads_df$paper_SV[i])
      print(suppressWarnings(vmatchPattern(strep_reads_df$SV[i], paper_sig_asvs[j], fixed = FALSE)))
    }
  }
}

# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "strepASV1"
# [1] "ASV2"
#$`AY612844.1 Streptococcus pseudopneumoniae strain ATCC BAA-960 16S ribosomal RNA gene, partial sequence`

# Did I find that to be significant?
"GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG" %in% SVs_sig_diff_on_t_test
# prints False, so no





# [1] "SV sequence:"
# [1] "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG"
# [1] "strepASV2"
# [1] "ASV3"
# $`HG315101.1 Streptococcus dentisani partial 16S rRNA gene, strain DSM 27088, isolate 7747`

"GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG" %in% SVs_sig_diff_on_t_test
# true
# I found it went up with antibiotics
# They found ASV3 went up with antibiotics
# I think what I am calling strepASV2 is what they are calling ASV3








