queso <- list.files(getwd())
read1fns <- queso[grep("_R1.fastq", queso)]
read2fns <- queso[grep("_R2.fastq", queso)]
# 
dir.create("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw")
dir.create("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw")
file.copy(read1fns, "/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw")
file.copy(read2fns, "/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw")

length(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw"))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw"))))
length(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw")))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw"))))

ena_metadata <- read_delim("filereport_read_run_PRJEB46580_tsv.txt", delim = "\t")
tail(str_split(ena_metadata$submitted_aspera, "/"), n = 1)
ena_paths_split <- sapply(ena_metadata$submitted_aspera, str_split, "/")
accessions <- as.data.frame(ena_paths_split)[10,]
fns <- as.data.frame(ena_paths_split)[11,]
# terrible way to do this lol yolo
key <- data.frame(matrix(nrow = 93, ncol = 2))
for (i in 1:nrow(key)) {
  key[i,1] <- accessions[[i]]
  key[i,2] <- fns[[i]]
}
colnames(key) <- c("enaID","file_name")
key$sampleID <- str_remove_all(key$file_name, "_R2.fastq.gz")
key$sampleID <- str_remove_all(key$sampleID, ".R2.fastq.gz")




file_info <- metadata_import(types = c("Forward", "Reverse")) # saves a copy as an rds
files_i_got <- str_remove_all(file_info$sample_names[[1]], "_R1")
missed <- key[which(!(key$sampleID %in% files_i_got)),]

read1files <- list.files(getwd(), pattern = "_1.fastq.gz")
read2files <- list.files(getwd(), pattern = "_2.fastq.gz")

for (i in 1:length(read1files)) {
  from_name <- str_remove(read1files[i], "_1.fastq.gz")
  to_name <- key$sampleID[match(from_name, key$enaID)]
  to_name <- paste0(to_name, "_R1.fastq.gz")
  file.rename(read1files[i], to_name)
}
for (i in 1:length(read2files)) {
  from_name <- str_remove(read2files[i], "_2.fastq.gz")
  to_name <- key$sampleID[match(from_name, key$enaID)]
  to_name <- paste0(to_name, "_R2.fastq.gz")
  file.rename(read2files[i], to_name)
}

# repeat the move with the new files
queso <- list.files(getwd(), pattern = ".gz")
read1fns <- queso[grep("_R1.fastq", queso)]
read2fns <- queso[grep("_R2.fastq", queso)]
# 

file.copy(read1fns, "/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw")
file.copy(read2fns, "/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw")

length(list.files("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw"))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Forward_raw"))))
length(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw")))
length(unique(list.files(("/Users/emilykibbler/Desktop/projects/R/AVS_554/nasal/Reverse_raw"))))

match(str_remove_all(read1fns, "_R1.fastq.gz"), str_remove_all(read2fns, "_R2.fastq.gz"))
read1fns[18]
