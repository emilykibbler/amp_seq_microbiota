# Introduction (Lab 1) -----------------------------
# This script is compilation of approaches for processing amplicon sequencing data, and is intended for student training purposes.  
# It is tailored for 16S rRNA sequencing data, and could be modified for ITS or other genes.
# This script was compiled by Sue Ishaq for AVS 454-554 DNA Sequencing Data Analysis Lab at the University of Maine.  
# Please cite this course, and the packages used here, in manuscripts generated using this walk-through.


  ## Basic R etiquette (Lab 1) -----------------------------
# 1. To make notes and annotate code, put a hash mark (#) in front of a line or specific text you want R to ignore.
  # Like this
  print('hello') # or like this to add notes in-line


# 2. Make collapsible headings and levels by adding one or more #, then the header name, then 3+ of --- in a row.
# For example: # Level 1----; ## Level 2-----; etc
# You can show the outline/table of contents by clicking the icon of horizontal lines in the upper right of the RStudio Source Code panel

# 3. I add the word "BOOKMARK" to my code where I want to remind myself to go back later.  Then, I do a word search for it to find the places I want to revisit.

# 4. Always add an introduction to your code explaining what it is.

# 5. Always add 'housekeeping' items at the top of your code, like package installation and setting directories.

# 6. You can highlight a command to run it, or hit Control+Enter if your cursor is on that command in your code. 
# You can highlight a variable and hit enter to have R tell you what is in it, and can be a good way to see that your command worked as intended.

# 7. To remove info from your working environment and to clear up the memory that R is using to store that info:
  rm(name_of_the_thing_you_want_R_to_remove) #check the Environment tab to see what Data you have loaded that you can jettison

# save files as .RDS to save them in a R-friendly file format: 
  saveRDS('name_of_thing_to_save', 'file_name.RDS')

# 8. To set up your panels, font options, saving setting, or word-wrapping in the workflow, go to the top menu: Tools > Global Options


# Install packages into R (Lab 1) -----------------------------
# 'Packages' are expansion packs.  They are packets of code specific to certain functions, and need to be downloaded locally (to the machine running R) prior to using them.  
# You only need to install a package once, not each time you use R. Packages need to be updated periodically, and re-downloaded when you install new versions of R.
# Install packages one at a time in case there are any screen prompts. You can also "daisy chain" and install multiple packages with one command, but this can be confusing if there are errors in any one of them.
# You can see which packages you have installed by looking at the list in the "Packages" tab in one of your RStudio panels
  
  
  ## install these packages prior to the first day of class ----
 # installation manager package to help you download/install other packages, plus dada2, https://benjjneb.github.io/dada2/dada-installation.html:
source("AVS554_packages.R")
install()


# These packages are optional to install
install_optional()



# Set directory and load packages (EACH SESSION) ----------------------------------------------

# Set working directory (wd) to be the folder on your computer you want R to look for and save your files. Set this every session.
# setwd('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/') # CHANGE ME
# setwd('D:/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/')

#example syntax for a mac: 
#setwd('/Users/Sue/Desktop/')


# Load the 'packages' as 'libraries' each time you open a new R session. You can also "daisy chain" and load multiple libraries with one command.
# You can see which packages you have loaded by checking the "Packages" tab in one of your RStudio panels, and you can manually load them by checking the boxes.  You can also update some or all packages automatically.

# load these packages for the first class
library(dada2); packageVersion('dada2')
library(beepr) # if you want it to beep at you when it is done running


# these are not needed the first day of class, skip for now
library(phyloseq); packageVersion("phyloseq")
library(vegan)
library(ape); packageVersion("ape")
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(DESeq2)
library(permute)
library(randomForest)
library(rfPermute)
library(ggplot2)
library(RColorBrewer)
library(corrplot)

# optional
library(DECIPHER)
library(phangorn)

# this one is just needed for certain packages that add on to DADA2
# BiocManager::install("GenomeInfoDbData", force = TRUE) # used by DADA 2 for adding species level taxonomy
library(GenomeInfoDbData) # often R says it can't find package it just downloaded.  You can ether 1) manually check the box to enable this in the list of packages in R studio, or 2) manually enter the file location it gave you installed the package:
#library("GenomeInfoDbData", lib.loc="~/R/win-library/4.0") #CHANGE ME


# Load raw data (Lab 2) ---------------------------
# When you first get a new dataset, you should play around and find out how many sequences, samples, and sequence lengths you have, as well as look at the quality score of the dataset. This will help you plan your quality control steps. 
# This only needs to be done once, after which you will use the Filtered Data you create in Lab 1.
# Note: If you are using data from Roche 454 pyrosequencing you may want to use Pyronoise or "shh.flows" in MOTHUR to denoise your raw data and create cleaner fastq files.  This is not needed for Illumina data.



# 1. Specify the folder path to the raw (unfiltered) sequence files by assigning the folder path to a variable.  Use fastq.qz (preferred) or fastq formats.  For organization, I keep forwards (R1) and reverses (R2) in separate folders. If you have  "Undetermined" sequence files, remove this from the folder - these are sequences that your sequencer could not identify the barcode and assign it to a sample.
path_raw_F <- '/Users/emilykibbler/Desktop/projects/R/AVS_554/dat_folder_1/scallop_Fb1_raw' # 
path_raw_R <- '/Users/emilykibbler/Desktop/projects/R/AVS_554/dat_folder_1/scallop_Fb2_raw' #CHANGE ME

# example for a Mac
#path_raw_F <- '/Users/Sue/Desktop/AVS 454-554 DNA Sequencing Analysis Lab/Data/raw_forward' #CHANGE ME
#path_raw_R <- '/Users/Sue/Desktop/AVS 454-554 DNA Sequencing Analysis Lab/Data/raw_reverse' #CHANGE ME


# 2. Specify what file names (fastq files) to look for in that folder based on the type of file extension. As written, this looks for only zipped (gz) files.
fnsF <- list.files(path_raw_F)
fastqsF <- fnsF[grepl('.gz$', fnsF)] # CHANGE ME inside the '' if different file extensions, for example, might be '.fastq$'

fnsR <- list.files(path_raw_R)
fastqsR <- fnsR[grepl('.gz$', fnsR)] # CHANGE ME inside the '' if you have different file extensions

# Example if forward and reverse are in same folder which are not zipped:
#fnsALL <- list.files(path_raw_ALL)
#fastqsF <- fnsF[grepl('for.fastq$', fnsALL)] # or might look like: 'L001_R1.fastq'
#fastqsR <- fnsF[grepl('rev.fastq$', fnsALL)] # or might look like: 'L001_R2.fastq'


# 3. Fully specify the path, and file names (fns) for the reads.This will come in handy later today.
fnsF <- file.path(path_raw_F, fastqsF)
fnsR <- file.path(path_raw_R, fastqsR)

# Highlight fnsF in the previous command and hit enter to make sure that the correct info was captured.  If the variable is empty, there is and error in your code such that R couldn't do what you asked. 

# 4. Pull the sample names by reading the forward read file names, cutting it, and taking the sample name chunk.
sample.names <- sapply(strsplit(basename(fastqsF), '_for'), `[`, 1) # CHANGE ME. As written, assumes filename = samplename_XXX.fastq.gz, splits the string based on the first "_" it finds, and takes chunk 1
  # example, take the back section: sample.names <- sapply(strsplit(basename(fastqsF), '_'), `[`, 2)
  # example, take the middle section: sample.names <- sapply(strsplit(basename(fastqsF), '_', "end.fastq"), `[`, 2)


# 5. Assign the sample names to the list of Read 1 forward fastq files
names(fastqsF) <- sample.names 


  ## Check sequence quality (Lab 2) ---------------------------
# Typically in paired-end datasets, R1 (forward reads) looks great and R2 (reverse reads) less so. It may be better to use only R1 rather than merge and trim off a lot of R2.


# 1.  Randomly select a few samples to look at quality profile plot. 
# This will show the quality scores for every base in the amplicon for all the sequences. Quality score is the sequencer's confidence score in that base call, and corresponds to how many errors can be expected. >30 is good, < 20 and you should trim it out (next section). 
# How to read the plot: The green line is the average quality score for all sequences in that sample, black/grey are relative abundance of that score at each base.

plotQualityProfile(fnsF[[9]]); beep() # change the number in the brackets to look at different forward sample. Number is order of how they were loaded, not your sample number

plotQualityProfile(fnsR[[10]]); beep() # change the number in the brackets to look at different reverse sample. Number is order of how they were loaded, not your sample number


# 2. Look at the whole dataset to aggregate all samples onto one graph. This is computationally more challenging so your computer will need a few more minutes to process.

plotQualityProfile(fnsF, aggregate=T); beep() # aggregates all forward reads/read1
ggsave("fwdreads_qplot.png")

plotQualityProfile(fnsR, aggregate=T); beep() # aggregates all reverse reads/read2
ggsave("revreads_qplot.png")

# 3. Save the forward and reverse quality plots.  In the plot viewer panel, click 'export' and save as a file format of your choice.  I suggest PNG, JPEG, or PDF.



# 4. (Optional) Get a full quality report of your samples


# more info on this package: https://www.rdocumentation.org/packages/Rqc/versions/1.6.2

# install devtools which will allow you to install the package from GitHub
# install.packages("devtools")
library(devtools)
# install_github("labbcb/Rqc")
library(Rqc)


# run the command to create a quality control report, be sure to update the folder path to your raw files
rqc(path = "/Users/emilykibbler/Desktop/projects/R/AVS_554/dat_folder_1", #CHANGE ME
    sample = TRUE, n = 1e+06, group = NULL, top = 10, pair = NULL, pattern = ".fastq.gz",
    outdir = tempdir(), file = "rqc_report", openBrowser = TRUE, workers = multicoreWorkers())

# when it is done running, a webpage will open up with your sample report


# 5. Optional: search for strings in your data
library(Biostrings)

## specify one file name to look through
DNA_to_read <- readDNAStringSet("Mock_S50_ExtNeg1_for.fastq", format="fastq", with.qualities=TRUE)

# count the number of different nucleotides
alphabetFrequency(DNA_to_read)

# look for ambiguous bases
vmatchPattern("N", DNA_to_read)

# look for a homopolymer. You can change the AAAs to other letters or number of letters. you can also use this to search for indexes, primers, adaptors, or other strings
vmatchPattern("AAAAAAAAAAA", DNA_to_read, fixed=FALSE) 



    ### Summary of quality ----
# CHANGE ME: This dataset contains 34 number of samples, and 2560717 (paired forward/reverse or forward only) reads which are 300 bases long
# CHANGE ME: The R1 quality goes below 35 at cycle 200 but going to trim at 275 anyway. Yes the quality is good enough to use all R1.
# CHANGE ME: The R2 quality goes below 35 at cycle 20. No the quality not is good enough to use all R2, seqs trash along entire length

# You must include a Forward aggregate quality plot, and a Reverse aggregate quality plot if you have Revese reads, plus the Summary of Quality in an assignment on Brightspace and in your Manuscript. All other quality stats from this lecture are optional to include.



  ## Filter and trim raw sequences (Lab 2) ---------------------------
# Filer and trim sequences based on various parameters.  There are two ways to do this, on the Forward, R1 reads only (they tend to be better quality than Reverse, R2s), or you can do it on Forward/Reverse at the same time with the intent to merge them together into contigs later.

# This section uses these packages. Hash them in and load them if they aren't already loaded.
# library(dada2); packageVersion('dada2')

# 1. Create folders on your computer to place your filtered data into. You only do this once (unless you delete the folders) 
dir.create('C:/Users/sueisis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered/') #CHANGE ME
dir.create('C:/Users/sueisis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered/filtered_forward') #CHANGE ME
dir.create('C:/Users/sueisis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered/filtered_reverse') #CHANGE ME

# 2. Specify the folder path to place filtered files for reads 1 and 2.
path_filt_F <- 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered_forward' #CHANGE ME
path_filt_R <- 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered_reverse' #CHANGE ME

# 3. Choose if you want to use only the Forward R1 reads or the Forward R1 + Reverse R2, then choose the appropriate section below to Filter and Trim.

    ### Filter and trim Forward only (option 1, used)------
# Process only the Forward R1 files if the quality of your R2 was low.
# The example for processing Forward R1 and Reverse R2 together is below.
# All options can be found here: https://rdrr.io/bioc/dada2/man/filterAndTrim.html

# run this whole chunk of code from filtoutput through verbose = TRUE)
  filtoutput <- filterAndTrim(
    file.path(path_raw_F,fastqsF), # path to and file names of the raw reads
    file.path(path_filt_F, paste0(sample.names, "_F_filt.fastq.gz")), # set filtered file paths and filtered file names to create
    trimLeft=10, # cuts off the first XX bases from the F reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
    trimRight=25, # cuts of last XX bases, or hash this out and use truncLen
    #truncLen=280, # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
    #minLen=10, # optional: remove any sequences longer than max length, for use with 454 if don't want to truncate
    maxEE=2, # the maximum number of expected errors allowed in a read, always >1
    maxN= 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any!!
    rm.phix=TRUE, # remove any PhiX DNA (used as positive control),
    verbose=TRUE, matchIDs = TRUE); beep() # verbose = print to screen

# save the filtered variable for sample totals.  Your filtered files are already saved, this saves the summary R item to work with later without filtering again so you can track read numbers at each step.
saveRDS(filtoutput, "filtered_output_EXAMPLE.rds") 

# check the dimensions of the variable you created, outputs two numbers: rows (# samples), columns (# of info it added)
dim(filtoutput) 

# take a look at the counts
filtoutput

# order the info by the first column (reads.in)
filtoutput[order(filtoutput[,1],decreasing=FALSE),]

# order the info by the second column (reads.out) which are filtered reads
filtoutput[order(filtoutput[,2],decreasing=FALSE),]

# get a sum total for raw reads in and filtered reads out
colSums(filtoutput)

# look at trends
library(ggplot2)
ggplot(as.data.frame(filtoutput)) + geom_point(aes(row.names(filtoutput), reads.in),color = "blue") + geom_point(aes(row.names(filtoutput), reads.out), color = "orange")

    ### Summary of filtering ----
# CHANGE ME: Which parameters did you use and why? 
# CHANGE ME: What was the largest and small number of raw reads?
# CHANGE ME: What was the largest and small number of filtered reads?
# CHANGE ME: How many total raw reads?
# CHANGE ME: How many total filtered reads?


    ### Filter and trim Forward/Reverse (option 2, not used) -------
# Alternatively, process the Forward R1 and Reverse R2 files together if the quality of your R2 was OK. This allow you to merge F/R into 'contigs' or longer sequences later.
# The example for processing Forward R1 only is shown above.
# All options can be found here: https://rdrr.io/bioc/dada2/man/filterAndTrim.html

# specify where to send the filtered files, and what to name them
filtFs <- file.path(path_filt_F, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_filt_R, paste0(sample.names, "_R_filt.fastq.gz"))

# associate the same names with the filtered files
names(filtFs) <- sample.names
names(filtRs) <- sample.names
	  
# filter that sample pair and output as new same names, run this whole chunk from filtoutput through verbose=TRUE)	  
    filtoutput <- filterAndTrim( #filter command
      fnsF, filtFs, fnsR, filtRs, # specifies in/out variables
	  	trimLeft=c(10, 10), # cuts off the first XX bases from the F,R reads. Trim 10 for Illumina, 15 for 454 pyrosequencing.
	  	trimRight=c(25,25), # cuts of last XX bases from F,R reads, or hash this out and use truncLen
	  	#truncLen=c(330, 330), # optional: cuts off end of F reads by trimming all to same length, use instead of trimRight
	  	#maxLen=XXX, # optional: remove any sequences longer than max length, for use with 454 if don't want to truncate
	  	maxEE=c(2,3), ## The maxEE parameter sets the maximum number of expected errors allowed in a read. Always >1
	  	maxN= 0, # max number of ambiguous bases (N) allowed, DADA2 doesn't allow any!!
	  	rm.phix=TRUE, # remove any PhiX DNA (used as positive control)
	  	verbose=TRUE); beep() # verbose = print to screen
    
# save the filtered variable for sample totals.  Your filtered files are already saved, this saves the summary R item to work with later without filtering again so you can track read numbers at each step.
  saveRDS(filtoutput, "filtered_output_EXAMPLE.rds") 
    
# check the dimensions of the variable you created, outputs two numbers: rows (# samples), columns (# of info it added)
  dim(filtoutput) 
  
# take a look at the counts
  filtoutput
    
# order the info by the first column (reads.in)
  filtoutput[order(filtoutput[,1],decreasing=FALSE),]
    
# order the info by the second column (reads.out) which are filtered reads
  filtoutput[order(filtoutput[,2],decreasing=FALSE),]
    
# get a sum total for raw reads in and filtered reads out
  colSums(filtoutput)
    
    ### Summary of filtering ----
# CHANGE ME: Which parameters did you use and why? 
# CHANGE ME: What was the largest and small number of raw reads?
# CHANGE ME: What was the largest and small number of filtered reads?
# CHANGE ME: How many total raw reads?
# CHANGE ME: How many total filtered reads?
  


# DADA2 learn error rates (Lab 3) ----------------------------
# DADA2 uses 4 methods to estimate base call error rates, using a random subset of bases from the data, DADA2 pipeline recommends learning the error rates prior to dereplication (combining identical sequences into one).
# You should do this for each sequencing run separately, because some error is inherent to the run.  You can combine sequencing runs later.

# This section uses these packages. Hash them in and load them if they aren't already loaded.
# library(dada2); packageVersion('dada2')

# 1. Set the seed for using random number generators, so it is more reproducible because you start in the same place.
  set.seed(1234) 

# 2. Specify the folder path (or paths if you used R1 + R2 to Filter and Trim) to your filtered files.
  path_filt_b1 <- '/Users/emilykibbler/Desktop/projects/R/AVS_554/scallop_data/filtered_b1'
#  path_filt_R <- 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered_reverse'

# 3. Set the names for Forward R1 filtered files. If your file names have a different patten to the ones used in the example, you may need to update the strsplit portion.
  fnsB1 <- list.files(path_filt_b1, full.names = TRUE)
  filtsB1 <- fnsB1[grepl('.gz$', fnsB1)] # CHANGE if different file extensions
  sample.namesB1 <- sapply(strsplit(basename(filtsB1), '_B1_filt'), `[`, 1) # Assumes filename = samplename_XXX.fastq
  names(filtsB1) <- sample.namesB1

# Set the  names for Reverse R2 filtered files if you created them.
  fnsR <- list.files(path_filt_R, full.names = TRUE)
  filtsR <- fnsR[grepl('.gz$', fnsR)] # CHANGE if different file extensions
  sample.namesR <- sapply(strsplit(basename(filtsR), '_R_filt'), `[`, 1) # Assumes filename = samplename_XXX.fastq
  names(filtsR) <- sample.namesR


# 4. Learn the error rates from your sequencing run. Some error is inherent to sequencing run, so this should be done on each run separately, and you can combine multiple datasets later at the sequence table step. Be sure to include controls or low quality samples for learning error. When using smaller datasets, DADA2 recommended nbases=1e6 reads, but for the big data workflow, it recommended nbases=2e6. Multithread = TRUE if you are using a mac that can handle it.
  errB1 <- learnErrors(filtsB1, nbases = 1e6, multithread=TRUE, randomize=TRUE); beep() #I can use multithread=TRUE because I'm on a mac
  # errR <- learnErrors(filtsR, nbases = 1e6, multithread=FALSE, randomize=TRUE)

# 5. Save the error profiles to output files
  saveRDS(errB1, paste0("/Users/emilykibbler/Desktop/projects/R/AVS_554/scallop_data/", 'error_profile_B1.RDS')) #CHANGE ME so 'example' is the name of your dataset
 # saveRDS(errR, 'error_profile_R_EXAMPLE.RDS') #CHANGE ME so 'example' is the name of your dataset

# 6. Optional, look at error plots to view the prediction rates of different sequencing errors.
# From the function info page: "This function plots the observed frequency of each transition (eg. A->C) as a function of the associated quality score. It also plots the final estimated error rates (if they exist). The initial input rates and the expected error rates under the nominal definition of quality scores can also be shown."
# Thread on plot interpretation: https://github.com/benjjneb/dada2/issues/975 
  plotErrors(errB1, nominalQ=TRUE) +
    ggtitle("Scallop, batch 1, filtered")
  ggsave("/Users/emilykibbler/Desktop/projects/R/AVS_554/scallop_data/batch1_errors.png")
 # plotErrors(errR, nominalQ=TRUE)


 
  ## DADA2 pick sequence variants (Lab 3) ----------------------
# Do not pool different sequencing runs prior to this step as they have their own inherent error rates. Can set the error rate as a matrix that has already been learned, or learn now as selfConsist=TRUE

# Reload as needed if you are picking this up after closing your session.
# errF <- readRDS('error_profile_F_EXAMPLE.RDS')
# errR <- readRDS('error_profile_R_EXAMPLE.RDS')
 
# Reload as needed if you are picking this up after closing your session.
# path_filt_F <- 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered_forward'
# path_filt_R <- 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/filtered_reverse'
 
# Reload as needed if you are picking this up after closing your session. If your file names have a different patten to the ones used in the example, you may need to update the strsplit portion.
 # fnsF <- list.files(path_filt_F, full.names = TRUE)
 # filtsF <- fnsF[grepl('.gz$', fnsF)] # CHANGE if different file extensions
 # sample.namesF <- sapply(strsplit(basename(filtsF), '_filt'), `[`, 1) # Assumes filename = samplename_XXX.fastq
 # names(filtsF) <- sample.namesF 
  
 # fnsR <- list.files(path_filt_R, full.names = TRUE)
 # filtsR <- fnsR[grepl('.gz$', fnsR)] # CHANGE if different file extensions
 # sample.namesR <- sapply(strsplit(basename(filtsR), '_'), `[`, 1) # Assumes filename = samplename_XXX.fastq
 # names(filtsR) <- sample.namesR
 
# 1. Decide if you are running the Forward R1 only or the Forward R1 + Reverse R2, and choose the appropriate section below.

    ### Sample inference of FORWARD reads only or F/R -----------------------  
# This step will dereplicate data (combine identical sequences while keeping count of them), and use the error rates to decide which is a sequence variant and which sequences should be combined because they are a few base errors away from being identical. 
  
  
# 1. Sample inference of FORWARD reads only using previously made error profile. This may take an hour if you have a slow computer and/or a lot of data.
  dadaFs <- dada(filtsF, err=errF, multithread=FALSE); beep()

# 1.1 OPTIONAL: Sample inference of REVERSE reads only. Do step 1 and 1.1 if you plan to merge F/R reads.
#  dadaRs <- dada(filtsR, err=errR, multithread=TRUE); beep()

# 2. Construct a sequence table. DADA2 may drop samples that didn't have any sequences which passed the quality control steps.
 seqtab <- makeSequenceTable(dadaFs); beep() #if you used F only
 
# 2.1 or construct this one if you are using F/R
#  mergers <- mergePairs(dadaFs, filtsF, dadaRs, filtsR, verbose=TRUE); beep()
#  seqtab <- makeSequenceTable(mergers)
  
# 3. Get the dimensions of your table, which is a great way to count how many samples are and SVs included.
 dim(seqtab) 
 # CHANGE ME: 34 # of samples (first number, row count) and 13,483# of SVs (second number, column count)
 
# 4. Save the seqtab before you quit!
 saveRDS(seqtab, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/seqtab_EXAMPLE.rds') #CHANGE ME
 
# DADA2 Remove chimeras from seqtab (Lab 4) ----------------------------	
# Chimeras are accidentally created sequences during lab protocols. Remove them.	

# This section uses these packages. Hash them in and load them if they aren't already loaded.
# library(dada2); packageVersion('dada2')

# Reload your seqtab if you are starting a new session. You don't need any other previous variables or data.
seqtab <- readRDS('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/seqtab_EXAMPLE.rds') #CHANGE ME
	
	
# 1. Remove chimeras. Leave the method as consensus. multithread processing can be used with mac or linux, and verbose means it will print results to the screen		
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE) 
# CHANGE ME: My samples have 1519 bimeras out of 13483 total sequences

# 2. Check dimensions of your cleaned sequence table to see how many samples and sequences now remain.
dim(seqtab.nochim) 
# CHANGE ME: 34 # of samples (first number, row count) and 11964# of SVs (second number, column count)

# Calculate the percentage of chimeras identified out of the total
sum(seqtab.nochim)/sum(seqtab)

# 3. Save your chimera-free sequence table.
	saveRDS(seqtab.nochim, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/seqtab.nochim_EXAMPLE.rds') #CHANGE ME
	
	#reload as needed
	seqtab.nochim <- readRDS('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/seqtab.nochim_EXAMPLE.rds')
	

## Workflow verification steps  (Lab 4) --------------------------
# This section can be optional, but as you are learning, or playing with new data, it is helpful to run some internal checks to assess whether you like how your quality control steps worked.

# This section uses these packages. Hash them in and load them if they aren't already loaded.
require(tidyr)
library(ggplot2)
		
# 1. Track reads through the analysis here to see how many were lost at each QC step, in case you were too stringent. Only need to count F reads, even if you did paired F/R
getN <- function(x) sum(getUniques(x)) #set function to count the unique reads
	
# load the filtered output file if it isn't already in your environment. This should be the file you created in the filterAndTrim step
filtoutput <- readRDS("filtoutput_EXAMPLE.rds") #CHANGE file name
	
# 2. Load the metadata file that goes with your sequencing data so you can match factors to seq data. This uses a datatable made in Excel and saved as a .csv file
meta <- read.csv('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/Particle_Size_metadata.csv', # CHANGE file name
                 header = TRUE, # specifies that the top row is column names
                 row.names = 1) # row.names specifies which column to call the samples by (sample name), doesn't have to be the first column


# 2.5 check the dimensions of the three data files you need for this to make sure the number of rows matches in each. If they do not, you may need to add/remove rows from your metadata file in case samples were removed/retained from your dataset.
dim(filtoutput) # 34 rows
dim(seqtab.nochim) # 34 rows
dim(meta) # 34 rows

# If they don't match, check that you aren't missing any samples you should have. If you have too many in the meta file because some samples failed sequencing, you can select which ones to keep based on matching to the list of names in your seqtab

meta_cleaned <- meta[row.names(meta) %in% row.names(seqtab.nochim), ]

# 3. Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
track <- cbind(filtoutput, rowSums(seqtab.nochim), meta$Sample_type) #CHANGE Sample_type to be the column name in your metadata file of the groups in your data

# 4. Assign column names to that new variable
colnames(track) <- c("reads.in","filtered", "nonchimeras", "Sample_type") # CHANGE "Sample_type" to be the column name in your metadata file of the groups in your data

# 5. Assign rownames for the samples in your new variable
rownames(track) <- rownames(meta)

# 6. Look at the header of that variable to check that it looks the way you want it.
head(track)

# 7. Save the tracking variable and data as an R file
saveRDS(track, 'tracked_seqs_EXAMPLE.RDS') #CHANGE file name to whatever you want

# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
plotData <- as.data.frame(track) %>% gather(type,totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))
                               
                               
# plot with Sample_type along the X axis
ggplot(plotData,aes(x=Sample_type,y=as.numeric(totals))) + geom_jitter(aes(color=type)) + #CHANGE Sample_type to be the factor you previously chose
  ylab("Sequences") + xlab("Sample_type") 
  
  # or, plot with QA stage along the X axis
  ggplot(plotData,aes(x=type,y=as.numeric(totals))) + geom_jitter(aes(color=Sample_type)) + #CHANGE Sample_type to be the factor you previously chose
    ylab("Sequences") + xlab("QA stage") 

  
# save this plot  
  
  
  
  

# 9. Randomly select a negative control sample to see what's in it
unqs.NC <- seqtab.nochim["Mock_S13_PCRneg",] #CHANGE ME, update " " with sample name of a negative control, such as Mock_S13_PCRneg

  # Sort by # seqs and drop SVs absent in the negative control sample
  unqs.NC <- sort(unqs.NC[unqs.NC>0], decreasing=TRUE) 

  # Print out how many reads are inferred in the negative control
  cat("DADA2 inferred", length(unqs.NC), "sequence variants present in the selected sample.\n")
	# CHANGE ME: There are 81 SVs in this PCR negative control

  # Plot the number of sequences in the SVs found in the negative control sample.
  plot(unqs.NC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Negative Control") #CHANGE ME

  # OPTIONAL: add taxonomy and see what's in the negative control, for example if you want to know if you have the same contaminant showing up in negative controls over time.
  taxa.NC <- assignTaxonomy(unqs.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
  write.csv(taxa.NC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_NC_EXAMPLE.csv') 

	
  
  
# 10. Randomly select a positive control to see what's in it
unqs.PC <- seqtab.nochim["Mock_S192",] #CHANGE ME, update with sample name of a positive control, such as Mock_S192

  # Sort by # seqs and drop SVs absent in the positive control sample
  unqs.PC <- sort(unqs.PC[unqs.PC>0], decreasing=TRUE) 

  # Print out how many reads are inferred in the positive control
  cat("DADA2 inferred", length(unqs.PC), "sequence variants present in the selected sample.\n")
  # CHANGE ME: There are 34 SVs in this positive control. We put 32 taxa into the mock community for this sample.

  # Plot the number of sequences in the SVs found in the positive control
  plot(unqs.PC, ylab="Number of seqs/SV, Neg Control", xlab="SVs", main="Number of sequences/SV in Positive Control") #CHANGE ME
	
# OPTIONAL: add taxonomy and see what's in the positive control, for example if you want to know if you only got back what you put in.
  taxa.PC <- assignTaxonomy(unqs.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz', minBoot = 75)
  write.csv(taxa.PC, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_PC_EXAMPLE.csv') 
	

  
  

# Phyloseq First look (Lab 5) ------------------------------
library(phyloseq)
  
# load seqtab as needed
seqtab.nochim <- readRDS('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/seqtab-nochim-EXAMPLE.RDS') #CHANGE ME


# load metadata table as needed
meta <- read.csv("C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/Particle_Size_metadata.csv", header=T, row.names=1) #CHANGE ME


# what to do if my sample names don't match. check the sample sames and see if they match
row.names(meta)
row.names(seqtab.nochim)

# OPTIONAL, if they don't match: rerun the code after changing the sample.names steps from earlier, or fix the names in your metadata file to match what is in the seqtab, or use one set of names to rename the other (make sure they are in the same order before you do this)
   meta <- meta[order(row.names(meta)),] # this orders them alphabetically
   seqtab.nochim <- seqtab.nochim[order(row.names(seqtab.nochim)),] # this orders them alphabetically
   row.names(seqtab.nochim) <- row.names(meta) # this copies and pastes the names
   saveRDS(seqtab.nochim, 'seqtab-nochim-EXAMPLE.RDS') # save so you don't have to do it again

## create a phyloseq object with all samples (subset out later)
EX_ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta))

EX_ps
# how many samples made it this far? CHANGE ME, 34 samples, 30 variables, 11964 taxa



# If you have sequenced controls, it's a good idea to look at your data in comparison to them.

library(ggplot2)

# Alpha diversity peek 
 plot_richness(EX_ps, x="Sample_type", # CHANGE the x-axis to a factor of your choice
               measures=c("Observed","Chao1", "Shannon"), # these are some of the alpha diversity measures you can use 
               color="Treatment") + # CHANGE the color to a factor of your choice
               theme_bw() # a ggplot theme to make the graph look nice
 
## Note: hopefully, diversity/richness is lower in the Negative controls than real samples
## SUMMARY: are the negative controls much higher or lower than your real samples? 
 
 
 
# Plot the taxa sums to see how populated each taxa is (do you have many rare taxa?
	plot(sort(taxa_sums(EX_ps), TRUE), 
	     type="h", 
	     ylim=c(0, 20)) #limit the y-axis to better see the long tail of rare taxa

	# optional to save a copy of this figure
	
	
# Create a simple ordination to look for clustering by extraction batch or confounding variables in your metadata 
 EX.ord <- ordinate(EX_ps, #calculate similarities
                    method ="PCoA", #ordination type
                    "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)
 
 
 plot_ordination(EX_ps, EX.ord, type="samples", color="Sample_type")
 # Initial clustering by extraction/sequencing batch or confounding variables implies contamination issues (see next section)
 # Horse-shoe patterns indicate underlying patterns to the data

  # save this graph for later
 
 
# Removing contamination using negative controls (Lab 5)-----------------------------------
 ### This section can only be accomplished if you have sequenced negative controls and/or DNA quantification data from when the pool was created.  If you have these, CHOOSE ONE method to follow.  If you don't, skip this section and consider including sequencing controls in your next run.
 
 # Dr. Ishaq's method: https://github.com/SueIshaq/Examples-DADA2-Phyloseq
 # decontam in R: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
 
 ## ------ Dr Ishaq's method -------------
 # Dr. Ishaq's method creates vectors out of the SV table data for negative controls, and subtracts those SVs from the sample data.  Depending on the type of negative control, these are removed from the whole data set or from subsets of batches. Remove PCR and sampling materials negative control SVs fully from all samples, and remove extraction kit SVs fully from each dna_extraction_batch, respectively.
 
 # 1, Subset the controls that would pertain to all of the samples, like ones from reagents, PCR, or other whole-batch effects. prune to be only those taxa	
  EX_ps_controls = subset_samples(EX_ps, Sample_type == "Mock" | Sample_type == "NegCon_PCR") # CHANGE ME to the column name that holds the variables associated with being a negative control
 
  EX_ps_controls <- prune_taxa(taxa_sums(EX_ps_controls) > 0, EX_ps_controls)
 #6368 taxa in the 4 mock samples
 
 # 2, Make the taxa names into a vector so you can remove them. 
 control_vec <- as.vector(taxa_names(EX_ps_controls))
 vec <- as.vector(taxa_names(EX_ps))
 keep <- setdiff(vec, control_vec)
 
 # 3. Use the keep vector for the prune taxa argument, because it wants the argument to be true (matching).
 EX_ps_NCclean <- prune_taxa(keep, EX_ps)
 
 # 4. Then remove the samples which are now empty, namely the NegCon_PCR and NegCon_swab.
 EX_ps_NCclean <- prune_samples(sample_sums(EX_ps_NCclean) > 0, EX_ps_NCclean)
 
 
 # NOTE: If you have multiple batches of DNA extractions, you'll want to copy the code from steps 5,6,7 below and have a set for each batch. If you don't have negative controls for each DNA extraction, jump to step 10.
 
 ## 5. subset out each dna_extraction_batch 
 batch1 = subset_samples(EX_ps_NCclean, DNA_extraction_batch == "1") #16 samples
 batch2 = subset_samples(EX_ps_NCclean, DNA_extraction_batch == "2") #16 samples
 
 
 ## 6. subset the controls and prune to be only those taxa
 batch1_kit = subset_samples(batch1, Sample_type == "NegCon_kit")
 batch2_kit = subset_samples(batch2, Sample_type == "NegCon_kit")
 
 batch1_kit <- prune_taxa(taxa_sums(batch1_kit) > 0, batch1_kit) # 1 sample
 batch2_kit <- prune_taxa(taxa_sums(batch2_kit) > 0, batch2_kit) # 1 sample
 
 
 # 7. make the taxa names into a vector so you can remove them 
 b1_control_vec <- as.vector(taxa_names(batch1_kit))
 b1_vec <- as.vector(taxa_names(batch1))
 b1_keep <- setdiff(b1_vec, b1_control_vec)
 ## then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching)
 b1_clean <- prune_taxa(b1_keep, batch1)
 
 # make the taxa names into a vector so you can remove them 
 b2_control_vec <- as.vector(taxa_names(batch2_kit))
 b2_vec <- as.vector(taxa_names(batch2))
 b2_keep <- setdiff(b2_vec, b2_control_vec)
 ## then, use the keep vector for the prune taxa argument, because it wants the argument to be true (matching)
 b2_clean <- prune_taxa(b2_keep, batch2)
 

 
 
 # 8. Merge the phyloseq objects back together, then remove any blank taxa or samples
 EX_ps_NC_batch_clean <- merge_phyloseq(b1_clean, b2_clean)
 
 ## 9. clean out taxa/SV columns that are no longer present
 EX_ps_NC_batch_clean <- prune_taxa(taxa_sums(EX_ps_NC_batch_clean) > 0, EX_ps_NC_batch_clean)
 EX_ps_NC_batch_clean <- prune_samples(sample_sums(EX_ps_NC_batch_clean) > 0, EX_ps_NC_batch_clean)
 # 30 number of samples and 11793 SVs remaining CHANGE ME
 
 # 10. Did it work? Check your ordination again
 EX_cleaner.ord <- ordinate(EX_ps_NC_batch_clean, #calculate similarities
                            method ="PCoA", #ordination type
                            "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)
 
 plot_ordination(EX_ps_NC_batch_clean, EX_cleaner.ord, type="samples", color="DNA_extraction_batch", title="After cleaning out negative controls")
 
 # save this graph for later
 
  ### Summary of decontamination ----
 # CHANGE ME: Which types of negative controls did you use and why? 
 # CHANGE ME: Did you notice clustering by batch before and after removing contaminants?

 
 
 ## ------ Decontam method -------------
 
 if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("decontam")
 
 library(decontam) # This packages identify contaminants by frequency of SVs.
 
 # The first contaminant identification method uses the distribution of the frequency of each sequence as a function of the input DNA concentration. Essentially, is it too rare to be real?
 
 # 1. Look for contaminants using DNA quantification data from your metadata file
 contamdf.freq <- isContaminant(EX_ps, method="frequency", conc="quant_reading") #CHANGE conc= "quant_reading" to the column heading that holds the concentration information in your metadata
 
 # 2. See what it looks like 
 head(contamdf.freq) 
 
 # 3. Make it into a table
 table(contamdf.freq$contaminant)
 
 # 4. See which SVs are categorized as contaminants based on frequency
 head(which(contamdf.freq$contaminant)) 
 
 # 5. Get rid of the contaminants 
 EX_ps.noncontamfreq <- prune_taxa(!contamdf.freq$contaminant, EX_ps)
 
 # 6. How much is left
 EX_ps.noncontamfreg
 # CHANGE ME: number of samples and SVs left 
 
 
 
 # The second contaminant identification method uses prevalence (presence/absence) across samples as compared to the prevalence in negative controls. Essentially, is it in most of your samples but not most of your controls?
 
 
 # 7. Identify negative controls by indicating which column/factor in your metadata and which variable indicate a negative control
 sample_data(EX_ps)$is.neg <- sample_data(EX_ps)$Sample_or_Control == "Control Sample"
 
 # 8. Calculate prevalence of SVs in samples versus controls
 contamdf.prev05 <- isContaminant(EX_ps, method="prevalence", neg="is.neg", threshold=0.5)
 
 # 9. Make a table
 table(contamdf.prev05$contaminant)
 
 # 10. Look at it
 head(which(contamdf.prev05$contaminant))
 
 # 11. get rid of the contaminants 
 EX_ps_NC_decontam_clean <- prune_taxa(!contamdf.prev05$contaminant, EX_ps.noncontamfreq) # remove them from the original phyloseq object, or the one with cleaned out SVs by frequency
 
 # 12. Check how many are left
 EX_ps_NC_decontam_clean
 # CHANGE ME: number of samples and SVs left
 
 
 # 13. Did it work? Check your ordination again
 EX_cleaner.ord <- ordinate(EX_ps_NC_decontam_clean, #calculate similarities
                            method ="PCoA", #ordination type
                            "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)
 
 plot_ordination(EX_ps_NC_decontam_clean, EX_cleaner.ord, type="samples", color="dna_extraction_batch", title="After cleaning out negative controls")
 
 # save this graph for later
 
  ### Summary of decontamination ----
 # CHANGE ME: Which types of negative controls did you use and why? 
 # CHANGE ME: Did you notice clustering by batch before and after removing contaminants?
 
 
 
 
 # DADA2 assigning taxonomy from a reference database file (Lab 6) ---------------------------------
 
 
 # Be sure to have a taxonomy database file downloaded from here: https://benjjneb.github.io/dada2/training.html. 
 
 # The Silva file shown below is used for 16S rRNA (prokaryotes) and nicely formatted versions can be downloaded from the DADA2 website, which also contains some (18S) eukaryotic versions.  Let me know if you need a custom one for your dataset.
 
 ## Option 1: If you were not able to use any of the decontamination steps, be sure to use the seqtab.nochim variable, and reload it if you need to.----
  # 1. reload as needed
 seqtab.nochim <- readRDS('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/seqtab.nochim_EXAMPLE.rds')
 
 #2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops. If you have a lot of samples, or a slow laptop, ask me about subsetting this into smaller chunks.
 
 all.taxa <- assignTaxonomy(seqtab.nochim, 
                            'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz',        # CHANGE file path
                            tryRC = TRUE, # didn't get any IDs the first time? Use this to try the reverse complement of your sequences instead
                            minBoot = 75, verbose = TRUE); beep() 
 
 # 3. saveRDS(all.taxa, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_EXAMPLE.rds') #CHANGE ME
 write.csv(all.taxa, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_EXAMPLE.csv') #CHANGE ME
 
 
 
 ## Option 2: If you were able to run a decontamination steps, you'll need to grab the cleaned sequence table from your phyloseq object. -----

 # 1. take the SVs (otu table) from the phyloseq object and make it into a dataframe again
 seqtab.nochim.decontam.sd = as(otu_table(EX_ps_NC_batch_clean), "matrix")

  seqtab.nochim.decontam = as.matrix(seqtab.nochim.decontam)
  
 ## 2. assign taxonomy. this may be memory intensive, and may take a few hours on slow laptops.
 all.taxa <- assignTaxonomy(seqtab.nochim.decontam, 
                            'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_nr99_v138_train_set.fa.gz',        # CHANGE file path
                            tryRC = TRUE, # didn't get any IDs the first time? Use this to try the reverse complement of your sequences instead
                            minBoot = 75, verbose = TRUE); beep() 
 
 
 
# 3. saveRDS(all.taxa, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_EXAMPLE.rds') #CHANGE ME
 write.csv(all.taxa, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_silva_EXAMPLE.csv') #CHANGE ME
 
 
 
 # OPTIONAL: try adding species designation to the table. this may be memory intensive, and may take a few hours on slow laptops.You can also try this after the Lab 6 step to remove certain taxa by name, in case you need to save on computational space.
 
 all.taxa.sp <- addSpecies(all.taxa, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/silva_species_assignment_v138.fa.gz', allowMultiple = FALSE, verbose = FALSE); beep() #CHANGE file path
 
 saveRDS(all.taxa.sp, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_species_silva_EXAMPLE.rds') #CHANGE ME
 write.csv(all.taxa.sp, 'C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_species_silva_EXAMPLE.csv') #CHANGE ME
 

 
 
 
# 4. Remake the phyloseq object with the new taxonomy file 
 
 # reload metadata table as needed
 meta <- read.csv("C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/Particle_Size_metadata.csv", header=T, row.names=1) #CHANGE ME
 
 #reload taxa table as needed
 all.taxa.sp <- readRDS('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/Data/taxa_species_silva_EXAMPLE.rds')
 
 otu_t <- otu_table(EX_ps_NC_batch_clean)
 
 ## create a phyloseq object with all samples
 EX_ps_NC_batch_clean <- phyloseq(otu_table(otu_t, taxa_are_rows=FALSE), 
                   sample_data(meta),
                   tax_table(all.taxa.sp))
 
 
 saveRDS(EX_ps_NC_batch_clean, "Example_phyloseq_object_cleaned.rds")
 
 
 
 
# Clean out unwanted taxa (Lab 6) ------------------------------------
# Regardless of whether you included sequenced negative controls, you can remove taxa which are of no interest to you

# clean out chloroplast and mitochondria. can also elect to remove other contaminating domains or kingdoms as needed. 
#== it means to keep that, or set as equal to it. if you have != it means not that (not matching). & adds more sections to this.
  require(dplyr)
 
 # Optional: explore your taxonomy before filtering. Use the tax table you made
 df <- as.data.frame(EX_ps)
 table(df$Kingdom)
 table(df$Phylum)
 table(df$Class)
 table(df$Order)
 
 
  EX_ps_clean <- EX_ps_NC_batch_clean %>% # CHANGE ME to your phyloseq object name 
   subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast") # CHANGE ME to taxa you want to remove
 

## clean out taxa/SV columns that are no longer present
  EX_ps_clean <- prune_taxa(taxa_sums(EX_ps_clean) > 0, EX_ps_clean)
  EX_ps_clean <- prune_samples(sample_sums(EX_ps_clean) > 0, EX_ps_clean)
  
  EX_ps_clean
  # CHANGE ME: 30 number of samples and 8091 SVs left 


# Save your clean phyloseq object
saveRDS(EX_ps_clean, 'EX_ps_clean_phyloseq_object.RDS') #CHANGE ME

# Reload as needed
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS') 








# Write out a description of experimental design (Homework)

# Rarefaction (Lab 7) --------------------------------------------------
# To compare sequence data accurately, it is often necessary to rarefy/normalize SVs to even depth across all samples
# Rarefaction is not subject to library effect sizes, and reportedly works best (compared to logUQ, CSS, DESeqVS, edgeR-TMM): https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
 
# Reload as needed
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS') 


# make a rarefaction curve to see if your samples have enough coverage. To make it prettier, check out this tutorial: https://fromthebottomoftheheap.net/2015/04/16/drawing-rarefaction-curves-with-custom-colours/
library(vegan)
EX-rarec <- rarecurve(otu_table(EX_ps_clean), step = 10, cex=0.5, label = TRUE) # step is major tick marks on x axis (x 1000), cex is text size, label is sample name
# optional to save this plot

# take a look at rowsums, or total sequences per sample
sort(rowSums(otu_table(EX_ps_clean)))
# CHANGE ME: smallest reads (sequences) in a sample ______
# CHANGE ME: largest in a sample ______
# CHANGE ME: number of samples with <5000 ______

EX_ps_clean.rar <- rarefy_even_depth(EX_ps_clean, 
                                     sample.size=3592, # CHANGE ME to the sequences/sample you want. 5-10k is a good amount, more is better
                                     replace=FALSE, #sampling with or without replacement
                                     trimOTUs=TRUE, #remove SVs left empty (called OTUs here but really they are SVs) 
                                     rngseed=711, 
                                     verbose=TRUE)

# CHANGE ME: how many 0 samples and 23932 SVs were removed

# Save your cleaned, rarefied phyloseq object  
saveRDS(EX_ps_clean.rar, 'EX_ps_clean_rarefied_phyloseq_object.RDS')

# Helpful to have an SV table from the clean, rarefied phyloseq
write.csv(otu_table(EX_ps_clean.rar), 'EX_ps_clean_rarefied_phyloseq_object.csv')





# Pick up the statistical analysis here (Lab 7+) -------------------------------------
# Once you are satisfied with your QC, you can start from here with your phyloseq object as needed

library(phyloseq)

# Non-rarefied version
EX_ps_clean <-readRDS('EX_ps_clean_phyloseq_object.RDS')

# Rarefied version
EX_ps_clean.rar <-readRDS('EX_ps_clean_rarefied_phyloseq_object.RDS')


# if you have calendar dates that you want to be treated as time:
require(lubridate)

# use as.Date to import date into a usable format for continuous variable. #"%m/%d/%y" is the format of how your date is written in the metadata, and $Sample_date is the column of your time data. Repeat as needed for any dates you have.
sample_data(EX_ps_clean.rar)$Sample_date <- as.Date(sample_data(EX_ps_clean.rar)$Sample_date, "%m.%d.%Y") # CHANGE ME


# Alpha diversity (Lab 7 and 8) ---------------------------------
  ## Alpha diversity graphics (Lab 7) ----------------------

# load ggplot and other packages to make pretty graphs
library(ggplot2)
library(RColorBrewer)

# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
sample_data(EX_ps_clean.rar)$Diet <- factor(sample_data(EX_ps_clean.rar)$Diet, levels=c("LooseAlfalfa", "PelletedAlfalfa"), labels=c("Loose Alfalfa", "Pelleted Alfalfa")) #CHANGE ME

# plot alpha diversity with phyloseq: https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/plot_richness. 
# measures include c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

plot_richness(EX_ps_clean.rar, 
              x="Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures="Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x=element_blank()) + #example, get rid of x axis labels
  geom_violin(trim=TRUE, aes(fill=Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  facet_grid(.~Week, space="free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  xlab("Week")


# Use geoms to add different shapes/styles of plots: https://ggplot2-book.org/individual-geoms.html#:~:text=3.1%20Basic%20plot%20types,plot%20has%20a%20special%20name.

# Use multiple geoms to create gaph overlays in a single plot, and rearrange their order to move to the foreground or background.

# want to map factors to geoms to get it to set colors by factor? use geom_point(aes(color=Factor)). otherwise, use geom_point(color=Black) to make them all the same.


# two different richness plots grouped
library(ggplot2)
library(ggpubr)

plot1 <- plot_richness(EX_ps_clean.rar, 
                       x="Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures="Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)")

plot2 <- plot_richness(EX_ps_clean.rar, 
                       x="Treatment", #CHANGE ME, A is whatever factor you want on x-axis
                       measures="Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)")


# plot together
ggarrange(plot1,
          ggarrange(plot2, nrow=1, labels = c("B")),
          ncol = 2, labels ="A", common.legend = TRUE)




# richness plot with significance added
library(ggplot2)
library(ggsignif)


# richness plot observed SVs with lines to fit the view screen
plot_richness(EX_ps_clean.rar, 
              x="Treatment", #CHANGE ME, A is whatever factor you want on x-axis
              measures="Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim=TRUE, aes(fill=Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group=Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") + 
  ylim(0,1500) + #define the y axis min/max
  geom_segment(aes(x = 1, y = 1200, xend = 2, yend = 1200)) +  geom_text(x = 1.5, y = 1250, label = "***") # add a drawn in line and significance tag, adjusting the x and y coordinates till it fits where you want in the window.  Add another for each line to add. As written, this will fit the view window you have, if you adjust that your segments will not adjust with it.




  ## Alpha diversity plotted against other metadata -----

# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon")) #change to whatever measures you want

# OPTIONAL: use phyloseq to calculate Faith's Diversity metric, https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
install.packages("remotes")
remotes::install_github("twbattaglia/btools")
EX_faith <- estimate_pd(EX_ps_clean.rar)

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)

# Coerce to data.frame and add the metadata for these samples
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

dim(EX_ps_clean.rar.rich.df)

# make a graph using that dataframe. this is a generic example, you will want to personalize this.
ggplot(data=EX_ps_clean.rar.rich.df, aes(x=Diet, y=Observed)) + 
    theme_minimal() + 
    geom_point(aes(color=Treatment), size = 3) + # sets colors by group and a set point size 
    xlab("") + 
    ylab("Bacterial Richness (SVs)") + 
    theme(text = element_text(size = 20)) # increases font size


# make that same graph but drop any samples that lack data for that FactorA (replace FactorA in the code with your factor name).
ggplot(data=subset(EX_ps_clean.rar.rich.df, !is.na(FactorA)), aes(x=Temperature, y=Observed)) + 
  theme_minimal() + 
  geom_point(aes(color=Group), size = 3) + # sets colors by group and a set point size 
  xlab("Temperature of Ocean Water (C)") + 
  ylab("Bacterial Richness (SVs)") + 
  theme(text = element_text(size = 20)) # increases font size



# Make an alpha diversity table: https://www.geeksforgeeks.org/create-table-from-dataframe-in-r/
EX_table = as.table(table(EX_ps_clean.rar.rich.df$Diet, EX_ps_clean.rar.rich.df$Week)) 
EX_table


 -----------
# Write out your research questions, make them specific. 
# Write a description of your experimental design, including notes about replication, repeated measures, nested design, or other factors which might impact your statistical models.




  # Alpha diversity metrics statistics (Lab 8)--------------
# Phyloseq can measure and visualize alpha diversity: https://joey711.github.io/phyloseq/plot_richness-examples.html
# phyloseq doesn't do stats or more complex graphs

#library(phyloseq)

# use phyloseq to measure alpha diversity
EX_ps_clean.rar.rich <- estimate_richness(EX_ps_clean.rar, measure=c("Observed", "Shannon")) #change to whatever measures you want

# use phyloseq to calculate Faith's Diversity metric (optional), https://rdrr.io/github/twbattaglia/btools/man/estimate_pd.html
EX_faith <- estimate_pd(EX_ps_clean.rar)

# OPTIONAL measure evenness for each SV individually
library(asbio)
library(microbiome)
EX_ps_clean.rar.even_SV <- evenness(otu_table(EX_ps_clean.rar))

# measure evenness for each sample
EX_ps_clean.rar.even <- EX_ps_clean.rar.rich$Shannon/log(EX_ps_clean.rar.rich$Observed)


# Coerce to data.frame and add the metadata for these samples
EX_ps_clean.rar.sd = as(sample_data(EX_ps_clean.rar), "matrix")
EX_ps_clean.rar.sd = as.data.frame(EX_ps_clean.rar.sd)
EX_ps_clean.rar.rich.df <- cbind(EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# make a histogram to look at the shape of the data (bell curve? skew?). You can save this graph for your own benefit if you want.
hist(EX_ps_clean.rar.rich.df$Observed)

# Want to know how much skew there is? Measure Kurtosis (a.k.a. tailedness). You can save this graph for your own benefit if you want.
install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
kurtosis(EX_ps_clean.rar.rich.df$Observed)
# example value output is 0.0006177001, looking for whether +/- and if a large/small number



head(EX_ps_clean.rar.rich.df)

#check the distribution of your data, which will change your approach
shapiro.test(EX_ps_clean.rar.rich.df$Shannon)
# W = 0.67622, p-value = 6.987e-06: my shannon diversity metric is not normally distributed (p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution)

shapiro.test(EX_ps_clean.rar.rich.df$Observed)
#W = 0.89057, p-value = 0.004976, my observed data are not normal. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution

shapiro.test(EX_ps_clean.rar.rich.df$EX_ps_clean.rar.even)
# W = 0.58972, p-value = 6.973e-07, not normally distributed. : p-value > 0.05 reveals that the distribution of the data are not significantly different from normal distribution



  ## If your alpha diversity is normally distributed (parametric)-----------------
library(vegan)

# pick your anova model based on your study design: https://www.statmethods.net/stats/anova.html. Replace "y" with the alpha diversity measure you want, such as observed, evvenness, shannon, etc.
# Replace "A and B" with the factors you want from your own data, and and replace x with numerical data

    ### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html

library(lme4)
library(lmerTest)
library(emmeans)

# you might run each factor separately to find out if they are significant on their own. 
summary(lm(Observed ~ as.numeric(PHOS), data= EX_ps_clean.rar.rich.df)) # PHOS is the phosphorus values from the lamb rumen data

# if you have multiple groups, this will give you the pairwise comparisons
emmeans(lm,pairwise ~ Diet) 
# paste the output here

# what if you want to add log transformation to a numerical variable in your metadata? Check out the syntax here: https://rpubs.com/marvinlemos/log-transformation


    ### Wilcoxon (a.k.a. Mann-Whitney) for pairwise numeric or categorical factors -------
# Wilcoxon in R (a.k.a. Mann-Whitney): https://data.library.virginia.edu/the-wilcoxon-rank-sum-test/
# Two data samples are independent if they come from distinct populations and the samples do not affect each other. Using the Mann-Whitney-Wilcoxon Test, we can decide whether the population distributions are identical while assuming them to follow the normal distribution. We can use this on non-normal data if we have many more than 30 samples per group and large variance between groups. Can only do pairwise comparisons (your factor can only have two variables in it)

# Replace "y" with the alpha diversity measure you want, "A" with the factor you want (can only do one as a time),
wilcox.test(y ~ A, data=EX_ps_clean.rar.rich.df)







    ### ANOVA for multiple categorical factors -------
# ANOVA is used when your factor has 2+ variables in it. If you have 3+ variables, you should run a post-hoc test to correct the p-value for multiple group comparisons

# One Way Anova (Completely Randomized Design)
fit <- aov(y ~ A, data=EX_ps_clean.rar.rich.df)

# Randomized Block Design (B is the blocking factor)
fit <- aov(y ~ A + B, data=EX_ps_clean.rar.rich.df)

# Two Way Factorial Design
fit <- aov(y ~ A*B, data=EX_ps_clean.rar.rich.df) 

# Analysis of Covariance
fit <- aov(y ~ A + x, data=EX_ps_clean.rar.rich.df) # x is a series of numerical values, like pH)


#example of one way anova from class
example1 <- lm(Observed ~ Diet, data=EX_ps_clean.rar.rich.df)
example2 <- lm(Observed ~ Week, data=EX_ps_clean.rar.rich.df)

#which is model better? Find the AIC value
anova(example1,example2)




#example of two way from class
summary(aov(Observed ~ Diet * Week, data=EX_ps_clean.rar.rich.df))

#For within subjects designs, the data frame has to be rearranged so that each measurement on a subject is a separate observation. See R and Analysis of Variance.

# One Within Factor
fit <- aov(y~A+Error(Subject/A),data=EX_ps_clean.rar.rich.df) #Subject would be a Factor in your metadata, y is the value of interest, A is the Factor of interest

# Two Within Factors W1 W2, Two Between Factors B1 B2
fit <- aov(y~(W1*W2*B1*B2)+Error(Subject/(W1*W2))+(B1*B2),
           data=EX_ps_clean.rar.rich.df)


# once you have chosen a model, add a correction if you have multiple group comparisons (multiple factor levels)
# Tukey Honestly Significant Differences
TukeyHSD(fit) # where fit comes from aov()




  ## If your alpha diversity is not normally distributed (non-parametric) -----------------
    ### Kruskal-Wallis Test -----
# K-W is the non-parametric version of ANOVA: http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
kruskal.test(Observed ~ Diet, data = EX_ps_clean.rar.rich.df)
# data:  Observed by Diet
# Kruskal-Wallis chi-squared = 2.2931, df = 1, p-value = 0.13


# Follow it up with a Conover Test if you have multiple comparisons. Note, code changed recently
# https://rdrr.io/cran/DescTools/man/ConoverTest.html 
install.packages("conover.test")
library(conover.test)

conover.test(EX_ps_clean.rar.rich.df$Observed, EX_ps_clean.rar.rich.df$Week,
            method = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) # method of correction, OK to pick just one
          
# output will look like this:
# Kruskal-Wallis rank sum test
# data: x and group
# Kruskal-Wallis chi-squared = 5.2867, df = 2, p-value = 0.07
# 
# Comparison of x by group                            
# (No adjustment)                                
# Col Mean-|
# Row Mean |          0          1
 ---------+----------------------
#       1 |  -1.818304
#         |     0.0417
#         |
#       2 |  -2.467699  -0.726045
#         |    0.0111*     0.2379
# 
# alpha = 0.05
# Reject Ho if p <= alpha/2

    ### Linear model for numeric factors -----
# Interpret linear model in R https://feliperego.github.io/blog/2015/10/23/Interpreting-Model-Output-In-R 
# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova
# https://boostedml.com/2019/06/linear-regression-in-r-interpreting-summarylm.html 



library(lme4)
library(lmerTest)
library(emmeans)

# you might run each factor separately to find out if they are significant on their own. 
summary(lm(Observed ~ as.numeric(PHOS), data= EX_ps_clean.rar.rich.df)) #CHANGE ME to update factor names or the alpha diversity metric you want


# you might run each factor separately on just a subset of the data.
summary(lm(Observed ~ as.numeric(PHOS), data=subset(EX_ps_clean.rar.rich.df, Week=="2"))) #CHANGE ME to update factor names or the alpha diversity metric you want

# what if you want to add log transformation to a numerical variable in your metadata? Check out the syntax here: https://rpubs.com/marvinlemos/log-transformation


    ### Linear mixed effects models for complicated experimental designs -----
# If you have more complicated experimental designs, you might need lmer or glmer models to accommodate that
# does not have a normal distribution, so compare means using glmer and poisson distribution (number of events within a time interval): https://towardsdatascience.com/the-poisson-distribution-and-poisson-process-explained-4e2cb17d459

lm <- (glmer(Observed ~ A + C + (1|Year_collected), family=poisson, data=EX_ps_clean.rar.rich.df))

summary(lm)
# paste the output here


# if you have multiple groups, this will give you the pairwise comparisons
emmeans(lm,pairwise ~ A) 
# paste the output here

emmeans(lm,pairwise ~ B) 
# paste the output here





# More ways to graph alpha diversity (Lab 9) -------
# Resources for visualization
# https://www.data-to-viz.com/
# https://r-graph-gallery.com/index.html

  ## Heatmaps --------------
# base plot
plot_heatmap(EX_ps_clean.rar)

# follow the tutorial to make prettier versions in phyloseq: https://joey711.github.io/phyloseq/plot_heatmap-examples.html

# example from class
EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")

plot_heatmap(EX_ps_clean.rar, fill="Phylum") +   facet_grid(~Week, space="free", scales="free") + theme(legend.position = "bottom") #, axis.text.x = element_blank()) 

#example from class, but looks like trash without a lot of fancying up
heatmap(otu_table(EX_ps_clean.rar))


  ## Correlogram --------------
# This makes a correlation matrix plot: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

library(RColorBrewer)
library(corrplot)


# corrplot requires one dataframe of metadata and/or seqtab data.  Any columns in the dataframe will be correlated against all others.  Too many columns makes the graph unreadable, try to keep it to <50.

# select to top 15 most abundant SVs from your phyloseq object and extract to a dataframe
#take out top 100 taxa based on abundance
topNOTUs = names(sort(taxa_sums(EX_ps_clean.rar), TRUE)[1:15])
top15 = prune_taxa(topNOTUs, EX_ps_clean.rar)

# combine with your metadata and create one dataframe object. you can include other info that you created for a previous dataframe, as long as those objects are still in your R environment. Reminder, you can only use numeric data in a correlation matrix, so you will have to drop certain columns or make them numeric instead.
# Coerce to data.frame and add the metadata for these samples
top15.sd = as(otu_table(top15), "matrix")
top15.sd = as.data.frame(top15.sd)

# add Genus names in place of the full sequence name that is in the SV columns
top15.tax <- as.data.frame(tax_table(top15))
## if the Genus is empty, replace with the Family
top15.tax$Genus[is.na(top15.tax$Genus)] <- top15.tax$Family[is.na(top15.tax$Genus)]
colnames(top15.sd) = top15.tax$Genus

# paste all the components together
EX_ps_clean.rar.corr.df <- cbind(top15.sd, EX_ps_clean.rar.rich, EX_ps_clean.rar.even, EX_ps_clean.rar.sd)

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# change any column factor names to make them prettier
names(EX_ps_clean.rar.corr.df)[names(EX_ps_clean.rar.corr.df) == "EX_ps_clean.rar.even"] <- "Evenness"

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# for example, drop a column with factorial data
EX_ps_clean.rar.corr.df <- subset(EX_ps_clean.rar.corr.df, select = -c(Treatment, Diet, Sheep_ID, ICTERIC_INDEX, LIPEMIC_INDEX))

# check header to make sure it looks good
head(EX_ps_clean.rar.corr.df)

# check that all remaining columns are numeric instead of factor or character
str(EX_ps_clean.rar.corr.df)

# clean up any columns which are not registering as numeric
EX_ps_clean.rar.corr.df <- sapply(EX_ps_clean.rar.corr.df, as.numeric)



# alternatively, load a premade dataframe containing your chosen SVs from your sequence table (otu)table in phyloseq) and the metadata you want to include
# EX_ps_clean.rar.corr.df <- read.csv("example_correlogram_dataframe.csv", check.names = FALSE, header=T, row.names=1)



# run correlations
ex_corr_calc <- cor(EX_ps_clean.rar.corr.df, 
                    use="complete.obs", # use=Specifies the handling of missing data. Options are all.obs (assumes no missing data - missing data will produce an error), complete.obs (listwise deletion), and pairwise.complete.obs (pairwise deletion)
                    method="spearman") # correlation method=pearson, spearman, or kendall

### Note, if you have too few samples, you may receive an error about too few observations. You may ignore the error message, or remove than column



### test significance of correlations, copy and paste this whole chunk together
cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(EX_ps_clean.rar.corr.df,0.95)
res2 <- cor.mtest(EX_ps_clean.rar.corr.df,0.99)
# Note: if you come up with a warning message, it means that one or more of your columns generated correlations of value 0, which makes it angry.  Visualize the plot, and "?" will come up in those columns.  To fix, remove the column or add a data transformation, to the dataframe you created to make this.  Set the corr.mtest function and run corr.mtest again.


## plot only those correlations with a significant p-value <0.05, using hierarchical clustering
corrplot(ex_corr_calc, 
         type="lower", #shape of the plot itself: full, upper, lower
         method="circle", #shape of the data: circle, square, ellipse, number, shade, color, pie 
         order="hclust", #how to cluster samples: AOE, hclust, FPC, alphabet, or leave blank
         p.mat = res1[[1]], #which correlations to use
         sig.level=0.05, #sets significance cutoff
         insig="blank", #leaves p > 0.05 blank
         tl.col="black", # text color
         tl.cex=.9, #text size
         col=brewer.pal(n=10, name="RdYlBu")) #specify a color palette and number of shades needed





  ## Barplots------------------
# can add ggplot components to make it pretty
plot_bar(EX_ps_clean.rar, fill="Phylum") #don't recommend using genus here, it make crash R


# agglomerate SVs by taxa level to get rid of all the blacklines on the graph. can be done for any level of taxonomy
EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Genus")

# WARNING: if you group along the x-axis phyloseq just stacks seqence counts, it DOES NOT calculate show average for that group
plot_bar(EX_ps_clean.rar.glom, fill="Genus")

# example from class
EX_ps_clean.rar.glom = tax_glom(EX_ps_clean.rar, "Phylum")
plot_bar(EX_ps_clean.rar.glom, fill="Phylum") +   facet_grid(~Diet, space="free", scales="free") + theme(legend.position = "bottom", axis.text.x = element_blank()) 


# make a stacked 100% bar chart in phyloseq

EX_ps_clean.rar.stacked = transform_sample_counts(EX_ps_clean.rar, function(x) x / sum(x) )

plot_bar(EX_ps_clean.rar.stacked, fill="Phylum") 

# to filter by abundance and pool low abundance groups: https://github.com/joey711/phyloseq/issues/901


  ## Core community members -----------
#code from this website: https://microbiome.github.io/tutorials/CoremicrobiotaAmplicon.html

#load libraries as needed
library(dplyr)
library(microbiome)


# grab the taxonomy to clean up first
tax.df <- data.frame(tax_table(EX_ps_clean.rar), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# put the new taxa back in
tax_table(EX_ps_clean.rar) <- as.matrix(tax.df)

# check it worked
as.data.frame(tax_table(EX_ps_clean.rar))$Genus


#calculate compositional version of the data (relative abundances)
core_biomeW.rel <- microbiome::transform(EX_ps_clean.rar, "compositional") #CHANGE ME to your phyloseq object

#This returns the taxa that exceed the given prevalence and minimum abundance detection thresholds. Set your preferred thresholds.
core_taxa_standardW <- core_members(core_biomeW.rel, detection = 1/10000, prevalence = 60/100) #CHANGE ME

#A full phyloseq object of the core microbiota at those limits is obtained as follows
phylo.coreW <- core(core_biomeW.rel, detection = 1/10000, prevalence = .6)

###retrieving the associated taxa names from the phyloseq object and add it to what you just made.
core.taxaW <- taxa(phylo.coreW)
class(core.taxaW)

# get the taxonomy data and assign it to what you just made.
tax.matW <- tax_table(phylo.coreW)
tax.dfW <- as.data.frame(tax.matW)

# add the SVs to last column of what you just made.
tax.dfW$SV <- rownames(tax.dfW)

## CHANGE ME, write down how many taxa are shared between these samples.

# select taxonomy of only those OTUs that are core members based on the thresholds that were used.
core.taxa.classW <- dplyr::filter(tax.dfW, rownames(tax.dfW) %in% core.taxaW)
knitr::kable(head(core.taxa.classW))

# save the list so you can access later, can report just the list
write.csv(core.taxa.classW, "~/core.taxa.example.csv")


# graph the abundance of those shared taxa, here are some example: https://microbiome.github.io/tutorials/Core.html

#aggregate the genera so we don't get a lot of lines separating all the SVs
plot.gen <- aggregate_taxa(phylo.coreW, "Genus")

# load libraries as needed
library(ggplot2)
library(RColorBrewer)
library(viridis)

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)

plot_core(plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = .5) + #CHANGE min prevalence
  xlab("Detection Threshold (Relative Abundance (%))") + ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis()



# Comparing changes in taxonomy (Lab 10) ---------------------------------
  ## Focusing on a single taxon -------
library(phyloseq) 
library(dplyr) 
library(tidyr)  
library(ggplot2)

# transform to relative abundance
relabun.ps <- transform_sample_counts(EX_ps_clean.rar,function(x) x / sum(x)) 

# glom ASVs by genus
ps_genus <- tax_glom(relabun.ps, taxrank = "Genus", NArm = FALSE) 

#subset by the genus of choice. Taxon is not present if you get error: "length of 'dimnames' [1] not equal to array extent" 
ps_genusP <- subset_taxa(ps_genus, Genus %in% c("Streptococcus"))                      

# melt the data into a different configuration
genus.df <- psmelt(ps_genusP) 

# grab that abundance data. CHANGE ME where it says Factor1, Factor2
MySummary <- genus.df %>% group_by(Diet, Week, Genus) %>% summarize(mean_abund = mean(Abundance, na.rm=TRUE)) 

#check it
head(MySummary)

# graph it
ggplot(data=MySummary, aes(x = Week, y = mean_abund)) +
  geom_point(aes(color=Diet)) + ylab("Mean relative abundance of reads")



  ## Differential Abundance with DESeq ---------------------------------
# DESeq only does pairwise comparisons. To make a multifactorial comparison and graph, use "DESeq_and_ternary_plot_example.R".  You can also subset your data
library(DESeq2)
packageVersion("DESeq2")

# OPTIONAL if you need to subset 
EX_ps_clean_subsetA = subset_samples(EX_ps_clean, FACTOR == "A") #CHANGE ME 
EX_ps_clean_subsetA <- prune_samples(sample_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)
EX_ps_clean_subsetA <- prune_taxa(taxa_sums(EX_ps_clean_subsetA) > 0, EX_ps_clean_subsetA)


# grab phyloseq data for use in deseq
diagdds = phyloseq_to_deseq2(EX_ps_clean, ~ Diet) #CHANGE ME to an factor with only 2 levels

# calculate differential abundance
gm_mean =  function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

# calculate significance for those abundance calculations
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(EX_ps_clean)[rownames(sigtab), ], "matrix")) #CHANGE ME if you didn't subset your data

head(sigtab)

dim(sigtab)
#CHANGE ME number of SVs that were different between them


# calculate log changes and set
 sigtab = sigtab[sigtab[, "log2FoldChange"] < 0, ] #CHANGE ME. use this to select only positive (or negative) log changes
sigtab = sigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")] #CHANGE ME add Order or Species - if you have it

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))


## if the Genus is empty, replace with the Family NOTE: this used to work but the syntax is broken for unknown reasons
# sigtab.df$Genus[is.na(sigtab.df$Genus)] <- sigtab.df$Family[is.na(sigtab.df$Genus)]

# if the Genus is empty, replace with the Family NOTE: works as of Feb 2020
sigtab$Genus = ifelse(is.na(sigtab$Genus), paste(sigtab$Family),paste(sigtab$Genus));sigtab

# OPTIONAL bind Genus and Species together - only works if you had species-level taxonomy AND you added it a few steps prior
# sigtab$Genus.species <- paste(sigtab$Genus, sigtab$Species)





library("ggplot2")

## graph differential abundance
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(aes(size=baseMean)) + #scale size by mean relative abundance
  theme(axis.text.x = element_text(hjust = 0, vjust=0.5, size=10), axis.text.y = element_text(size=10)) + xlab("log2 Fold Change") + labs(size = "Mean Sequence Abundance") + theme_minimal()



  ## Feature Prediction (Differential Abundance) with Random Forest ---------------------------------
#  https://rpubs.com/michberr/randomforestmicrobe
# https://www.stat.berkeley.edu/~breiman/RandomForests/cc_home.htm
# https://www.rdocumentation.org/packages/randomForest/versions/4.6-12/topics/randomForest

library(phyloseq)
library(vegan)
library(plyr)
library(dplyr)
require(magrittr)
require(scales)
require(grid)
require(reshape2)
require(knitr)
library(randomForest)
library(rfPermute)
library(ggplot2)
library(RColorBrewer)

# Make a dataframe of training data with OTUs as column and samples as rows, which is the phyloseq OTU table
predictors <- otu_table(EX_ps_clean)

dim(predictors)
#CHANGE ME how many samples and SVs includes

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)


# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

# clean up some of the other taxa info
colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))

### start here when choosing factors, can reuse the above lines as needed. one example for factorial data, and one for numeric data is provided. select as needed.

# Make one column for our outcome/response variable. Choose which one applies to the thing you want to test, and then follow factorial or numeric through the rest of the code section.
response <- as.factor(sample_data(EX_ps_clean)$Factorial_Data) #CHANGE ME
response <- as.numeric(sample_data(EX_ps_clean)$Numeric_Data) #CHANGE ME

# Combine response and SVs into data frame
rf.data <- data.frame(response, predictors)


# set seed for random number generation reproducibility
set.seed(2)

# classify for factorial data
response.pf <- rfPermute(response ~. , data = rf.data, na.action = na.omit, ntree= 500, nrep = 100) #na.omit ignores NAs in data (not tolerated). ntrees is how many forests to build, nreps generates p-value

# or this way for numeric data
response.pf <- rfPermute(as.numeric(response) ~. , data = rf.data, na.action = na.omit, ntree= 500, nrep = 100)

print(response.pf)
# paste the print out here, especially the OOB error. 1-(Out-of-the-box error) = accuracy of your model


# grab which features were labeled "important"
imp <- importance(response.pf, scale = TRUE)

# Make a data frame with predictor names and their importance
imp.df <- data.frame(predictors = rownames(imp), imp) 


# For factorial data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, MeanDecreaseAccuracy.pval <= 0.05) 
print(dim(imp.sig))
# CHANGE ME, how many SVs (sig features) were identified

# For numeric data, grab only those features with p-value < 0.05
imp.sig <- subset(imp.df, IncNodePurity.pval <= 0.05)
print(dim(imp.sig))
# CHANGE ME, how many SVs (sig features) were identified


# or For factorial data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$MeanDecreaseAccuracy),]


# For numeric data, sort by importance amount
imp.sort <- imp.sig[order(imp.sig$IncNodePurity),]



#create levels to the factor based on SV table
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

# Select the top so many predictors (more than 50 is a crowded graph)
imp.top <- imp.sort[1:50, ]

# figure out what they are and name them, otherwise will just be named the full SV
otunames <- imp.top$predictors


# need to grab the abundance data out of the otu_table and make a new data.frame, and add the taxa names back in
# grab the column names from the otu table that match those in the forest set
pred.abun.colnum <- which(colnames(rf.data) %in% otunames)

# when you find a match, grad the abudnance data
pred.abun <- rf.data[,sort(c(pred.abun.colnum))]

# make this into a dataframe for manipulation
pred.abun.df <- data.frame(pred.abun, stringsAsFactors=FALSE)

# use the row.names (sample names) from the phyloseq object to name the samples in your forest
row.names(pred.abun.df) <- row.names(sample_data(EX_ps_clean))


# set color palette  
col.bro <- (rainbow(6))

# add white to that color list
col.bro <- append(col.bro, "#ffffff", after = 6)

# add some factors that you can use to make your graph pretty, as many as you want
pred.abun.df$Sample <- row.names(sample_data(EX_ps_clean)) #always grab the sample names
pred.abun.df$FactorA <- sample_data(EX_ps_clean)$FactorA #CHANGE ME
pred.abun.df$FactorB <- sample_data(EX_ps_clean)$FactorB #CHANGE ME


# optional, if you want factors to graph in a specific order, you can set that manually, and relabel them so they are more readable in the graph
pred.abun.df$Date <- factor(pred.abun.df$Date, levels=c("21_Apr", "12_May", "01_Jun", "22_Jun", "25_Jul"), labels=c("21 Apr", "12 May", "01 Jun", "22 Jun", "25 Jul")) #CHANGE ME

# reload these packages in this order, because sometimes the ddply function breaks
 library(plyr)
 library(dplyr)

# melt and transform the data using ALL the factors you added
m <- melt(pred.abun.df, id.vars=c("Sample", "FactorA", "FactorB"))
m <- ddply(m, .(variable), transform, rescale = log(1 + value)) #NOTE: may need to change it to as.numeric(value) if it tells you an error about non-binary operators. if you go to graph it and it doesn't recognize 'rescale', it means this piece failed.


# adjusted plot for factorial data, recommend using sample ID as 'factor A'
ggplot(m, aes(as.factor(FactorA), variable)) + 
  theme_minimal() + 
  facet_grid(.~FactorA, space="free", scale="free") + #set up graph facets to separate out levels of FactorA
  geom_tile(aes(fill = rescale), color="gray") + #add the heatmap coloring 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + #use the preset color palette
  labs(fill="Log abundance") + #rename legend heading
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 0, size = 8),
        axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 6)) +
  ylab('Predictor Taxa') +
  xlab('Factor A') 


# adjusted plot for numeric data, recommend using sample ID as 'factor A'
for_moist <-ggplot(m, aes(FactorA, variable)) + #can use as.factor(FactorA) if don't have consistent values across gradient and will cut out white space on x-axis
  theme_minimal() + 
  geom_tile(aes(fill = rescale), color="gray") + 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + 
  labs(fill="Log abundance") +
  theme(legend.position = 'bottom', # could be none, top, right, left
        axis.text.x = element_text(angle = 0, size = 10),
        #axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 10))+
 # scale_x_discrete(breaks=c(0,24.4, 50.3, 74.8, 98.3), labels=c("0" = "0", "24.4" ="25", "50.3" = "50", "74.8" = "75", "98.3" ="100")) + #optional, may want to add x-axis scale breaks
  ylab('Predictor Taxa') +
  xlab('Factor A') 



  ## LEFSe ---------------------------------
# can be done within R using the lefse package

# can transform phyloseq data to lefse-ready format with the command that you create, phyloseq_to_lefse
# uses code by seahorse001x: https://github.com/seashore001x/Rrumen/blob/master/phyloseq2lefse.R

# require(dplyr)
# require(tibble)
# 
# # this script defines a function to convert phyloseq object into lefse recognized file format. no need to change anything in this bit.
# phyloseq_to_lefse <- function(physeq){
#   # aggregate to genus level
#   ps_lefse <- physeq %>% tax_glom(taxrank = 'Genus', NArm = F)
#   
#   # extract taxonomic data from phyloseq object and then stored in a vector called lefse_tax
#   lefse_tax <- ps_lefse %>% tax_table %>% data.frame(stringsAsFactors=FALSE)
#   lefse_tax <- replace(lefse_tax, is.na(lefse_tax), 'Unclassified')
#   lefse_tax <- lefse_tax %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% mutate(id = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = "|")) %>% ungroup %>% pull(id)
#   
#   # extract otu abundance matrix from phyloseq object and annotated with tax information
#   lefse_matrix <- otu_table(ps_lefse) %>% data.frame(stringsAsFactors = F) %>% t %>% data.frame
#    
#   # bookmark this is what is throwing the error
# #  colnames(lefse_matrix) <- lefse_tax
#   
# #  row.names(lefse_matrix) <- lefse_tax
#   
#   # extract sample metadata and order sample same in lefse_matrix
#   lefse_sample <- sample_data(ps_lefse)
#   
#   
#   # convert factor in lefse_sample into character in order to combine sample and abundance data
#   lefse_sample_isfactor <- sapply(lefse_sample, is.factor)
#   lefse_sample[,lefse_sample_isfactor] <- lefse_sample[,lefse_sample_isfactor] %>% lapply(as.character)
#   lefse_sample <- lefse_sample %>% data.frame
#   
#   lefse_table <- full_join(rownames_to_column(lefse_sample), rownames_to_column(lefse_matrix), by = ("rowname" = "rowname")) %>% t
#   
#   return(lefse_table)
# }
# 
# 
# EX_clean_for_lefse <- phyloseq_to_lefse(EX_ps_clean)




predictors <- otu_table(EX_ps_clean)

dim(predictors)
#CHANGE ME how many samples and SVs includes

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)


# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(predictors) <- tax.df$Genus.species

# clean up some of the other taxa info
colnames(predictors) = gsub("_unclassified", "", colnames(predictors))
colnames(predictors) = gsub("_Intercertae_Sedis", "", colnames(predictors))

### start here when choosing factors, can reuse the above lines as needed. one example for factorial data, and one for numeric data is provided. select as needed.

# Make one column for our outcome/response variable. Choose which one applies to the thing you want to test, and then follow factorial or numeric through the rest of the code section.
response.class <- as.factor(sample_data(EX_ps_clean)$Diet) #CHANGE ME to the main factor you want to test
response.subclass <- as.factor(sample_data(EX_ps_clean)$Week)

subjects <- as.factor(sample_data(EX_ps_clean)$Sheep_ID) # change this to be ID of the individual if animals/environments have repeated measures sampling

# Combine response and SVs into data frame
Ex_ps_for_lefse.df <- data.frame(response.class, response.subclass, subjects, predictors)


# save this output and upload it to the LEFSe web version on Galaxy: http://huttenhower.sph.harvard.edu/galaxy/
write.table(Ex_ps_for_lefse.df, file = "Ex_ps_for_lefse.df.txt", 
            append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)










#  Beta diversity ordinations and stats (Lab 11 and 12) ---------------------------------
# ordinations in phyloseq: https://joey711.github.io/phyloseq/plot_ordination-examples.html

  ## PCA  -----------------
# currently, phyloseq doesn't run a PCA, but you can manually perform one using this tutorial: https://www.datacamp.com/community/tutorials/pca-analysis-r

install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot") #takes a few minutes
library(ggbiplot)

# calculate the components
EX_clean.rar.pca <- prcomp(otu_table(EX_ps_clean.rar), center = TRUE, scale = TRUE) #CHANGE ME

# take a look
summary(EX_clean.rar.pca) #CHANGE ME

# graph it. add ggplot2 text to make it pretty
ggbiplot(EX_clean.rar.pca)




  ## PCoA in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_pcoa <- ordinate(EX_ps_clean.rar, #calculate similarities
                   method ="PCoA", #ordination type
                   "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples", color="Diet") #CHANGE ME but leave type as "samples"

# fancy ordination
library(ggplot2) # can add components on top of the phyloseq graph
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_pcoa, type="samples") + #CHANGE ME
  geom_point(aes(size = as.factor(Week), color=as.factor(Diet), shape=as.factor(Week), alpha = as.numeric(Weight_kg))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  theme_minimal()+
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=as.factor(Week), linetype=as.factor(Week))) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Diet", size = "Week", shape = "Week", alpha = "Weight (kg)", linetype="Week") #rename the headers in the legend


  ## NMDS in phyloseq -----------------

# use phyloseq to calculate ordination for use in the plot
EX_uJ_nmds <- ordinate(EX_ps_clean.rar, #calculate similarities
                    method ="NMDS", #ordination type
                    "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted) but is usually run as binary=FALSE.

# report the amount of stress at the end of the NMDS calculation, as it is typical to report that in your manuscript. >0.2 is not great, and >0.3 plot is meaningless
# Paste output here: Run 20 stress 9.04618e-05 #CHANGE ME

# simple ordination
plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", color="Diet") #CHANGE ME but leave type as "samples"


# fancy ordination
library(ggplot2)
install.packages("viridis")
library(viridis) #add a cool color palette

plot_ordination(EX_ps_clean.rar, EX_uJ_nmds, type="samples", shape="Diet") + #CHANGE ME
  geom_point(aes(size = as.numeric(Weight_kg), color=as.factor(Week), alpha = as.factor(Week))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group=Diet, linetype=Diet)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Weeks 0 to 2", size = "Weight in kg", shape = "Diet", linetype="Diet", alpha= "Week") #rename the headers in the legend




### ----- permanova stats for ordinations --------

# Example basic permanova test using Adonis in vegan
adonis2(EX_uJ_nmds ~ Diet * Week, as(sample_data(EX_ps_clean.rar), "data.frame"), permutations=9999, na.action=na.omit) #CHANGE ME so variable and factor reflect your data
# paste output here

# Example permanova output
#             Df  SumsOfSqs MeanSqs F.Model R2      Pr(>F)   
#   Diet       1    0.6479 0.64792  1.7560 0.07303 0.0047 **
#   Week       1    0.7136 0.71361  1.9340 0.08044 0.0019 **
#   Diet:Week  1    0.4997 0.49967  1.3542 0.05632 0.0495 * 
#   Residuals 19    7.0106 0.36898         0.79021          
#   Total     22    8.8718                 1.00000          
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# interpretation: The F-model/F statistic is a measure of importance, R2 may not be relevant for factorial data



# Or, run a more complicated permanova test using Adonis in vegan if you have a complicated experimental design

# if want to add repeated sampling use this line
ctrl <- with(as.data.frame(sample_data(EX_ps_clean.rar)), how(blocks = Sheep_ID, nperm = 999))

adonis2(EX_uJ_nmds ~ Diet * Week * as.numeric(Weight_kg), data = data.frame(sample_data(EX_ps_clean.rar)), # CHANGE ME so 'Factor' is the name of your Factor variables
        permutations = ctrl, #for repeated measures
       # strata = as.data.frame(sample_data(EX_ps_clean.rar))$Blocking_factor #CHANGE ME for replicates/blocks
       ) 


### NOTE: if adonis is not working and says your data is of the wrong type:
# 1) check to make sure you don't have NAs in what you are trying use
# 2) update the data.table package, then reload the vegan package
# 3) Repeat the ordination calculation by trying the block of code below:

# Code to use Vegan to calculation ordination data and run stats, in case phyloseq and vegan aren't sharing code well.
# function and tutorial from https://rpubs.com/DKCH2020/587758
library(vegan)

veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}
# export data from phyloseq to vegan-compatible object
EX_vegan <- veganotu(EX_ps_clean.rar)

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
sampledf <- data.frame(sample_data(EX_ps_clean.rar))

# run the ordination calculation, change the variable name from EX_uJ to reflect the calculation method you choose (EX_wBC)
EX_uJ <- vegdist(wisconsin(sqrt(EX_vegan)), method = "jaccard", binary=TRUE) #CHANGE ME to be the method and weighted (binary = false) or unweighted version (binary = true)

adonis2(EX_uJ ~ Diet * Week, as(sample_data(EX_ps_clean.rar), "data.frame"), permutations=9999, na.action=na.omit)


#### run betadispersion test to see how tight/loose the clusters are ----
betadisp_EX <- betadisper(EX_uJ_nmds, sampledf$Diet) # CHANGE me so variable and factor reflect your data

betadisp_EX
# paste output here

# Example Call: betadisper(d = EX_BC, group = sampledf$Diet)

# No. of Positive Eigenvalues: 23
# No. of Negative Eigenvalues: 0
# Average distance to median:
# LooseAlfalfa PelletedAlfalfa 
# 0.6472          0.5767 
# Eigenvalues for PCoA axes: (Showing 8 of 23 eigenvalues)
# PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
# 1.3428 0.7263 0.6223 0.5387 0.5006 0.4863 0.4460 0.4361 


# run ANOVA to see if clusters overlap or not
anova(betadisp_EX)
# paste output here

# Example output Analysis of Variance Table
# Response: Distances
# Df   Sum Sq   Mean Sq F value    Pr(>F)    
# Groups      2 0.043240 0.0216200  94.901 < 2.2e-16 ***
#   Residuals 109 0.024832 0.0002278                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


permutest(betadisp_EX)
# paste output here

# Example output Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
# Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups      2 0.043240 0.0216200 94.901    999  0.001 ***
#   Residuals 109 0.024832 0.0002278                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


#  correct for multiple comparisons
TukeyHSD(betadisp_EX)
# paste output here

# Example output Tukey multiple comparisons of means
# 95% family-wise confidence level
# Fit: aov(formula = distances ~ group, data = df)
# $group
# diff          lwr       upr     p adj
# wild veligers-veligers 0.099436111  0.079792204 0.1190800 0.0000000
# tank-veligers          0.106213720  0.087878344 0.1245491 0.0000000
# tank-wild veligers     0.006777609 -0.002106682 0.0156619 0.1702206




  ## 3D scatterplot --------

# Use the vegan3d package and any distance data to graph it in 3D
# https://rdrr.io/cran/vegan3d/man/ordiplot3d.html
install.packages('vegan3d')
library(vegan3d)
install.packages('scatterplot3d')
library(scatterplot3d)
# library(vegan)


# basic version of the plot, using the object made in the Betadispersion section above
ordiplot3d(EX_uJ)



# Beta diversity components (Lab 13) ---------------------------------

# Component analysis perform ordination calculations on your sample data (metadata or environmental data) and on your sequence table data (SVs). Having missing or NA in the datasets may throw an error. Having outliers may distort the plots. 

# optional, if you need to drop samples before creating the plots below, use this:
# EX_ps_clean.rar_altered = subset_samples(EX_ps_clean.rar, Sample_type == "Mock") # CHANGE ME to the sample to drop
# EX_ps_clean.rar_altered <- prune_taxa(taxa_sums(EX_ps_clean.rar_altered) > 0, EX_ps_clean.rar_altered)

# optional, use this if you want to convert N/A in your metadata to something else
# sample_data(EX_ps_clean.rar)$Factor[is.na(sample_data(EX_ps_clean.rar)$Factor)] <- "Something_else" # CHANGE ME so factor and replacement reflect the column you want to search through and what you want to replace NA with 


  ## CCA in phyloseq -------------------------------------------
# correspondence analysis, or optionally, constrained correspondence analysis (a.k.a. canonical correspondence analysis), usually used with bell-curve or unimodal relationships hypothesis-based testing does not need normally distributed data, but if you have a lot of outliers you may want to consider adding a data transformation.

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.

# reload as needed
library(phyloseq)
library(vegan)
library(ggplot2)

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray") #CHANGE ME

cca_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGE ME
  method = "CCA",
  distance = bray_not_na,
  formula = ~ Diet * Week +  # CHANGE ME add factors
    Condition(Sheep_ID) #CHANGE ME use for repeated measures
)

cca_ord
# paste the output here

# Call: cca(formula = OTU ~ Diet + Week, data = data)
# 
# Inertia Proportion Rank
# Total          7.2413     1.0000     
# Constrained    0.9964     0.1376    2
# Unconstrained  6.2450     0.8624   20
# Inertia is scaled Chi-square 
# 
# Eigenvalues for constrained axes:
#   CCA1   CCA2 
# 0.6864 0.3100 
# 
# Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
# 0.9639 0.6452 0.5080 0.4325 0.4097 0.3276 0.3207 0.3114 
# (Showing 8 of 20 unconstrained eigenvalues)


# NOTE: do you have CCA eigenvalues for each of the factors you put in? CCA1 = x axis, and CCA2 = y axis, so if you don't have both of those it will use CA which is the unconstrained axis instead.



# anova of whole model
anova(cca_ord, permu=1000)
# paste the output here    
# Permutation test for cca under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = OTU ~ Diet + Week, data = data)
#           Df  ChiSquare F       Pr(>F)    
# Model     2    0.9964 1.5955  0.001 ***
# Residual 20    6.2450                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# anova of the factors (terms) you specified
anova(cca_ord, by="terms", permu=1000)
# Permutation test for cca under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: cca(formula = OTU ~ Diet + Week, data = data)
#           Df ChiSquare    F   Pr(>F)    
# Diet      1    0.4152 1.3297  0.058 .  
# Week      1    0.5812 1.8613  0.001 ***
# Residual 20    6.2450                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1             


# cca plot with FACTORA in the model
cca_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = cca_ord, 
  color = "Diet", #CHANGE ME
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Week)) + #CHANGE ME
  geom_point(aes(colour = Diet), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cca_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CCA1, 
                 yend = CCA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CCA1, 
                 y = 1.3 * CCA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cca_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 





  ## RDA in phyloseq -------------------------------------------
# Redundancy analysis for linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

rda_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGEME
  method = "RDA",
  distance = bray_not_na,
  formula = ~ Diet * Week + # CHANGE ME add factors
  Condition(Sheep_ID) #CHANGE ME use for repeated measures
)

rda_ord
# paste the output here

# NOTE: do you have RDA eigenvalues for each of the factors you put in? RDA1 = x axis, and RDA2 = y axis, so if you don't have both of those it will use PC which is the unconstrained axis instead.



# anova of whole model
anova(rda_ord, permu=1000)
# paste the output here    

# anova of the factors (terms) you specified
anova(rda_ord, by="terms", permu=1000)
# paste the output here              


# rda plot with FACTORA in the model
rda_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = rda_ord, 
  color = "Diet", #CHANGE ME
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Week)) + #CHANGE ME
  geom_point(aes(colour = Diet), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 



# Now add the environmental variables as arrows
arrowmat <- vegan::scores(rda_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = RDA1, 
                 yend = RDA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * RDA1, 
                 #y = 1.3 * RDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 



  ## dbRDA in phyloseq -------------------------------------------
# Constrained Analysis of Principal Coordinates or distance-based RDA, for non-linear relationships

# first create a distance ordination. BE SURE THERE ARE NO NAs in your factorial data!! If you need to alter something, see the "beta diversity components" section.
bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

cap_ord <- ordinate(
  physeq = EX_ps_clean.rar, #CHANGEME
  method = "CAP", #borrows capscale from the vegan package
  distance = bray_not_na,
  formula = ~ Diet * Week + # CHANGE ME add factors
    Condition(Sheep_ID) #CHANGE ME use for repeated measures
)

cap_ord
# paste the output here

# NOTE: do you have CAP eigenvalues for each of the factors you put in? CAP1 = x axis, and CAP2 = y axis, so if you don't have both of those it will use MDS which is the unconstrained axis instead.

# anova of whole model
anova(cap_ord, permu=1000)
# paste the output here    

# anova of the factors (terms) you specified
anova(cap_ord, by="terms", permu=1000)
# paste the output here              


# CAP plot with FACTORA in the model
cap_plot <- plot_ordination(
  physeq = EX_ps_clean.rar, #CHANGE ME
  ordination = cap_ord, 
  axes = c(1,2)) + 
  theme_minimal() +
  aes(shape = as.factor(Week)) + #CHANGE ME
  geom_point(aes(colour = Diet), size = 4) + #CHANGE ME
  geom_point(colour = "grey90", size = 1.5) 




# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    linewidth = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  ) 



  ## variance partitioning --------------
# How much did any one factor contribute to the sample clustering?
# http://www.hiercourse.com/docs/variationPartitioning.html 

library(vegan)

# create a dataframe from your SV table
EX.df <- as.data.frame(otu_table(EX_ps_clean.rar))

# extract your sample data from the phyloseq object
env.df <- as.data.frame(sample_data(EX_ps_clean.rar))


# calculate Principal Coordinates Of Neighborhood Matrix, diversity distance data transformed into rectangular format
EX.pcnm <- pcnm(dist(bray_not_na))

# environmental variables as predictors of community similarity
cap.env <- capscale(EX.df ~ ., data=env.df, distance='bray')

# calculate CCA 
cap.pcnm <- capscale(EX.df ~ ., data=as.data.frame(scores(EX.pcnm)), distance='bray')


# make a model using SV ordination and sample data
mod0.env <- capscale(EX.df ~ 1, data=env.df, distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# make a model using SV ordination and the component scores
mod0.pcnm <- capscale(EX.df ~ 1, data=as.data.frame(scores(EX.pcnm)), distance='bray') # add + Condition(SAMPLE_ID) for repeated measures

# select variables in each predictor table
step.env <- ordistep(mod0.env, scope=formula(cap.env))

# check variance inflation factors, higher number = data are redundant to another factor, 1= data are unique.  If factors are redundant (conflated) to each other (for example they basically report the same thing, go back and remove one and mention they were redundant/conflated in your manuscript Methods)
vif.cca(step.env)


step.pcnm <- ordistep(mod0.pcnm, scope=formula(cap.pcnm))

step.env$anova
# paste output here

step.pcnm$anova
# paste output here

EX.pcnm.sub <- scores(EX.pcnm, 
                            choices=c(1,3:16)) #CHANGE ME to include the significant pcnm numbers from previous command. If you want them all, write 1:16. to select only a few of them, write something like 1,3,6,7:9. 

# partition variation among four predictor tables:
EX.var <- varpart(EX.df, 
                        ~ FACTORA, #CHANGE ME
                        ~ FactorB, #CHANGE ME
                        EX.pcnm.sub, data=env.df)

#plot 
par(mfrow=c(1, 2))
showvarparts(4)
plot(EX.var)



EX.var
# paste output here

anova(rda(EX.df  ~ FACTORA + Condition(EX.pcnm.sub), data=env.df)) # add + Condition(env.df$SAMPLE_ID) for repeated measures
# paste output here











  ## Mantel test -------

bray_not_na <- phyloseq::distance(physeq = EX_ps_clean.rar, method = "bray")

# make a data frame that can be used in vegan from the sample_data in the phyloseq object
meta.df <- data.frame(sample_data(EX_ps_clean.rar))

# check the columns in your meta dataframe, need numeric data for mantel test
str(meta.df)

# make numeric any columns that are a number but as listed as "chr", for example 
meta.df$DNA_extraction_batch <- as.numeric(meta.df$DNA_extraction_batch) # CHANGE ME to any columns that should be numeric

# replace any binary factors with numbers, for example
meta.df["Diet"][meta.df["Diet"] == "LooseAlfalfa"] <- "1"
meta.df["Diet"][meta.df["Diet"] == "PelletedAlfalfa"] <- "2"
meta.df$Diet <- as.numeric(meta.df$Diet)

#remove any columns that you can't make numeric
meta.df <- subset(meta.df, select = -c(Sample_type, Sheep_ID, Treatment))

# Make a distance matrix for your meta data
bray.dist.meta <- vegdist(meta.df, method = "bray") 

# Mantel test to see if SV and meta distance matrixes correlate
bray.mantel <- mantel(bray_not_na, bray.dist.meta, permutations = 999)

# look at the output
bray.mantel

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = bray_not_na, ydis = bray.dist.meta, permutations = 999) 
# 
# Mantel statistic r: 0.06023 
#       Significance: 0.224 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.108 0.143 0.166 0.189 
# Permutation: free
# Number of permutations: 999


# Biogeography (Lab 14)-----

# Time series Generalized Additive Model (Lab 14)-----

#example alpha diversity plot code with GAM line added ## CHANGE ME
plot_richness(no_wild, x="Date_sample_collected", measures="Observed", color="Sample_type") +
  geom_point(aes(shape=Sample_type, size = as.numeric(Days_since_full_moon))) + 
  scale_color_manual(values = nowild_colors) + scale_shape_manual(values = shapes) + 
  ylim(0,300) + 
  geom_labelsmooth(label = "sin(Date)", text_smoothing = 30, method = "lm", formula = y ~ sin(x), col="darkgoldenrod", alpha = 0.85, linetype = 2, boxlinewidth = 0.3) +
  geom_labelsmooth(label = "GAM(Date)", text_smoothing = 30, method = "gam", formula = y ~ s(x, bs = "cs"), col="blueviolet", alpha = 0.75, linetype = 1, boxlinewidth = 0.3) +
  theme_minimal() + ylab("Observed Bacterial Richness (SVs)") + theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_blank()) +xlab("Date Collected in 2021") + labs(shape = "Sample Type", color = "Sample Type", size = "Days since full moon") + theme(legend.text = element_text(size=11)) 


# run a GAM model

library(mgcv)
library(lubridate)

# correct the time format on your dates so R can read them
Vibrio_tankonly.df$Date <- as.Date(Vibrio_tankonly.df$Date_sample_collected) #CHANGE ME


mm3 <- gam(Observed ~ Clean_or_dirty 
            + Tank_setup 
            + s(Days_from_start)
            + Days_since_full_moon
            + s(Lunar_phase_percent),
            data = subset(Vibrio_tankonly.df, Tank_setup != "n/a" & Sequencing_batch=="1"), na.action = na.omit)

summary(mm3)

plot.gam(mm3)



# Niche/neutral models (Lab 15)-----
# Check out the online version: https://web.rniapps.net/webglv/
# There is an extra code packet for the course to get you started if you want to perform this in R. 

# Prep your sequence table for online L-V
# Make a dataframe of training data with OTUs as column and samples as rows, which is the phyloseq OTU table
seqtab_LV <- otu_table(EX_ps_clean)

# output phyloseq tax table as a dataframe to make it manipulable
tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

## if the Genus is empty, replace with the Family
tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

# bind Genus and Species together
tax.df$Genus.species <- paste(tax.df$Genus, tax.df$Species)

# set column of combined genus and species names as the column names for the predictors, replacing the full SV
colnames(seqtab_LV) <- tax.df$Genus.species

# clean up some of the other taxa info
colnames(seqtab_LV) = gsub("_unclassified", "", colnames(seqtab_LV))
colnames(seqtab_LV) = gsub("_Intercertae_Sedis", "", colnames(seqtab_LV))

# Make labelled SVs into a datatable
seqtab_LV.data <- data.frame(seqtab_LV)

# write to a TSV file which you can upload online
write.table(seqtab_LV.data, file='seqtab_LV.tsv', quote=FALSE, sep='\t', col.names = NA)



  # Source tracking (Lab 16) -----
# Source Tracker was developed for QIIME, and can identify samples which could be the “source” of the microbes in another community, called the “sink”. There was a version of the software coded for R, and the original code and example files can be found here: https://github.com/danknights/sourcetracker.

# IMPORTANT! I could not get the code to run following the walkthrough, until I went into the source code and hashed out all lines which rarefied the data. I made a modified version of SourceTracker.R to share with the class, and the rest of the code below  is a modified version of this walkthrough to go from phyloseq to sourcetracker in R: https://mgaley-004.github.io/MiSeq-Analysis/Tutorials/SourceSink.html 

# download the SourceTracker_modified.R that Sue shared

library(ape)
library(plyr)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(tidyr)
library(vegan)
library(scales)
library(grid)
library(reshape2)
library(gridExtra)
library(phyloseq)

# need to grab the OTU_table (SV table) from your phyloseq object
ex_seqtab <- as.data.frame(otu_table(EX_ps_clean))

# take the taxa names from the tax_table and paste onto the colnames for the seqtab (otherwise they will just be the full sequences)
  # output phyloseq tax table as a dataframe to make it manipulable
  tax.df <- data.frame(tax_table(EX_ps_clean), stringsAsFactors=FALSE)

  ## if the Genus is empty, replace with the Family
  tax.df$Genus[is.na(tax.df$Genus)] <- tax.df$Family[is.na(tax.df$Genus)]

  # bind Genus and Species together
  tax.df$Genus.species <- paste0(tax.df$Genus, "_", tax.df$Species)

  # set column of combined genus and species names as the column names for the predictors, replacing the full SV
  colnames(ex_seqtab) <- tax.df$Genus.species

# go from dataframe format to sourcetracker preferred format
ex_seqtab_st <- t(as.matrix(ex_seqtab))

# grab the meta
ex_meta <- as.data.frame(sample_data(EX_ps_clean))

# create a column that indicates 'source' and 'sink' if you don't already have it in there
ex_meta$SourceSink <- NA

# designate week 0 as 'source' and weeks 1 and 2 as 'sinks' using and if/else statement. If Week matches 0, will add 'source' to the new column, otherwise it will add 'sink'
ex_meta$SourceSink <-  ifelse(ex_meta$Week == "0", "source", "sink") #CHANGE ME

# indicate which samples you think are the source (the donor) and which you think are the sink (the receiver)
train.ix <- which(ex_meta$SourceSink=='source')
test.ix <- which(ex_meta$SourceSink=='sink')

envs <- ex_meta$Diet # usually a factor that indicates location 
desc <- ex_meta$Diet # usually a factor that indicates location 

# Change this line to reflect wherever you downloaded the modified sourcetracker code to
source('C:/Users/sueis/OneDrive/Documents/Teaching/AVS 454-554 DNA Sequencing Analysis Lab/R/SourceTracker_modified.r') #CHANGE ME

# tune the alpha values using cross-validation (this is slow! like very slow!) Use 0.001 if pressed for time. NOTE: I could not get this next command to work so I set my alpha values below.
# tune.results <- tune.st(ex_seqtab_st[train.ix,], envs[train.ix])

alpha1 <- 0.05 #tune.results$best.alpha1
alpha2 <- 0.05 #tune.results$best.alpha2

# Train sourcetracker using the source samples. Be unafraid of an error message that says ‘the condition has length > 1 and only the first element will be used.’ Seems to be an internal bug that doesn’t affect the outcome.
st <- sourcetracker(ex_seqtab_st[train.ix,], envs[train.ix])

#Estimate source proportions in sink samples. #BOOKMARK, sometimes throws errors
results <- predict(st, ex_seqtab_st[test.ix,], alpha1=alpha1, alpha2=alpha2)

# This part validates the 'source' models we made using the original training data.
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results with legend
plot(results, type='pie', include.legend=TRUE)
plot(results.train, type='pie', include.legend=TRUE)

# Instead of using those plots, you can extract the results like so and then make all your plots in either base R or ggplot just like you’ve already learned to do.

downstream <- data.frame(results$proportions)
downstream$id <- row.names(downstream)
meltdown <- melt(downstream, id.vars=c("id"))

ex_meta

## BOOKMARK did not finish fancying this up
colnames(meltdown) <- c("Bacterial_SV", "Diet", "value")
meltdown$Diet <- ex_meta[as.vector(meltdown$id), "Diet"]
#meltdown$Week <- ex_meta[as.vector(meltdown$id), "Week"]

meltdown$Month <- as.factor(sapply(meltdown$Date, function(x) strsplit(x, "/")[[1]][1]))
levels(meltdown$Month) <- c("June", "July", "August")
meltdown$Day<- as.factor(sapply(meltdown$Date, function(x) strsplit(x, "/")[[1]][2]))

meltdown$Rep <- as.factor(sapply(meltdown$id, function(x) strsplit(x, "_")[[1]][3]))
levels(meltdown$Rep) <- c("A", "AA", "B")

meltdown$Category <- as.factor(meltdown$Category)
levels(meltdown$Category) <- c("B", "E")

meltdown$SID <- as.factor(sapply(meltdown$id, function(x) strsplit(x, "_")[[1]][5]))

meltdown$TypeRep <- paste(meltdown$Category, meltdown$Rep, meltdown$Day, sep=".")
theme_set(theme_classic(base_size=10, base_family="Avenir"))

ggplot(meltdown, aes(x=Week, y=value, fill=Diet)) + geom_bar(stat="identity") + facet_wrap(~Downstream.Site+Month, scales="free_x") + xlab("")










# Dendogram (optional) ---------------------------------------
## Note: I typically use the UPGMA algorothm for Dendograms when I make them as stand-alone.  Some other graph functions will intergrate dendograms themselves.

#calculate dissimilarity based on your sequence data in base R
dist.object <- dist(as.data.frame(otu_table(EX_ps_clean.rar)), method = "euclidean", diag = FALSE, upper = FALSE, p = 2) # method = "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

# cluster the dissimiliary info with base R: https://www.rdocumentation.org/packages/stats/versions/3.2.1/topics/hclust
EX_ps.clust <- hclust(dist.object, method = "average", members = NULL) #method = average for UPGMA

# visualization in phyloseq: https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-analysis.html#hierarchical-clustering
plot(EX_ps.clust)


