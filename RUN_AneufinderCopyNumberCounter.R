###-------------AneufinderCopyNumberCounter-----------###

#----Thomas van Ravesteyn
#----Kops Group
#----Hubrecht Institute

###--------------- Load dependencies -----------------###
###---------------------------------------------------###

#load libraries
library(AneuFinder)
library(tidyr)
library(openxlsx)

#Set source  file that contains AneufinderCopyNumberCounter functions
source("FUNC_AneufinderCopyNumberCounter.R")


###------------- Input / Output folders --------------###
###---------------------------------------------------###

#set input directory, each subfolder representing a sample which contains the Aneufinder Output
input_dir <- c("Input/...")

#set output directory
output_dir <- c("Output/...")


###--------------Settings before start----------------###
###---------------------------------------------------###

#Determine presence of CNAs irrespective of exact deviation from baseline ("relative") or 
#include deviation from copy number base level and count absolute number of alterations ("absolute")
analysis.type <- "relative"

#Set model from which output should be analyzed ("dnaCopy" / "hmm" / "edivisive") ! case-sensitive
selected.model <- "dnaCopy"

#set copy number baseline (e.g. for diploid cells "2", for haploid cells "1")
chr_n <- c(2)
 
#set copy number baseline for the X chromosome (e.g. for male derived cells "1", for female "2")
chrX_n <- c(2)

#[optional] set divergent copy number baseline for any one or more other chromosomes (e.g. chr1_n <- c(3))
#chr1_n <- c(3)
#chr11_n <- c(5)
#chr15_n <- c(3)


#set minimum segment size (number of basepairs) that should be included in the counting
segment_min_size <- 15e6


###--------------- Collect files to analyze ----------------###
###---------------------------------------------------------###

#Collect all sample names from input folder
sampleIDs <- list.files(path = input_dir, pattern = "")

#create empty objects
hmmFiles <- list()
dnaFiles <- list()
edivisiveFiles <- list()

###-------------------- CHOOSE TYPE OF FILE SOURCE ------------------------###

###-------- Collect Files from -- unfiltered -- Aneufinder output----------###
for(sample in sampleIDs) {
  hmmFiles[[sample]] <- list.files(paste0(input_dir, "/",sample, "/MODELS/method-HMM/"), full = TRUE)
  dnaFiles[[sample]] <- list.files(paste0(input_dir, "/",sample, "/MODELS/method-dnacopy/"), full = TRUE)
  edivisiveFiles[[sample]] <- list.files(paste0(input_dir, "/",sample, "/MODELS/method-edivisive/"), full = TRUE)
}

###--------------------------------- OR -----------------------------------###


###-------- Collect Files from -- AneufinderFileFilter -- output-----------###
for(sample in sampleIDs) {
  hmmFiles[[sample]] <- list.files(paste0(input_dir, "/", sample, "/QC_FILES_hmm/selected_model-files"), full = TRUE)
  dnaFiles[[sample]] <- list.files(paste0(input_dir, "/", sample, "/QC_FILES_dnaCopy/selected_model-files"), full = TRUE)
  edivisiveFiles[[sample]] <- list.files(paste0(input_dir, "/",sample, "/QC_FILES_edivisive/selected_model-files"), full = TRUE)
}



###----------Run Copy number counter-----------###
###--------------------------------------------###

AneufinderCopyNumberCounter(sampleIDs)
