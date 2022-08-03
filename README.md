# AneufinderCopyNumberCounter

AneufinderCopyNumberCounter contains scripts to quickly fetch the number of whole chromosome and partial copy number alterations from Aneufinder model files. It exports the counts to .xlxs files and prints the results in basic graphs to a PDF file. Copy number counting is based on the deviation from a predefined copy number baseline, e.g. n=2 for diploid cells and n=1 for haploid. If a single chromosome contains multiple segments, the copy number of each segment is considered separately. Small partial alterations can be excluded from counting by setting a minimum segment length.


## Scripts

To use the code you need two included scripts:

1. __RUN_AneufinderCopyNumberCounter__  
_Set input/output folders and choose correct settings_

2. __FUNC_AneufinderCopyNumberCounter__  
_Contains the code that is used to count the copy number alterations for each cell_


## Required R packages

The script makes use of the following R packages:

1. __[Aneufinder](https://bioconductor.org/packages/release/bioc/html/AneuFinder.html)__
1. __[tidyr](https://tidyr.tidyverse.org/)__
1. __[xlsx](https://github.com/colearendt/xlsx)__

## Input / Output

__Input:__  
* Aneufinder Model files  

__Output:__
* for each sample a .xlxs file with total number of whole chromosome and partial copy number alterations for each chromosome
* .pdf with total number of copy number alterations represented in a basic bar graph

## Instructions for use

To get started I would advise users to make use of R studio and create a new project in R and name it 'AneufinderCopyNumberCounter'. Then download both scripts from GitHub and place these within the project folder of your new AneufinderCopyNumberCounter project.

In general there's no need to open and/or adjust the function script, this is only needed if you like to make adjustments to the code that performs the actual counting or the code that creates the output files. You only need to make sure that the RUN script contains the correct source-path to the FUNC script.

The 'RUN_AneufinderCopyNumberCounter'-script is subdivided in multiple sections to create a good overview of the different settings. Prior to each run you probably like to give your project a new name, assign the correct input folder and check the settings.

After making all required adjustments, run the code line-by-line. The actual counting is initiated at the end of the RUN script by running  __AneufinderCopyNumberCounter(sampleIDs)__.

## Options

1. Set copy number base level for all and/or indvidual autosomes.
2. Set copy number base level for X chromosome.
3. Set minimum segment size for partial copy number alterations
4. Determine presence/absence of copy number alterations irrespective of exact deviation from baseline or include deviation from copy number base level and thereby count absolute number of alterations
5. Collect files from Aneufinder _or_ after filtering with [AneufinderFileFilter](https://github.com/TWvR/AneufinderFileFilter)

## Final comments

If you have any questions or need help with running the script, please don't hesitate to send me a message.

Thomas van Ravesteyn
