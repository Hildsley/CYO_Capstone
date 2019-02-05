library(tidyverse)
library(caret)
library(R.utils)


download.file(url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00401/TCGA-PANCAN-HiSeq-801x20531.tar.gz",destfile = "Data/RNA_SeqData.tar.gz")

untar(tarfile = "Data/RNA_SeqData.tar.gz",exdir = "Data/")

myfile_data <- read_csv(file = "Data/TCGA-PANCAN-HiSeq-801x20531/data.csv") # reads the csv file from the data file in the folder
myfile_labels_rows <- read_csv(file = "Data/TCGA-PANCAN-HiSeq-801x20531/labels.csv") # Reads the csv file that contains the information regarding the types of samples

colnames(myfile_data)[1] <- "sample_number" # Changes the X1 column name to sample_number
colnames(myfile_data[1:5]) # Show the change

rna_seq_dat <- myfile_data # Change the name of the file to rna_seq_dat, indicates the type of information

dim(rna_seq_dat) # show that there are 801 rows and 20532 columns

head(rownames(rna_seq_dat), n = 10)
head(colnames(rna_seq_dat), n = 10)






























