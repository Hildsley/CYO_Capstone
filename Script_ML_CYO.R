library(tidyverse)
library(caret)
library(R.utils)

#############################################################
###   Creating the dataset for the Capstone CYO Project   ###
#############################################################

download.file(url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00401/TCGA-PANCAN-HiSeq-801x20531.tar.gz",destfile = "Data/RNA_SeqData.tar.gz")  # Downloads the dataset

untar(tarfile = "Data/RNA_SeqData.tar.gz",exdir = "Data/") #Unzips the dataset as was downloaded from uci

myfile_data <- read_csv(file = "Data/TCGA-PANCAN-HiSeq-801x20531/data.csv") # reads the csv file from the data file in the folder, contains rna seq data
myfile_labels_rows <- read_csv(file = "Data/TCGA-PANCAN-HiSeq-801x20531/labels.csv") # Reads the csv file that contains the information regarding the types of samples

colnames(myfile_data)[1] <- "sample_number" # Changes the X1 column name to sample_number
colnames(myfile_data[1:5]) # Show the change

colnames(myfile_labels_rows)[1] <- "sample_number" # Changes the X1 column name to sample_number, same as myfile_data for tidying data later
colnames(myfile_labels_rows[1:2]) # Show the Change


# Change the name of the file to rna_seq_dat, indicates the type of information
# Add the Class of the labels to the data.csv
rna_seq_dat <- data.frame(myfile_labels_rows) %>% left_join(myfile_data, by = "sample_number") # Adds the classes of the samples together with the rna seq data into one data frame
dim(rna_seq_dat) # show that there are 801 rows and 20533 columns, 20531 genes is present

str(rna_seq_dat) # Gives the structure of the dataset
head(rownames(rna_seq_dat), n = 5) # Shows the rows' names of first 10

head(colnames(rna_seq_dat), n = 5) # Shows first 10 columns' names

head(rna_seq_dat$Class) # Shows the classes of samples present in the dataset

levels(rna_seq_dat$Class) # This code shows us that the Class column which consists of BRCA, PRAD, LUAD, KIRC and COAD do not have those levels

rna_seq_dat$Class <- as.factor(x = rna_seq_dat$Class) # Will add the Levels present in the Class column

levels(rna_seq_dat$Class) # We now see the present levels, this will be the predictions the ML algorithm would need to predict given the predictors, in this case
                          # the genes' relative up- or downregulated readings.

##############################################################################
###  Splitting the dataset into two groups, a training and validation set.  ##
##############################################################################
set.seed(2019)
ind <- createDataPartition(y = rna_seq_dat$Class, times = 1,p = 0.2, list = FALSE) # creates the indexes for creating the training
                                                                                  # and validation datasets, 20 % of the original dataset
                                                                                  # will be allocated towards the validation set of samples.

rna_seq_val <- rna_seq_dat[ind,]  # Splices out the validation samples from the original dataset 
rna_seq_train <- rna_seq_dat[-ind,] # Splices out the training samples for ML from the original dataset

nrow(rna_seq_val) # Shows rows for the different spliced out sets
nrow(rna_seq_train)

nrow(rna_seq_val) + nrow(rna_seq_train) == nrow(rna_seq_dat) # Ensures no samples are lost during the data partitioning step


#################################### 
###   Exploratory Data analysis   ##
####################################


head(rna_seq_train$gene_44, n = 7) # The following code shows how the readings of the genes look like for different samples
head(rna_seq_train$gene_1838, n =7)

#Plot showing the points and the the five number summary of the readings of gene_1, which may have some predictive power
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_1)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 1") +
  ggtitle(label = "Average Reading For The Different Samples For Gene 1" ) +
  theme(legend.title.align = 0.5) +
  geom_jitter(alpha = 0.3, color = "dark blue")
  

# Plot showing how some genes' readings may have no predictive power 
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_0)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 0") +
  ggtitle(label = "Average Reading For The Different Samples For Gene 0" )
  
#Plot showing the points and the the five number summary of the readings of gene_1072, which may have a good predictive power
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_1072)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 1072") +
  ggtitle(label = "Average Reading For The Different Samples For Gene 1072" ) +
  theme(legend.title.align = 0.5) +
  geom_jitter(alpha = 0.3, color = "dark blue")

 






















