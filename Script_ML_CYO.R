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
# Outliers are values that are further than 1.5 * ICR (Inter-quartile range) from the closest hinge
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

# As senn from the plots and the dataset, there is a lot of predictors, where many may not be predictive in any way
# Therefore a dimensionality reduction method will follow to remove those genes that may not be predictive or 
# add to the variance of the data much. 

pca_rna_seq_train <- prcomp(x = rna_seq_train[,3:ncol(rna_seq_train)], center = TRUE) # Performs a PCA analysis on all the genes present in the dataset

str(pca_rna_seq_train) # Gives and indication of the structure and values after the pca analysis
dim(pca_rna_seq_train$rotation) # Dimensions shows the amount of Principal Components, 20531 genes and 638 principal components

pca_rna_seq_train$rotation[1:5,1:7] # Shows the first few PCs

plot(pca_rna_seq_train$x[],pca_rna_seq_train$x[])

pca_var <- pca_rna_seq_train$sdev^2 # Computes the variance of each PC
pca_var[1:10] # Shows the top 10 PC's variance

# For dimensionality reduction, we are interested in those PC's that inherently explains the most variance in the dataset
# therefore we should determine a cutoff point for the amount of variance included by the PC's

prop_var <- pca_var/sum(pca_var) # calculates the proportion each PC adds to the total variance

cumsum(prop_var[1:10]) # Shows the variance explained by the first ten PC's
                      # We already see that the first 10 PC's explain about 55 % of the variance in the dataset

plot(x = 1:length(prop_var), y = cumsum(prop_var), main = "Proportion Of Variance Explained By The Principal Components",
                                                  xlab = "Principal Component Number",
                                                  ylab = "Proportion Of Variance" ) # Plot showing how the different Principal components add to the total variance

min(which(cumsum(prop_var) > 0.8)) # The PC number where the total variance explained equals 80%

min(which(cumsum(prop_var) > 0.9)) # The PC number where the total variance explained equals 90%

min(which(cumsum(prop_var) > 0.95)) # The PC number where the total variance explained equals 95 %

ncol(pca_rna_seq_train$rotation) # Shows the total amount of PC's

min(which(cumsum(prop_var) > 0.8))/ncol(pca_rna_seq_train$rotation) # proportion of predictors that explain 80 % of variance

min(which(cumsum(prop_var) > 0.9))/ncol(pca_rna_seq_train$rotation) # proportion of predictors that explain 90% of variance

min(which(cumsum(prop_var) > 0.95))/ncol(pca_rna_seq_train$rotation) # proportion of predictors that explain 95 % of variance

# The cutoff value of the total variance included by the PC's should be carefully decided.
# We would like to include as much of the predictors as possible, while being effecient in using computer processing power
# For 90% of the variance explained we reduce the predictors by aprox. a half while
# for 80 % of the variance explained the predictors reduced to aprox. a fifth 

# For training the ML algorithm, both the 80% and 90 % variance cutoff values will be used and 
# accuracy measured.








