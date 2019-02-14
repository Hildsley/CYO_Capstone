if (!require(devtools))install.packages("devtools")
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(caret)) install.packages("caret")
if (!require(R.utils)) install.packages("R.utils")
library(devtools)
install_github("eddelbuettel/rbenchmark")

library(tidyverse)
library(caret)
library(R.utils)
library(rbenchmark)
#############################################################
###   Creating the dataset for the Capstone CYO Project   ###
#############################################################

download.file(url = "https://archive.ics.uci.edu/ml/machine-learning-databases/00401/TCGA-PANCAN-HiSeq-801x20531.tar.gz",
              destfile = "RNA_SeqData.tar.gz")  # Downloads the dataset

untar(tarfile = "RNA_SeqData.tar.gz") #Unzips the dataset as was downloaded from uci

myfile_data <- read_csv(file = "TCGA-PANCAN-HiSeq-801x20531/data.csv") # reads the csv file from the data file in the folder, contains rna seq data
myfile_labels_rows <- read_csv(file = "TCGA-PANCAN-HiSeq-801x20531/labels.csv") # Reads the csv file that contains the information regarding the types of samples

colnames(myfile_data)[1] <- "sample_number" # Changes the X1 column name to sample_number
colnames(myfile_data[1:5]) # Show the change

colnames(myfile_labels_rows)[1] <- "sample_number" # Changes the X1 column name to sample_number, same as myfile_data for tidying data later
colnames(myfile_labels_rows[1:2]) # Show the Change

unlink(x = "TCGA-PANCAN-HiSeq-801x20531",recursive = TRUE) # Remove the directory not needed anymore
unlink("RNA_SeqData.tar.gz") # remove file not needed anymore
# Change the name of the file to rna_seq_dat, indicates the type of information
# Add the Class of the labels to the data.csv

rna_seq_dat <- data.frame(myfile_labels_rows) %>% left_join(myfile_data, by = "sample_number") # Adds the classes of the samples together with the rna seq data into one data frame
dim(rna_seq_dat) # show that there are 801 rows and 20533 columns, 20531 genes is present

str(rna_seq_dat) # Gives the structure of the dataset
head(rownames(rna_seq_dat), n = 5) # Shows the rows' names of first 10

rm(myfile_data,myfile_labels_rows) # removes the objects not needed any more

head(colnames(rna_seq_dat), n = 5) # Shows first 10 columns' names

head(rna_seq_dat$Class) # Shows the classes of samples present in the dataset

levels(rna_seq_dat$Class) # This code shows us that the Class column which consists of BRCA, PRAD, LUAD, KIRC and COAD do not have those levels

rna_seq_dat$Class <- as.factor(x = rna_seq_dat$Class) # Will add the Levels present in the Class column

levels(rna_seq_dat$Class) # We now see the present levels, this will be the predictions the ML algorithm would need to predict given the predictors, in this case
                          # the genes' relative up- or downregulated readings.

#############################################################################
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

rna_seq_train[1:6,1:6] %>% knitr::kable() # Shows the first 6 rows and column of the training set

#Plot showing the points and the the five number summary of the readings of gene_1, which may have some predictive power
# Outliers are values that are further than 1.5 * ICR (Inter-quartile range) from the closest hinge
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_1)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 1") +
  ggtitle(label = "Average Reading The Different Samples For Gene 1" ) +
  theme(legend.title.align = 0.5) +
  geom_jitter(alpha = 0.3, color = "dark blue")
  

# Plot showing how some genes' readings may have no predictive power 
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_0)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 0") +
  ggtitle(label = "Average Reading For Different Samples For Gene 0" )
  
#Plot showing the points and the the five number summary of the readings of gene_1072, which may have a good predictive power
rna_seq_train %>% group_by(Class) %>% 
  ggplot(aes(Class,gene_1072)) + geom_boxplot (color = "red", fill = "dark red") + 
  xlab(label = "Sample Type") +
  ylab(label = "Average Reading For Gene 1072") +
  ggtitle(label = "Average Reading For Different Samples For Gene 1072" ) +
  theme(legend.title.align = 0.5) +
  geom_jitter(alpha = 0.3, color = "dark blue")

# As seen from the plots and the dataset, there is a lot of predictors, where many may not have great predictive power 

##########################################################################################################################
# Therefore a dimensionality reduction method will follow to remove those genes that may not be predictive or 
# add to the variance of the data much. 


pca_rna_seq_train <- prcomp(x = rna_seq_train[,3:ncol(rna_seq_train)], center = TRUE) # Performs a PCA analysis on all the genes present in the dataset

str(pca_rna_seq_train) # Gives an indication of the structure and values after the pca analysis
dim(pca_rna_seq_train$rotation) # Dimensions shows the amount of Principal Components, 20531 genes and 638 principal components

pca_rna_seq_train$x[1:5,1:7] # Shows the first few PCs

pca_ggplot <- data.frame(rna_seq_train[,2] , pca_rna_seq_train$x)
colnames(pca_ggplot)[1] <- "Cancer_Type"
pca_ggplot[,1] <- as.factor(pca_ggplot[,1])

pca_ggplot %>% ggplot(aes(x = PC1, y = PC2, color = Cancer_Type)) +
  geom_point() +
  ggtitle(label = "Plot Of PC1 Versus PC2 For Training Dataset" ) +
  xlab(label ="Principal Component 1") +
  ylab(label = "Principal Component 2") # Plot showing PC1 and PC2 of classes with good clustering


pca_var <- pca_rna_seq_train$sdev^2 # Computes the variance of each PC
pca_var[1:10] # Shows the top 10 PC's variance

#########################################################################################################################
# For dimensionality reduction, we are interested in those PC's that inherently explains the most variance in the dataset
# therefore we should determine a cutoff point for the amount of variance included by the PC's

prop_var <- pca_var/sum(pca_var) # calculates the proportion each PC adds to the total variance

cumsum(prop_var[1:10]) # Shows the variance explained by the first ten PC's
                      # We already see that the first 10 PC's explain about 55 % of the variance in the dataset

plot(x = 1:length(prop_var), y = cumsum(prop_var), main = "Proportion Of Variance Explained",
                                                  xlab = "Principal Component Number",
                                                  ylab = "Proportion Of Variance" ) # Plot showing how the different Principal components add to the total variance

min(which(cumsum(prop_var) > 0.8)) # The PC number where the total variance explained equals 80%

min(which(cumsum(prop_var) > 0.9)) # The PC number where the total variance explained equals 90%

ncol(pca_rna_seq_train$rotation) # Shows the total amount of PC's

min(which(cumsum(prop_var) > 0.8))/ncol(pca_rna_seq_train$rotation) # proportion of predictors that explain 80 % of variance

min(which(cumsum(prop_var) > 0.9))/ncol(pca_rna_seq_train$rotation) # proportion of predictors that explain 90% of variance


########################################################################################################
# The cutoff value of the total variance included by the PC's should be carefully decided.
# We would like to include as much of the predictors as possible, while being effecient in using computer processing power
# For 90% of the variance explained we reduce the predictors by aprox. a half while
# for 80 % of the variance explained the predictors reduced to aprox. a fifth 

# For training the ML algorithm, both the 80% and 90 % variance cutoff values will be used and 
# accuracy measured.

pca_rna_seq_train_80 <- data.frame( rna_seq_train[,2],pca_rna_seq_train$x[,1:min(which(cumsum(prop_var)> 0.8))]) # Joins the classes of the samples with the respective
                                                                                                                # PCs after the pca, for the 80% variance PCs.

pca_rna_seq_train_90 <- data.frame( rna_seq_train[,2],pca_rna_seq_train$x[,1:min(which(cumsum(prop_var) > 0.9))]) # same joining as above for the 90% variance PCs

colnames(pca_rna_seq_train_80)[1] <- "sample_class" # change column name to sample class
colnames(pca_rna_seq_train_80)


colnames(pca_rna_seq_train_90)[1] <- "sample_class"
colnames(pca_rna_seq_train_90)

models <- c("bayesglm", "rpart","knn","svmLinear3","lda","naive_bayes","pls","snn","svmLinear","svmRadial") # Different models which will be trained

fitControl <- trainControl(method = "cv", number = 10 , p = 0.8) # Control ensuring a cross-validation of 10 times would be completed on the training set of 80%


fits_80 <- lapply(models, function(models){
  set.seed(2020)
  print(models)
  train(sample_class ~ . , data = pca_rna_seq_train_80, method = models , trControl = fitControl)
})

############################################################################################## 
# Loop printing the associated Accuracy values for the different models

method_80 <- matrix() # variable to store the methods
acc_80 <- matrix() # Variable to store the accuracies

for (x in 1:length(fits_80)) {
  method_80[x] <- (fits_80[[x]][1]$method)
  acc_80[x] <- (max(fits_80[[x]][4]$results$Accuracy))
}
acc_list_80 <- data.frame(method_80,acc_80)

acc_list_80 %>% knitr::kable() # Shows the accuracies for each model


# From the for loop it is seen that models svmLinear3 and lda had perfect accuracy
# while knn, svmlinear and avNNet had very high accuracies
# Therefore the svmlinear3 and lda models would be evaluated on the validation set

################# Same fitting is done for the 90%variance training set

fits_90 <- lapply(models, function(models){
  set.seed(2020)
  print(models)
  train(sample_class ~ . , data = pca_rna_seq_train_90, method = models , trControl = fitControl)
})


method_90 <- matrix()
acc_90 <- matrix() 

for (x in 1:length(fits_90)) {
  method_90[x] <- (fits_90[[x]][1]$method)
  acc_90[x] <- (max(fits_90[[x]][4]$results$Accuracy))
  }
acc_list_90 <- data.frame(method_90,acc_90)

acc_list_90 %>% knitr::kable() # Shows the accuracies for each model

# This piece of code evaluates if the more PC's dataset has better accuracies after the models have been trained.
sum(acc_list_80$acc_80 > acc_list_90$acc_90) # counts how many models have higher accuracies on the 80% variance set than the 90% variance set.

##########################################################################################################################################
# The accuracies between the 80% and 90% variance sets did not give any indication why the 90% variance dataset that contains more predictors 
# should be considered above the 80% variance training set. Therefore the 80% variance training models would be used for evaluating the accuracy
# on the validation set.

acc_list_80 %>% knitr::kable()

# The svmLinear 3 and lda models were perfectly accurate after training.
# The next consideration would be which model's computation time is the least.

time_svmL <- benchmark(train(sample_class ~ . , data = pca_rna_seq_train_80 , method = "svmLinear3", trControl = fitControl , tuneGrid = expand.grid(cost = 0.25, Loss = 1)),
                      columns = c("elapsed"), replications = 10) # runs the model fitting algorithm 10 times and records the time it takes to complete.

# The expand grid best tune variables was used as was determined in the original fits_80 piece of code 
# fits_80[[4]][6]

time_lda <- benchmark(train(sample_class ~ . , data = pca_rna_seq_train_80, method = "lda" , trControl = fitControl), columns = c("elapsed"),replications = 10 )

times <- data.frame(time_svmL) %>% rbind(time_lda)#  
rownames(times) <- c("svmL","lda")                # Shows that the computation time for lda is better
times                                             #
  
fit_lda <- train(sample_class ~ . , data = pca_rna_seq_train_80, method = "lda",trControl = fitControl) # fit the best model chosen


#############################################################################################################################################
# The validation dataset should first be transformed into principal components which contains the same rotations as the training dataset.


pca_rna_seq_val <- predict(pca_rna_seq_train,newdata = rna_seq_val) #Transform validation set into similar pca parameters of the training set

pca_rna_seq_val <- data.frame(pca_rna_seq_val[,1:min(which(cumsum(prop_var)> 0.8))] ) # removes the PCs we are not interested in and adds the sample's classes to the data frame

pca_rna_seq_val %>% ggplot(aes(x=PC1, y = PC2, color = rna_seq_val[,2])) +
  geom_point() +
  ggtitle(label = "PC1 versus PC2 Of The Validation Set") +
  xlab(label = "PC1") +
  ylab(label = "PC2")
  # Plot showing the similarity between the original plot of PCs of only the training dataset

pred_val <- predict(fit_lda,pca_rna_seq_val) # Predicts the classes of validation set given the model that is fitted to the training dataset

confusionMatrix(pred_val,rna_seq_val[,2]) # Confusion matrix showing how well the model performs

# The model perfecly predicted the classes of the validation set.




