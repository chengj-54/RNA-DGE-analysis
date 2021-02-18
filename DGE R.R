#install packages
install.packages("vegan")
install.packages("factoextra")
install.packages("DESeq2")
install.packages("pheatmap")
install.packages("dplyr")
install.packages("ggfortify")
install.packages("car")
#this package allows you to import data from xls or xlsx files
install.packages("readxl")
#load packages
library("vegan")
library("factoextra")
library("DESeq2")
library("pheatmap")
library("dplyr")
library("ggfortify")
library("car")
library("readxl")

#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
#download DESeq2 from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
#load DESeq2
library("DESeq2")

#start DGE analysis
#read in the xls/xlsx files from RNA seq data
test_data<- read_excel(file.choose())
view(test_data)
#create a copy of the file, do not change the original file data
test_data1<- read_excel(file.choose())
#create a race column and bind it to the dataset
race<-c("Race"," ", "white", "white", "white", "white", "white", "white", "black/African American", "white", "white", "black/African American")
race
data_race<-rbind(race, test_data1)
#reorder the data
data_race2 <- data_race[, c(1,2,3,4,5,6,7,8,10,11,9,12)]
#subset the data so that white patients undergo a separate DGE than the black patients
white_data<- subset(data_race2[,1:10])
black_data<- subset(data_race2[,11:12])
#remove first two columns
white_data2<- subset(white_data[,3:10])
#remove the first row
white_data3<- subset(white_data2[2:16,])
            
#make the dataframe into a matrix (should only have numbers)
white_data_matrix<-data.matrix(white_data3)
print(white_data_matrix)
#heatmap
white_data_heatmap <- pheatmap(white_data_matrix, annotation_names_col = F, margins=c(5,10), scale="column", fontsize=11, angle_col = 45, display_numbers = T, fontsize_number = 6, treeheight_row = 30, treeheight_col = 30)


