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
colnames(white_data3) = c(" ", "TCGA-E2-A15E",	"TCGA-D8-A1JM",	"TCGA-AN-A0AT",	"TCGA-EW-A3E8",	"TCGA-BH-A0E0",	"TCGA-D8-A1JS",	"TCGA-E2-A2P5",	"TCGA-D8-A1J9",	"TCGA-BH-A0H5",	"TCGA-E2-A56Z")

            
#make the dataframe into a matrix (should only have numbers)
white_data_matrix <-data.matrix(white_data3)
colnames(white_data_matrix) = c(" ", "TCGA-E2-A15E",	"TCGA-D8-A1JM",	"TCGA-AN-A0AT",	"TCGA-EW-A3E8",	"TCGA-BH-A0E0",	"TCGA-D8-A1JS",	"TCGA-E2-A2P5",	"TCGA-D8-A1J9",	"TCGA-BH-A0H5",	"TCGA-E2-A56Z")
print(white_data_matrix)
#create metaData


countData<- data.frame(white_data_matrix)
metaData<- data.frame(white_data_matrix)
dds1<- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              )
###################################potential new file

countData <- data.frame(read.csv('Sample RNAseq w Clinical Data.csv', header = T, sep = ','))
countData<- subset(countData, select=-c(1))
colnames(countData) = c(" ", "TCGA-E2-A15E",	"TCGA-D8-A1JM",	"TCGA-AN-A0AT",	"TCGA-EW-A3E8",	"TCGA-BH-A0E0",	"TCGA-D8-A1JS",	"TCGA-E2-A2P5",	"TCGA-D8-A1J9",	"TCGA-BH-A0H5",	"TCGA-E2-A56Z")
#missing gene name column when converted to matrix (not all numbers)
countData2<- subset(countData, select=-c(1))
countDataMatrix<- data.matrix(countData2)
metaData<- data.frame(read.csv('Just Clinical Data.csv', header = T, sep = ','))
metaData1 <- metaData[8,]
colnames(metaData1) = c("status", "status"," ", "TCGA-E2-A15E",	"TCGA-D8-A1JM",	"TCGA-AN-A0AT",	"TCGA-EW-A3E8",	"TCGA-BH-A0E0",	"TCGA-D8-A1JS",	"TCGA-E2-A2P5",	"TCGA-D8-A1J9",	"TCGA-BH-A0H5",	"TCGA-E2-A56Z")
metaData2 <- subset(metaData1, select=-c(1,2,3))
metaData3 <- t(metaData2)
colnames(metaData3)=c("condition")
#if you want to add the gene names back, include the lines below (adds a "buffer" row)
#buffer<- c(" ")
#metaData4<-rbind(buffer, metaData3)
#rownames(metaData4)=c(" ", "TCGA-E2-A15E",	"TCGA-D8-A1JM",	"TCGA-AN-A0AT",	"TCGA-EW-A3E8",	"TCGA-BH-A0E0",	"TCGA-D8-A1JS",	"TCGA-E2-A2P5",	"TCGA-D8-A1J9",	"TCGA-BH-A0H5",	"TCGA-E2-A56Z")
metaDataMatrix<- data.matrix(metaData3)
all(rownames(metaDataMatrix)==colnames(countDataMatrix))
idx<-match(colnames(countData), rownames(metaData3))
reordered_metaData<-(metaData4[idx,])
print(reordered_metaData)
all(rownames(reordered_metaData)==colnames(countData2))


dds <- DESeqDataSetFromMatrix(countData=countDataMatrix, 
                              colData=metaDataMatrix, 
                              design = ~ condition, tidy = T)
dds
#variance stabilizing transformation
#heatmap
white_data_heatmap <- pheatmap(white_data_matrix, annotation_names_col = T, margins=c(5,10), scale="row", fontsize=9, angle_col = 45, display_numbers = T, fontsize_number = 6, treeheight_row = 0, treeheight_col = 0, main = "Potential Molecular Targets for TNBC")
print(white_data_heatmap)

#MA plot
maplot <- function (white_data_matrix, thresh=0.05, labelsig=TRUE, textcx=1, ...)
with(white_data_matrix, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
with(subset(white_data_matrix, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
if (labelsig) {
  }require(calibrate)
  with(subset(white_data_matrix, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))  

png(filename = "maplot.png", 480, 480, pointsize=20)
maplot(white_data_matrix, main="MA Plot") dev.off()

#Mean and Variance plot
#PCA

  