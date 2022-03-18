rm(list = ls())
library(DESeq2)
library(stringr)
library(readr)
library(tximportData)
library(readxl)
library(dplyr)
library(ggcorrplot)
library(pheatmap)
library(ggplot2)
library(xlsx)
library(Rtsne)
library(RColorBrewer)
library(EnhancedVolcano)

####Read Genome Database####
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
RefeLibrary = EnsDb.Hsapiens.v86

####Define Function####
ff <- function(input, TargetFolder){
  CountsFile = paste(TargetFolder, input, sep="")
  data = as.matrix(read.table(CountsFile, header = TRUE, row.names = 'Geneid')[6])
}
ReadData <- function(BasicFolder, OriginFolder){
  FinalOriFolder = paste(BasicFolder,OriginFolder, sep="")
  TargetFolder = paste(BasicFolder,OriginFolder,'/FeatureCounts/', sep="")
  CountsFile = list.files(TargetFolder)
  temp = str_detect(CountsFile, "summary") 
  CountsFile = CountsFile[!temp]
  Tempfile = paste(TargetFolder, CountsFile[1], sep="")
  TempGene = as.character(read.table(Tempfile, header = TRUE)$Geneid)
  
  Symbol = ensembldb::select(RefeLibrary, keys = TempGene, column = "SYMBOL", keytype = "GENEID", multiVals = "first")
  Symbol = Symbol$SYMBOL
  print("Read File")
  a = lapply(CountsFile, ff, TargetFolder = TargetFolder)
  print("Done")
  
  return(list(a, CountsFile, Symbol))
}

####File Read####
BasicFolder = "./"
OriginFolder = "BRTissue1"
ReadData1=ReadData(BasicFolder,OriginFolder)

####File Read####
BasicFolder = "./"
OriginFolder = "BRTissue2"
ReadData2=ReadData(BasicFolder,OriginFolder)

####Merge the count file####
CountsFile = c(ReadData1[[2]],ReadData2[[2]])
a = c(ReadData1[[1]],ReadData2[[1]])
Symbol=ReadData1[[3]]

StrSplit = str_split(CountsFile, '-')
NameList = StrSplit[[1]][1]
for (i in 2:length(StrSplit)){
  NameList = c(NameList, StrSplit[[i]][1])
}
TypeList = StrSplit[[1]][2]
for (i in 2:length(StrSplit)){
  TypeList = c(TypeList, StrSplit[[i]][2])
}
GroupList = StrSplit[[1]][3]
for (i in 2:length(StrSplit)){
  GroupList = c(GroupList, StrSplit[[i]][3])
}

SelectIndex = c(1:96, 103:112)


TargetA = a[SelectIndex]
TargetCountsFile = CountsFile[SelectIndex]
colnames(TargetA[[1]]) = TargetCountsFile[[1]]
Count = TargetA[[1]]
for (i in 2:length(TargetA)){
  colnames(TargetA[[i]]) = TargetCountsFile[[i]]
  Count = cbind(Count, TargetA[[i]])
}
EntrezID = rownames(Count)
StrSplit = str_split(TargetCountsFile, '-')
TargetNameList = StrSplit[[1]][1]
for (i in 2:length(StrSplit)){
  TargetNameList = c(TargetNameList, StrSplit[[i]][1])
}
TargetTypeList = StrSplit[[1]][2]
for (i in 2:length(StrSplit)){
  TargetTypeList = c(TargetTypeList, StrSplit[[i]][2])
}
TargetGroupList = StrSplit[[1]][3]
for (i in 2:length(StrSplit)){
  TargetGroupList = c(TargetGroupList, StrSplit[[i]][3])
}

####Sort file order####
rankorder = str_order(colnames(Count), numeric = TRUE)
sortedCount = Count[,rankorder]
sortedNameList = TargetNameList[rankorder]
sortedTypeList = TargetTypeList[rankorder]
sortedGroupList = TargetGroupList[rankorder]
FinalQualityFilter = c(1:80, 94:106)
sortedCount = Count[,FinalQualityFilter]
sortedNameList = TargetNameList[FinalQualityFilter]
sortedTypeList = TargetTypeList[FinalQualityFilter]
sortedGroupList = TargetGroupList[FinalQualityFilter]
AtoIGroup = c(rep("None",93))
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_25")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_73")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_75")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_77")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_78")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_79")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_80")] = "AtoI"
AtoIGroup[stri_detect_fixed(sortedNameList, "200603_103")] = "AtoI"

####DESeq2####
Sample_condition = AtoIGroup
Sample_set_type = sortedGroupList
coldata = data.frame(condition = Sample_condition, type = Sample_set_type)
dds = DESeqDataSetFromMatrix(sortedCount, colData = coldata, design = ~ condition)
keep = colSums(counts(dds) > 0.5) > 1000
dds = dds[,keep]
tempGroupList = sortedTypeList[keep]
tempNameList = sortedNameList[keep]
tempSymbol = Symbol
dds = DESeq(dds)
Normalized_count = counts(dds, normalized=TRUE)
rownames(Normalized_count) = tempSymbol

res = results(dds, contrast=c("condition","AtoI","None"))
resLFC = lfcShrink(dds,contrast=c("condition","AtoI","None"),res=res)

resOrdered <- resLFC[order(abs(resLFC$log2FoldChange),resLFC$padj,resLFC$baseMean, decreasing = TRUE),]
resOrdered$Symbol = Symbol

EnhancedVolcano(res, lab=Symbol, x='log2FoldChange',y='pvalue')

