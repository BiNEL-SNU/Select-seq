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

####Read Genome Database####
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v75)
library(org.Hs.eg.db)
RefeLibrary = EnsDb.Hsapiens.v86

####Define targetgene####
targetGene = c("CD274","PDCD1LG2","SIGLEC10")
filename="/home/taliq/2020_RNA/210615_SIGLEC10_PDL1_PDL2.png"

#CD2, CD3 ,CD10==MME, CD16==FCGR, CD18==ITGB, CD31==PECAM, CD64==FCGR, and CD140b==PDGFR#

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
BasicFolder = "/home/taliq/2020_RNA/Sample/"
OriginFolder = "201227_BRTissue1"
ReadData1=ReadData(BasicFolder,OriginFolder)

####File Read####
BasicFolder = "/home/taliq/2020_RNA/Sample/"
OriginFolder = "201227_BRTissue2"
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

####Select 0610 Tissue####
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

####DESeq2####
Sample_condition = sortedTypeList
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

####Draw heatmap####
select <- tempSymbol %in% targetGene
existgene = targetGene[targetGene %in% tempSymbol]
Normalized_count_select = log10(Normalized_count[select,]+1)
Targeted_symbol = tempSymbol[select]
Targeted_group = tempGroupList
Targeted_name = tempNameList
tempcount = Normalized_count_select[Targeted_symbol == targetGene[1]]
for (i in 2:length(targetGene)){
  tempcount = rbind(tempcount,Normalized_count_select[Targeted_symbol == targetGene[i]])
}
bluecolor="#619DFF"
redcolor="#F9776E"
greencolor="#01BA38"
yellowcolor="#FBCF35"

rownames(tempcount) = targetGene
colnames(tempcount) = colnames(Normalized_count_select)
Target1 = tempcount[,Targeted_group == "CD_44"]
Target2 = tempcount[,Targeted_group == "ALDH1"]
Target3 = tempcount[,Targeted_group == "DN"]
Target4 = tempcount[,Targeted_group == "DP"]
finalcount = cbind(Target4, Target2, Target1, Target3)
finalgroup = c(rep("DP", length(colnames(Target4))),
               rep("ALDH1", length(colnames(Target2))),
               rep("CD_44", length(colnames(Target1))),
               rep("DN", length(colnames(Target3))))

annotation_df = data.frame(value=c(finalgroup), row.names = colnames(finalcount))
aka3 = list(value = c(CD_44 = greencolor, ALDH1=redcolor, DN=bluecolor, DP=yellowcolor))
png(filename, width=1000, height=400)
pheatmap(finalcount,cluster_rows = FALSE,cluster_cols=FALSE,show_rownames = TRUE, 
         annotation_col=annotation_df,annotation_colors = aka3)
dev.off()

