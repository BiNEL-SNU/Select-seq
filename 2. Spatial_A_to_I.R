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
library(pvclust)
library(EnhancedVolcano)
library(dendextend)

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

TargetSeqNumList_split = str_split(TargetNameList, '_')
TargetSeqNumList = TargetSeqNumList_split[[1]][3]
for (i in 2:length(TargetSeqNumList_split)){
  TargetSeqNumList = c(TargetSeqNumList, TargetSeqNumList_split[[i]][3])
}

Metadata_path1 = './Metadata/Metadata_2012_BRCA_1.xlsx'
Metadata1 = read_excel(Metadata_path1)
colnames(Metadata1) <- c('Seq_num',colnames(Metadata1)[2:10])

Metadata_path2 = './Metadata/Metadata_2012_BRCA,Bone_2.xlsx'
Metadata2 = read_excel(Metadata_path2)
colnames(Metadata2) <- c('Seq_num',colnames(Metadata2)[2:9])
Metadata2_select = Metadata2[1:10,]

Metadata = rbind(Metadata1[,1:7], Metadata2_select[,1:7])

TargetPosList = Metadata[(Metadata['Name'] == TargetNameList[1]),]['Location']
for (i in 2:length(TargetSeqNumList)){
  TargetPosList = c(TargetPosList, Metadata[(Metadata['Name'] == TargetNameList[i]),]['Location'])
}
TargetPosList = unlist(TargetPosList, use.names = FALSE)

FinalFilter = rep(TRUE, length(count))
Add_filter_count = Count

Spatial_location='./Position/CD44_Stem_cell_output.txt'
Spatial_location=read.table(Spatial_location)

X1 = (levels(Spatial_location$V2)[Spatial_location$V2])[2:length(Spatial_location$V2)]
X1 = as.numeric(X1)
X2 = (levels(Spatial_location$V4)[Spatial_location$V4])[2:length(Spatial_location$V4)]
X2 = as.numeric(X2)
MeanXList = cbind(X1, X2)
MeanXList = rowMeans(MeanXList)

X1 = (levels(Spatial_location$V3)[Spatial_location$V3])[2:length(Spatial_location$V3)]
X1 = as.numeric(X1)
X2 = (levels(Spatial_location$V5)[Spatial_location$V5])[2:length(Spatial_location$V5)]
X2 = as.numeric(X2)
MeanYList = cbind(X1, X2)
MeanYList = rowMeans(MeanYList)

SpatialIndex = as.numeric((levels(Spatial_location$V1)[Spatial_location$V1])[2:length(Spatial_location$V1)])

TargetXList = MeanXList[SpatialIndex == TargetPosList[1]]
TargetYList = MeanYList[SpatialIndex == TargetPosList[1]]
for (i in 2:length(TargetPosList)){
  TargetXList = c(TargetXList, MeanXList[SpatialIndex == TargetPosList[i]])
  TargetYList = c(TargetYList, MeanYList[SpatialIndex == TargetPosList[i]])
}
TargetXList = TargetXList[FinalFilter]
TargetYList = TargetYList[FinalFilter]

Filtercount = Add_filter_count[rowMeans(Add_filter_count > 0) > 0.5, ]
PCAcount = prcomp(Filtercount)
Filterpos = TargetPosList[FinalFilter]
PCAresults = PCAcount$rotation
PCAresults = as.data.frame(PCAresults)
PCAresults$X_pos = TargetXList
PCAresults$Y_pos = TargetYList

Rownames_target = str_split(rownames(PCAresults), '_')
rowname_select = paste(str_split(Rownames_target[[1]][3], '-')[[1]][1], 
                       str_split(Rownames_target[[1]][3], '-')[[1]][2],
                       sep = '_')
for (i in 2:length(rownames(PCAresults))){
  rowname_select = c(rowname_select, paste(str_split(Rownames_target[[i]][3], '-')[[1]][1], 
                                           str_split(Rownames_target[[i]][3], '-')[[1]][2],
                                           sep = '_'))
}
rownames(PCAresults) = rowname_select

PCAresults = as.data.frame(t(as.matrix(PCAresults)))
fit <- pvclust(PCAresults, method.hclust="ward",
               method.dist="euclidean")


dend1 <- color_branches(fit$hclust, k = 8)
dend2 <- color_labels(dend1, k=8)
pdf('./Sample/Tree_label.pdf', width=16, height=8)
plot(dend2)
dev.off()

defined_group_8 = as.data.frame(cutree(fit$hclust, k = 8))
defined_group_8$pos = Filterpos
colnames(defined_group_8) = c("Group", "Pos")

Total_A_to_I_tissue1 = './AtoI_calling/201227_BRTissue1_editing_sites_annotated_v2.csv'
Total_A_to_I_tissue1 = read.csv(Total_A_to_I_tissue1)

Total_A_to_I_tissue2 = './AtoI_calling/201227_BRTissue2_editing_sites_annotated_v2.csv'
Total_A_to_I_tissue2 = read.csv(Total_A_to_I_tissue2)
Select_tissue = str_detect(Total_A_to_I_tissue2$Name, "200603")
Total_A_to_I_tissue2_select = Total_A_to_I_tissue2[Select_tissue,]

Total_A_to_I = rbind(Total_A_to_I_tissue1, Total_A_to_I_tissue2_select)

A_to_I_include = Total_A_to_I$Name %in% TargetNameList
Total_A_to_I_select = Total_A_to_I[A_to_I_include,]

i=1
Selected_position = defined_group_8[defined_group_8$Group == i,]$Pos
Selected_A_to_I = Total_A_to_I[Total_A_to_I$Location %in% Selected_position,]
Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 =
  levels(Selected_A_to_I$Gene.wgEncodeGencodeBasicV34)[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34]
Selected_unique_gene = Selected_A_to_I$Gene.wgEncodeGencodeBasicV34
#Selected_unique_gene = levels(Selected_unique_gene)[Selected_unique_gene]
Selected_unique_gene = unique(Selected_unique_gene)

Selected_unique_gene_length = length(unique(Selected_A_to_I[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 == Selected_unique_gene[1], ]$Location)) / length(Selected_position)
for (j in 2:length(Selected_unique_gene)){
  Selected_unique_gene_length = c(Selected_unique_gene_length, 
                                  length(unique(Selected_A_to_I[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 == Selected_unique_gene[j], ]$Location)) / length(Selected_position))
}

Final_gene = Selected_unique_gene[Selected_unique_gene_length > 0.1]

Selected_group = rep(i, length(Final_gene))
Select_df = data.frame(Selected_group,
                       Final_gene)
colnames(Select_df)=c("Group","Gene")

for (i in 2:8){
  Selected_position = defined_group_8[defined_group_8$Group == i,]$Pos
  Selected_A_to_I = Total_A_to_I[Total_A_to_I$Location %in% Selected_position,]
  Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 =
    levels(Selected_A_to_I$Gene.wgEncodeGencodeBasicV34)[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34]
  Selected_unique_gene = Selected_A_to_I$Gene.wgEncodeGencodeBasicV34
  #Selected_unique_gene = levels(Selected_unique_gene)[Selected_unique_gene]
  Selected_unique_gene = unique(Selected_unique_gene)
  
  Selected_unique_gene_length = length(unique(Selected_A_to_I[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 == Selected_unique_gene[1], ]$Location)) / length(Selected_position)
  for (j in 2:length(Selected_unique_gene)){
    Selected_unique_gene_length = c(Selected_unique_gene_length, 
                                    length(unique(Selected_A_to_I[Selected_A_to_I$Gene.wgEncodeGencodeBasicV34 == Selected_unique_gene[j], ]$Location)) / length(Selected_position))
  }
  
  Final_gene = Selected_unique_gene[Selected_unique_gene_length > 0.1]
  
  Selected_group = rep(i, length(Final_gene))
  Select_df_temp = data.frame(Selected_group,
                         Final_gene)
  colnames(Select_df_temp)=c("Group","Gene")
  Select_df = rbind(Select_df, Select_df_temp)
}

Gene_type = unique(Select_df$Gene)
Gene_type = levels(Gene_type)[Gene_type]
Group_temp = length(Select_df[Select_df$Gene == Gene_type[1],]$Group) / 8
if (Group_temp <= 0.5){
  Last_gene = Gene_type[1]
}

for (i in 2:length(Gene_type)){
  Group_temp = length(Select_df[Select_df$Gene == Gene_type[i],]$Group) / 8
  if (Group_temp <= 0.5){
    Last_gene = c(Last_gene, Gene_type[i])
  }
}

Whole_a_to_I = matrix(0, nrow=length(unique(Last_gene)), ncol=length(rownames(defined_group_8)))
rownames(Whole_a_to_I) = unique(Last_gene)
Group_temp = Select_df[Select_df$Gene == Last_gene[1],]$Group
Group_temp = which(defined_group_8$Group %in% Group_temp)
Whole_a_to_I[1, Group_temp] = 1

for (i in 2:length(Last_gene)){
  Group_temp = Select_df[Select_df$Gene == Last_gene[i],]$Group
  Group_temp = which(defined_group_8$Group %in% Group_temp)
  Whole_a_to_I[i, Group_temp] = 1
}

colnames(Whole_a_to_I) = rownames(defined_group_8)

Target1 = Whole_a_to_I[,defined_group_8$Group == 1]
Target2 = Whole_a_to_I[,defined_group_8$Group == 2]
Target3 = Whole_a_to_I[,defined_group_8$Group == 3]
Target4 = Whole_a_to_I[,defined_group_8$Group == 4]
Target5 = Whole_a_to_I[,defined_group_8$Group == 5]
Target6 = Whole_a_to_I[,defined_group_8$Group == 6]
Target7 = Whole_a_to_I[,defined_group_8$Group == 7]
Target8 = Whole_a_to_I[,defined_group_8$Group == 8]

finalcount = cbind(Target1, Target2, Target3, Target4,
                   Target5, Target6, Target7, Target8)

finalgroup = c(rep("Group1", length(colnames(Target1))),
               rep("Group2", length(colnames(Target2))),
               rep("Group3", length(colnames(Target3))),
               rep("Group4", length(colnames(Target4))),
               rep("Group5", length(colnames(Target5))),
               rep("Group6", length(colnames(Target6))),
               rep("Group7", length(colnames(Target7))),
               rep("Group8", length(colnames(Target8))))

filename='./Sample/A-to-I_heatmap_long_with_border.png'
annotation_df = data.frame(value=c(finalgroup), row.names = colnames(finalcount))
png(filename, width=600, height=6000)
pheatmap(finalcount,cluster_rows = FALSE,cluster_cols=FALSE,show_colnames = FALSE, show_rownames = TRUE,
         border_color="black",
         annotation_col = annotation_df,color=brewer.pal(n=9, name = "Reds"), fontsize_row = 9)
dev.off()