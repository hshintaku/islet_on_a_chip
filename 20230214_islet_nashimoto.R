library(readxl)
library(Seurat)
library(stringr)
library(stringdist)
library(destiny)

#
# first run
#fstseq <- read.csv("/home/watson/sanger/shintaku/islet/nashimoto_2021.csv",sep = ";",row.names = 1)
fstseq<-Read10X(data.dir = "/home/samba/watson_SeqData/2021Islet_ProfNashimoto/20211124HiSeqX006_Islet/",gene.column = 1)
# detect rat genes
ratgene<-substr(fstseq@Dimnames[[1]],1,4)=="ENSR"
fstseq<-Read10X(data.dir = "/home/samba/watson_SeqData/2021Islet_ProfNashimoto/20211124HiSeqX006_Islet/")
fstmeta<- data.frame(read_xlsx("/home/samba/public/shintaku/islet/cellid_list_culture_1st_seq_230210mod.xlsx"))
rownames(fstmeta)<-fstmeta$cellid
# Create Seurat object with rat genes
fst.seurat <-CreateSeuratObject(fstseq[ratgene,fstmeta$cellid])

#
# second run
#sndseq <-read.csv("/home/watson/sanger/shintaku/islet/nashimoto_2022.csv",sep = ";",row.names = 1)
sndseq<-Read10X(data.dir="/home/samba/watson_SeqData/2021Islet_ProfNashimoto/20220922HiSeqX015_islet/",gene.column = 1)
# detect rat genes
ratgene<-substr(sndseq@Dimnames[[1]],1,4)=="ENSR"
sndseq<-Read10X(data.dir="/home/samba/watson_SeqData/2021Islet_ProfNashimoto/20220922HiSeqX015_islet/")
sndmeta<- data.frame(read_xlsx("/home/samba/public/shintaku/islet/cellid_list_culture_2nd_seq_230210.xlsx"))
#
# convert barcodes to numbers
rownames(sndmeta)<-sndmeta$cellid
cellidlist <- read.table("/home/samba/public/shintaku/github/hunter2/cell_id_list.txt")
source("/home/samba/public/shintaku/github/hunter2/util/whitelist_encode.R")
romin <- whitelist.umi_tools.encode(substr(colnames(sndseq),12,21),cellidlist$V1)
sndseq<-sndseq[,romin$value==0]
romin <- subset(romin ,subset=value==0)
colnames(sndseq)<-paste0("ENSR_",substr(colnames(sndseq),1,11),romin$index)
# Create Seurat object with rat genes
snd.seurat <-CreateSeuratObject(sndseq[ratgene,colnames(sndseq) %in% sndmeta$cellid])
#
# Add meta data to Seurat objects
fst.seurat[["batch"]]<-"first"
fst.seurat[["vascularization"]]<-fstmeta[colnames(fst.seurat),]$vascularization
fst.seurat[["lot"]]<-fstmeta[colnames(fst.seurat),]$lot
#fst.seurat[["species"]]<-fstmeta[colnames(fst.seurat),]$species
snd.seurat[["batch"]]<-"second"
snd.seurat[["vascularization"]]<-sndmeta[colnames(snd.seurat),]$vascularization
snd.seurat[["lot"]]<-sndmeta[colnames(snd.seurat),]$lot
#snd.seurat[["species"]]<-sndmeta[colnames(snd.seurat),]$species
#
# merge first and second runs (remove NA and HUVEC)
#islet <- merge(subset(fst.seurat,subset=vascularization=="Unknown" | vascularization=="minus" | vascularization=="plus"),
#               subset(snd.seurat,subset=vascularization=="Unknown" | vascularization=="minus" | vascularization=="plus"))
# merge first and second runs
islet <- merge(fst.seurat,snd.seurat)
# remove NA and HUVEC
islet <- subset(subset(subset(islet,subset=vascularization=="NA",invert=TRUE),
                subset=vascularization=="HUVEC",invert=TRUE), subset=lot=="Lot1",invert=TRUE)
# remove low quality data
FeatureScatter(islet,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
islet<- subset(islet, subset=nFeature_RNA > 3000)
FeatureScatter(islet,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
#Normalize
islet <- NormalizeData(islet,normalization.method = "LogNormalize")
islet<-FindVariableFeatures(islet, selection.method = "vst",nfeatures=500)
all.genes <- rownames(islet)
islet <- ScaleData(islet, features = all.genes)
islet <-RunPCA(islet,npcs=10, features = VariableFeatures(object = islet))
DimPlot(islet,group.by = "vascularization")+DimPlot(islet,group.by = "batch")

Idents(islet)<-islet[["vascularization"]]
#marker.genes <- FindAllMarkers(islet, only.pos = FALSE, )
marker.genes <- FindMarkers(islet, ident.1 = "plus",ident.2 = "minus")
# check if DEG with significance
View(marker.genes)

# anyways, UMAP
islet <- JackStraw(islet, num.replicate = 100)
islet <- ScoreJackStraw(islet, dims = 1:10)
JackStrawPlot(islet, dims = 1:10)
ElbowPlot(islet)
# perhaps informative up to PC3
islet <- FindNeighbors(islet, dims = 1:3)
islet <- FindClusters(islet)
# let's make UMAP just with PC1 and PC2
islet<-RunUMAP(islet, dims = 1:2)
#islet<-RunTSNE(islet, dims = 1:2)
DimPlot(islet)+DimPlot(islet,group.by = "vascularization")+DimPlot(islet,group.by = "batch")
# extract PC1 & 2 genes
pc1.gene<-PCASigGenes(islet,1,pval.cut = 0.01,use.full = FALSE,max.per.pc = NULL)
pc2.gene<-PCASigGenes(islet,2,pval.cut = 0.01,use.full = FALSE,max.per.pc = NULL)
# combine PC1 & 2 genes
pc.genes <- unique(c(pc1.gene,pc2.gene))
#
# order samples with diffusion map
islet.data<-data.frame(t(data.frame(islet[["RNA"]]@data)))
dm <- DiffusionMap(islet.data[,colnames(islet.data) %in% pc.genes])
# add diffusionmap result to a metadata
islet[["diffusion"]]<-dm$DC1
# check the result in umap
DimPlot(islet,group.by = "vascularization")+FeaturePlot(islet,features = "diffusion")
vascular <-islet[["vascularization"]]
rownames(vascular)<-gsub("-",".",rownames(vascular))
#
# Visualize PC genes with the order defined by the diffusion map
library(pheatmap)
pheatmap(islet.data[order(islet[["diffusion"]]),colnames(islet.data) %in% pc.genes],
         cluster_rows = FALSE,
         annotation_row= vascular)
#
#
# gsea analysis
Idents(islet)<-islet[["vascularization"]]
# now the analysis is done minus against all.
# if you want to compare like minus vs plus
plusfc<-FoldChange(islet,ident.1 = "plus",ident.2 ="minus")
#plusfc<-FoldChange(islet,ident.1 = "plus",ident.2 ="Unknown")
View(plusfc)
FeaturePlot(islet,features = rownames(plusfc[head(order(plusfc$avg_log2FC)),]))
FeaturePlot(islet,features = rownames(plusfc[head(order(plusfc$avg_log2FC,decreasing = TRUE)),]))
#
# GSEA analysis with cluster profiler
library("enrichplot")
library("pathview")
library("org.Rn.eg.db")
library("clusterProfiler")
db<-data.frame(org.Rn.egSYMBOL2EG)
# function for ordering genes
gene_list <- function(plusfc,subdb){
  plusfc_ref <-match(rownames(plusfc),subdb$symbol)
  plusfc$entrez <- subdb[plusfc_ref,]$gene_id
  plusfc <- plusfc[!is.na(plusfc$entrez),]
  gene_list_log2fc<- unlist(plusfc$avg_log2FC)
  names(gene_list_log2fc) <-as.character(plusfc$entrez)
  gene_list_log2fc <- gene_list_log2fc[order(gene_list_log2fc,decreasing = T)]
  return(gene_list_log2fc)
}
# GSEA
gene_list_log2fc <- gene_list(plusfc,db)
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Rn.eg.db,
                   ont          = "BP",
                   minGSSize    = 12,
                   pvalueCutoff = 0.4,
                   pAdjustMethod = "BH",
                   verbose      = FALSE)
# overview
ridgeplot(gse_result,showCategory = 30)
View(gse_result@result)
# extract with "angiogenesis"
angiogenesis.index <-grep("angiogenesis", gse_result$Description)
gse_result$Description[angiogenesis.index]
angiogenes <- strsplit(gse_result@result$core_enrichment[angiogenesis.index],"/")
db[db$gene_id %in% angiogenes[[1]],]$symbol
write.csv2(gse_result@result, "/home/samba/public/shintaku/islet/20230216_gse_result.csv")
# visualize the enrichment
gseaplot(gse_result, geneSetID =angiogenesis.index[1], title = gse_result$Description[angiogenesis.index[1]])
gseaplot(gse_result, geneSetID =angiogenesis.index[2], title = gse_result$Description[angiogenesis.index[2]])
gseaplot(gse_result, geneSetID =angiogenesis.index[3], title = gse_result$Description[angiogenesis.index[3]])
gseaplot(gse_result, geneSetID =angiogenesis.index[4], title = gse_result$Description[angiogenesis.index[4]])
gseaplot(gse_result, geneSetID =angiogenesis.index[5], title = gse_result$Description[angiogenesis.index[5]])
gseaplot(gse_result, geneSetID =angiogenesis.index[6], title = gse_result$Description[angiogenesis.index[6]])
pheatmap(islet.data[order(islet[["diffusion"]]),colnames(islet.data) %in% db[db$gene_id %in% angiogenes[[1]],]$symbol],
         cluster_rows = FALSE,
         annotation_row= vascular)
# extract with "insulin"
insulin.index <-grep("insulin", gse_result$Description)
gse_result$Description[insulin.index]
insulingenes <- strsplit(gse_result@result$core_enrichment[insulin.index],"/")
# visualize the enrichment
gseaplot(gse_result, geneSetID =insulin.index[1], title = gse_result$Description[insulin.index[1]])
gseaplot(gse_result, geneSetID =insulin.index[2], title = gse_result$Description[insulin.index[2]])
gseaplot(gse_result, geneSetID =insulin.index[3], title = gse_result$Description[insulin.index[3]])
gseaplot(gse_result, geneSetID =insulin.index[4], title = gse_result$Description[insulin.index[4]])
gseaplot(gse_result, geneSetID =insulin.index[5], title = gse_result$Description[insulin.index[5]])
gseaplot(gse_result, geneSetID =insulin.index[6], title = gse_result$Description[insulin.index[6]])
pheatmap(islet.data[order(islet[["diffusion"]]),colnames(islet.data) %in% db[db$gene_id %in% insulingenes[[1]],]$symbol],
         cluster_rows = FALSE,
         annotation_row= vascular)
#
#
# GO analysis with PC genes
subdb <- db[db$symbol %in% pc.genes,]
ego_result <- enrichGO(gene          = subdb$gene_id, 
                       OrgDb         = org.Rn.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE, showCategory=30)
View(ego_result@result)
angiogenesis.index <-grep("angiogenesis", ego_result$Description)
ego_result$Description[angiogenesis.index]
ego_result$geneID[angiogenesis.index]
# six genes as examples visualize with FeaturePlot 
FeaturePlot(islet,features = head(unlist(strsplit(ego_result$geneID[1],"/")),6))
#
# genes expressed in naive 
naivefc<-FoldChange(islet,ident.1 = "Unknown")
View(naivefc)
minusfc<-FoldChange(islet,ident.1 = "minus")
#
# optional kegg analysis
#
# kegg_organism = "rat"
# kk <- gseKEGG(geneList     = gene_list_log2fc,
#               organism     = kegg_organism,
#               nPerm=10000,
#               minGSSize    = 3,
#               maxGSSize    = 800,
#               pvalueCutoff = 0.05,
#               pAdjustMethod = "none",
#               keyType       = "ncbi-geneid")
# ridgeplot(kk)
# View(kk@result)
