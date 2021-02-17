library(Seurat)
s1<-readRDS("s1.rds")
s3<-readRDS("s3.rds")
s4<-readRDS("s4.rds")
library(dplyr)
rawdata<-readRDS("cds.experiment.2.counts.rds")
cellinfo<-readRDS("cellinfo.rds")
map<-read.table("Asuu-Cele.map.tsv",header=FALSE,stringsAsFactors=F)
map<-map %>% distinct(V1,.keep_all=TRUE)
map1<-map$V1
map<- as.vector(map$V2)
names(map)=map1
map<-gsub("_","-",map)
map<-map[intersect(names(map),row.names(rawdata))]
rawdata<-rawdata[names(map),]
map<-gsub("^","gene:",map)
map1=names(map)
names(map1)=map
exp2<-CreateSeuratObject(rawdata,project="L2",min.cells = 3, min.features = 50)
subgenes<-intersect(intersect(intersect(row.names(s1),row.names(s3)),row.names(s4)),map)
sub1<-s1[subgenes,]
sub3<-s3[subgenes,]
sub4<-s4[subgenes,]
sub1.count<-as.matrix(sub1[["RNA"]]@counts)
sub3.count<-as.matrix(sub3[["RNA"]]@counts)
sub4.count<-as.matrix(sub4[["RNA"]]@counts)
str(row.names(sub1.count))
row.names(sub1.count)<-map1[row.names(sub1.count)]
row.names(sub3.count)<-map1[row.names(sub3.count)]
row.names(sub4.count)<-map1[row.names(sub4.count)]
str(row.names(sub1.count))
rm(s1,s3,s4)
s1<-CreateSeuratObject(counts = sub1.count,project="S1",min.cells = 3, min.features = 50)
s3<-CreateSeuratObject(counts = sub3.count,project="S3",min.cells = 3, min.features = 50)
s4<-CreateSeuratObject(counts = sub4.count,project="S4",min.cells = 3, min.features = 50)
for(i in c("s1","s3","s4","exp2")){assign(i,SCTransform(get(i),verbose=FALSE))}
int.features<-SelectIntegrationFeatures(object.list = c(s1,s3,s4,exp2), nFeatures=2000)
int.list <- PrepSCTIntegration(object.list = c(s1,s3,s4,exp2), anchor.features = int.features)
library(future)
options(future.globals.maxSize=1200000000)
plan("multiprocess", workers = 8)
int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", anchor.features=int.features, reference=4)
options(future.globals.maxSize=2400000000)
all4<-IntegrateData(anchorset = int.anchors, normalization.method = "SCT",features.to.integrate = map1[subgenes])
all4<-RunPCA(all4)
all4<-RunTSNE(all4,dims=1:30)
all4<-RunUMAP(all4,dims=1:30)
dev<-all4[['pca']]@stdev
library(ggplot2)
p<-ggplot(data = data.frame(dims = 1:30, var=head((dev^2)*100/(sum(dev^2)),n=30))) +
geom_point(mapping = aes_string(x = "dims", y = "var"))
p+theme_classic()
dev.copy(pdf,"all4.elbowplot.variance.pdf")
dev.off()
DimPlot(all4,reduction="pca")
DimPlot(all4,reduction="tsne")
dev.copy(pdf,"all4.tsne.colorbysamples.pdf")
dev.off()
DimPlot(all4,reduction="umap")
dev.copy(pdf,"all4.umap.colorbysamples.pdf")
dev.off()
tail(row.names(all4@meta.data))
t1<-row.names(cellinfo)
str(grep("cele",row.names(all4@meta.data)))
str(row.names(all4@meta.data))
t2<-row.names(all4@meta.data)[15115:22439]
t1<-gsub("$","_4",t1)
sum(t1==t2)
all4@meta.data$cell.type=c(rep("Ascaris",15114),cellinfo$tissue)
str(Idents(all4))
all4.l2<-subset(all4,idents="L2")
DimPlot(all4,reduction="tsne",group.by="cell.type")
dev.copy(pdf,"all4.tsne.colorbycelltype.pdf")
dev.off()
DimPlot(all4,reduction="umap",group.by="cell.type")
dev.copy(pdf,"all4.umap.colorbycelltype.pdf")
dev.off()
DimPlot(all4,reduction="umap",group.by="cell.type",split.by="orig.ident",ncol=2)
dev.copy(pdf,"all4.umap.colorbycelltype.splitbysample.pdf")
dev.off()
DimPlot(all4,reduction="tsne",group.by="cell.type",split.by="orig.ident",ncol=2)
dev.copy(pdf,"all4.tsne.colorbycelltype.splitbysample.pdf")
dev.off()
rm(int.anchors,int.features,int.list)
rm(sub1,sub3,sub4,sub1.count,sub3.count,sub4.count)
rm(s1,s3,sd)
rm(s1,s3,s4)
all4<-FindNeighbors(all4,dims=1:30)
library(future)
options(future.globals.maxSize=2400000000)
plan("multiprocess", workers = 8)
all4<-FindClusters(all4,resolution= c(0.06,0.08,0.1,0.15),save.SNN=TRUE)
levels(all4@meta.data$integrated_snn_res.0.15)=as.character(1:9)
all4@meta.data$seurat_clusters=all4@meta.data$integrated_snn_res.0.15
DimPlot(all4,reduction="tsne",group.by="integrated_snn_res.0.15",label=TRUE,pt.size=0.1,label.size=4,split.by="orig.ident",ncol=2)+NoLegend()
dev.copy(pdf,"all4.tsne.9clusteres0.15.splitbysample.pdf")
dev.off()
library(reshape2)
library(clustree)
clustree(all4,prefix="integrated_snn_res.",node_colour = "sc3_stability")
dev.copy(pdf, "all4.clustree.pdf")
dev.off()
all4@meta.data$seurat_clusters=all4@meta.data$integrated_snn_res.0.15
saveRDS(all4,"all4.rds")
dat13<-melt(table(all4@meta.data$integrated_snn_res.0.15[all4@meta.data$orig.ident=="L2"],all4@meta.data$cell.type[all4@meta.data$orig.ident=="L2"],useNA="always"))
colnames(dat13)<-c("Cluster","Cell_Type","Cell_Count")
library(scales)
cols<-(hue_pal()(8))
names(cols)=unique(dat13$Cell_Type)
cols["Intestine"]="#CCCC00"
barplt13<-ggplot(dat13, aes(x = Cluster, y=Cell_Count, fill=Cell_Type)) +
geom_bar(stat = "identity") +
xlab("\nCluster") +
ylab("Cell Count\n") +
theme_bw()
barplt13<-barplt13+scale_fill_manual("legend", values = cols,na.value="darkgrey")
barplt13<-barplt13+NoLegend()
dat<-melt(table(all4@meta.data$orig.ident,all4@meta.data$integrated_snn_res.0.15))
colnames(dat) <- c("Sample","Cluster","Cell_Count")
barplt<-ggplot(dat, aes(x = Cluster, y=Cell_Count, fill=Sample)) +
geom_bar(stat = "identity") +
xlab("\nCluster") +
ylab("Cell Count\n") +
theme_bw()
fractionbars<-ggplot(dat, aes(x = Cluster, y=Cell_Count, fill=Sample)) +
geom_bar(position="fill",stat = "identity") +
xlab("\nCluster") +
ylab("Fraction cells\n") +
theme_bw()
tsne13<-DimPlot(all4,reduction="tsne",group.by="integrated_snn_res.0.15",label=TRUE,pt.size=0.1,label.size=5)+NoLegend()
barplt<-barplt+NoLegend()
fractionbars13<-fractionbars13+theme_bw()
tsne13|(barplt/fractionbars)|(barplt13/fractionbars13)
ggsave("all4.tsne.9clusteres0.15.barplot.fractionbarplot.bysample.byeleganscelltype.pdf",device="pdf",scale=2,width=12,height=4,units="in")
dev.off()
umap9<-DimPlot(all4,reduction="umap",group.by="integrated_snn_res.0.15",label=TRUE,pt.size=0.1,label.size=5)+NoLegend()
umap9|(barplt/fractionbars)|(barplt13/fractionbars13)
ggsave("all4.umap.9clusteres0.15.barplot.fractionbarplot.bysample.byeleganscelltype.pdf",device="pdf",scale=2,width=12,height=4,units="in")
dev.off()
rm(all4)
q()
