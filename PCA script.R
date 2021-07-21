#Navigate to workspace
setwd("C:/Users/") 

#Load packages 
library('devtools')
library('FactoMineR')
library('factoextra')
library('corrplot')
library("PerformanceAnalytics")

#Load transcript/multiomic/flux .csv data files manually (transcriptsnew/all_ATPTF/all_P1TF/all_P2TF/all_atp_flux/all_p1_flux/all_p2_flux)
transcripts <- read.csv(file='transcriptsnew.csv',head=FALSE,sep=",")
ATPTF <- read.csv(file='all_ATPTF.csv',head=FALSE,sep=",")
P1TF <- read.csv(file='all_P1TF.csv',head=FALSE,sep=",")
P2TF <- read.csv(file='all_P2TF.csv',head=FALSE,sep=",")
ATPflux <- read.csv(file='all_atp_flux.csv',head=FALSE,sep=",")
P1flux <- read.csv(file='all_p1_flux.csv',head=FALSE,sep=",")
P2flux <- read.csv(file='all_p2_flux.csv',head=FALSE,sep=",")

#Perform PCA for each dataset 
res_transcripts.pca <- PCA(transcripts)
res_ATPTF.pca <- PCA(ATPTF)
res_P1TF.pca <- PCA(P1TF
res_P2TF.pca <- PCA(P2TF)
res_ATPflux.pca <- PCA(ATPflux)
res_P1flux.pca <- PCA(P1flux)
res_P2flux.pca <- PCA(P2flux)

#Plot PCA graph of individuals colored according to cos2 values for each dataset
transcripts_PCA_plot <- fviz_pca_ind(res_transcripts.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

ATPTF_PCA_plot <- fviz_pca_ind(res_ATPTF.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

P1TF_PCA_plot <- fviz_pca_ind(res_P1TF.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

P2TF_PCA_plot <- fviz_pca_ind(res_P2TF.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

ATPflux_PCA_plot <- fviz_pca_ind(res_ATPflux.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

P1flux_PCA_plot <- fviz_pca_ind(res_P1flux.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

P2flux_PCA_plot <- fviz_pca_ind(res_P2flux.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
             )

#Save PCA plots to .pdf files in the workspace
pdf("transcripts_PCA.pdf")
print(transcripts_PCA_plot)
dev.off()
pdf("ATPTF_PCA.pdf")
print(ATPTF_PCA_plot)
dev.off()
pdf("P1TF_PCA.pdf")
print(P1TF_PCA_plot)
dev.off()
pdf("P2TF_PCA.pdf")
print(P2TF_PCA_plot)
dev.off()
pdf("ATPflux_PCA.pdf")
print(ATPflux_PCA_plot)
dev.off()
pdf("P1flux_PCA.pdf")
print(P1flux_PCA_plot)
dev.off()
pdf("P2flux_PCA.pdf")
print(P2flux_PCA_plot)
dev.off()

#Obtain contributions of principal component variables (genes/metabolic reactions/both) for each dataset
contributions_transcripts <-res_transcripts.pca$var$contrib
contributions_ATPTF <-res_ATPTF.pca$var$contrib
contributions_P1TF <-res_P1TF.pca$var$contrib
contributions_P2TF <-res_P2TF.pca$var$contrib
contributions_ATPflux <-res_ATPflux.pca$var$contrib
contributions_P1flux <-res_P1flux.pca$var$contrib
contributions_P2flux <-res_P2flux.pca$var$contrib

#Save contributions to new .csv files in the workspace
write.csv(contributions_transcripts, file = "contrib_transcripts.csv")
write.csv(contributions_ATPTF, file = "contrib_ATPTF.csv")
write.csv(contributions_P1TF, file = "contrib_P1TF.csv")
write.csv(contributions_P2TF, file = "contrib_P2TF.csv")
write.csv(contributions_ATPflux, file = "contrib_ATPflux.csv")
write.csv(contributions_P1flux, file = "contrib_P1flux.csv")
write.csv(contributions_P2flux, file = "contrib_P2flux.csv")

#Obtain principal component coordinates for individual growth conditions
ind_coord_transcripts <-res_transcripts.pca$ind$coord
ind_coord_ATPTF <-res_ATPTF.pca$ind$coord
ind_coord_P1TF <-res_P1TF.pca$ind$coord
ind_coord_P2TF <-res_P2TF.pca$ind$coord
ind_coord_ATPflux <-res_ATPflux.pca$ind$coord
ind_coord_P1flux <-res_P1flux.pca$ind$coord
ind_coord_P2flux <-res_P2flux.pca$ind$coord

#Save coordinates to new .csv files in the workspace
write.csv(ind_coord_transcripts, file = "ind_coord_transcripts.csv")
write.csv(ind_coord_ATPTF, file = "ind_coord_ATPTF.csv")
write.csv(ind_coord_P1TF, file = "ind_coord_P1TF.csv")
write.csv(ind_coord_P2TF, file = "ind_coord_P2TF.csv")
write.csv(ind_coord_ATPflux, file = "ind_coord_ATPflux.csv")
write.csv(ind_coord_P1flux, file = "ind_coord_P1flux.csv")
write.csv(ind_coord_P2flux, file = "ind_coord_P2flux.csv")

