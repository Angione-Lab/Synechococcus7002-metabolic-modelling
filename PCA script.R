setwd("C:/Users/") ##navigate to workspace
##Load packages 
library('devtools')
library('FactoMineR')
library('factoextra')
library('corrplot')
library("PerformanceAnalytics")
ATPflux <- read.csv(file='all_atp_flux.csv',head=FALSE,sep=",")
summary(ATPflux)
res_atp.pca <- PCA(ATPflux)

PCA(ATPflux, scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(res_atp.pca)
print(res_atp.pca)
eigenvalues_atp <- res_atp.pca$eig
variables_atp <- res_atp$var
contributions_atp <-res_atp$var$contrib
plot(res_atp,choix="ind")
plot(res_atp,choix="var")
plot(res_atp,select="contrib 10")
dimdesc(res_atp,axes = 1:2)
atpcontributionsort <- sort(contributions_atp, decreasing = TRUE)
write.csv(contributions_atp, file = "contrib_all_atp_flux.csv")

P1flux <- read.csv(file='all_p1_flux.csv',head=FALSE,sep=",")
summary(P1flux)
res_p1.pca <- PCA(P1flux)

PCA(P1flux, scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(res_p1.pca)
print(res_p1.pca)
eigenvalues_p1 <- res_p1.pca$eig
variables_p1 <- res_p1.pca$var
contributions_p1 <-res_p1.pca$var$contrib
plot(res_p1.pca,choix="ind")
plot(res_p1.pca,choix="var")
plot(res_p1.pca,select="contrib 10")
dimdesc(res_p1.pca,axes = 1:2)
p1contributionsort <- sort(contributions, decreasing = TRUE)
write.csv(contributions_p1, file = "contrib_all_p1_flux.csv")

P2flux <- read.csv(file='all_p2_flux.csv',head=FALSE,sep=",")
summary(P2flux)
res_p2.pca <- PCA(P2flux)

PCA(P2flux, scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(res_p2.pca)
print(res_p2.pca)
eigenvalues_p2 <- res_p2.pca$eig
variables_p2 <- res_p2.pca$var
contributions_p2 <-res_p2.pca$var$contrib
plot(res_p2.pca,choix="ind")
plot(res_p2.pca,choix="var")
plot(res_p2.pca,select="contrib 10")
dimdesc(res_p2.pca,axes = 1:2)
p2contributionsort <- sort(contributions, decreasing = TRUE)
write.csv(contributions_p2, file = "contrib_all_p2_flux.csv")