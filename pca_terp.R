#defining data sets
terp <- read.csv("~/Desktop/f2_terp.csv")
data = read.csv("f2_terp.csv")[,2:8]
X = cor(data,use="p") 

ag <- read.csv("/Users/admin/Desktop/ag_f2.csv")
data = read.csv("ag_f2.csv")[,2:9]
X = cor(data,use="p") 

cannabinoid <- read.csv("~/Desktop/cannabinoid_f2.csv")
data = read.csv("cannabinoid_f2.csv")[,2:10]
X = cor(data,use="p") 

alltraits <- read.csv("~/Desktop/alltraits_f2.csv")
data = read.csv("alltraits_f2.csv")[,2:25]
X = cor(data,use="p") 

#Unsupervised Machine Learning Plots: PCA, Hierarchical Clusters, and Graphical Model
# PCA
pca = prcomp(X)
biplot(pca,col=c(0,1))

#Hierarchical clusters
Dst = dist(t(data))
Clust = hclust(Dst)
plot(Clust)

#Graphical model (lambda 0.45) 
#install.packages("huge")
require(huge)
GGM = function(X,L) plot(graph.adjacency(as.matrix(huge(X,lambda=L,method="glasso",verbose=F)$path[[1]]),mode="undirected",diag=F),vertex.label = colnames(X),vertex.size=0)
GGM(X,0.45)

#pca of correlation of cavariance
#data = read.csv("Vg.csv")[,2:8]
X = cor(data,use="p")

#data = read.csv("Vg_cannabinoid.csv")[,2:10]
X = cor(data,use="p")

data = read.csv("Vg_ag.csv")[,2:9]
X = cor(data,use="p")

#Transpose original dataset and saving transposed dataset (only necessary if needed to have both datasets in equivalent orientation)
cannabinoidsf2 <- read.csv("~/Desktop/hd_markers_sr.csv")
hdmarker <- t(cannabinoidsf2)
write.csv(hdmarker, file = "hdmarkercannabinoids.csv")

#Reading in transposed dataset
cannabinoid_markers <- read.csv("~/Desktop/hdmarkercannabinoids.csv")

#Defining markers and saving file (should rename to marker data matrix as all phenotypic data from same mapping population)
terp_markers <- read.csv("hdmarker.csv", row.names=1, stringsAsFactors=FALSE)
terp_markers[terp_markers=="a"]=-1
terp_markers[terp_markers=="h"]=0
terp_markers[terp_markers=="b"]=1
terp_markers=data.matrix(terp_markers)
write.csv(terp_markers, file = "terp_markers2.csv")

cannabinoid_markers <- read.csv("hdmarkercannabinoids.csv", row.names=1, stringsAsFactors=FALSE)
cannabinoid_markers[cannabinoid_markers=="a"]=-1
cannabinoid_markers[cannabinoid_markers=="h"]=0
cannabinoid_markers[cannabinoid_markers=="b"]=1
cannabinoid_markers=data.matrix(cannabinoid_markers)
write.csv(cannabinoid_markers, file = "cannabinoid_markers.csv")

#Install NAM for analyses below
#install.packages("NAM")
library(NAM)

#IMP(X): Imputes missing points from matrix X with the average value of the column. This function does not return anything, rather it modifies X directly. Incase of missing data. 
IMP(terp_markers)

#GAU(X): Created a Gaussian kernel from matrix X.
G = GAU(terp_markers)
G = GAU(cannabinoid_markers)

#Define phenotypic datasets
terp_pheno <- read.csv("f2_terp.csv", row.names=1, stringsAsFactors=FALSE)
cannabinoid_pheno <- read.csv("cannabinoid_f2.csv", row.names=1, stringsAsFactors=FALSE)
ag_pheno <- read.csv("ag_f2.csv", row.names=1, stringsAsFactors=FALSE)
alltraits_pheno <- read.csv("alltraits_f2.csv", row.names=1, stringsAsFactors=FALSE)

#Dataset check should equal one
mean(rownames(terp_markers) == rownames(terp_pheno))
mean(rownames(cannabinoid_markers) == rownames(cannabinoid_pheno))
mean(rownames(ag_markers) == rownames(ag_pheno))
mean(rownames(alltraits_markers) == rownames(alltraits_pheno))


#Define y with phenotypic trait data
y = terp_pheno$a.pinene_area
#y = terp_pheno$a_guaiene_area
#y = terp_pheno$s_guaiene_area
#y = terp_pheno$Caryophyllene_area
#y = terp_pheno$Humulene_area
#y = terp_pheno$Limonene_area
#y = terp_pheno$myrcene_area

#Define K as a Guassian Kernel using market dataset. Computes a list of orthogonal kernels containing additive, dominant and first-order epistatic effects. 
library(NAM)
K = G2A_Kernels(terp_markers)
K = G2A_Kernels(cannabinoid_markers)

#Load the package that runs variance components, and reformat the Guassian Kernels from above
library(BGLR)
K2 = lapply(K, function(K) list(K=K,model="RKHS") )

#One trait at a time is currently broken
#One trait at a time (in this example, y is the phenotypes)
#FIT1 = BGLR(y, ETA = K2, verbose=FALSE)
#VC = unlist(c(lapply(FIT1$ETA, function(x) x$varU ),Ve=FIT2$varE))
#pie(VC)

#All traits at once (in this example, Y is the matrix of phenotypes)
Y = terp_pheno
Y = cannabinoid_pheno
Y = ag_pheno
Y = alltraits_pheno

##go back and add back in logthccbd removed it from all traits during code testing trying to fix!!!!

#Annotate
VC_ALL = apply(Y, 2, function(y){FIT = BGLR(y,ETA=K2,verbose=FALSE);return(unlist(c(lapply(FIT$ETA, function(x) x$varU ),Ve=FIT$varE)))})
VC_ALL
write.csv(VC_ALL, file = "VC_ALL.csv")
write.csv(VC_ALL, file = "VC_ALL_cannabinoid.csv")
write.csv(VC_ALL, file = "VC_ALL_ag.csv")
write.csv(VC_ALL, file = "VC_ALL_alltraits.csv")

#Annotate
VC_ALL_SCALE = apply(VC_ALL, 2, function(x) x/sum(x) )
VC_ALL_SCALE
write.csv(VC_ALL_SCALE, file = "VC_ALL_SCALE_cannabinoid.csv")
write.csv(VC_ALL_SCALE, file = "VC_ALL_SCALE_ag.csv")

#plot a pie chart of heritability
pie( VC_ALL[, "a.guaiene_area"]  )
pie( VC_ALL[, "s.guaiene_area"]  )
pie( VC_ALL[, "Caryophyllene_area"]  )
pie( VC_ALL[, "Limonene_area"]  )
pie( VC_ALL[, "Humulene_area"]  )
pie( VC_ALL[, "myrcene_area"]  )

pie( VC_ALL[, "THCV"]  )
pie( VC_ALL[, "THC"]  )
pie( VC_ALL[, "CBD"]  )
pie( VC_ALL[, "CBC"]  )
pie( VC_ALL[, "CBG"]  )
pie( VC_ALL[, "CBN"]  )
pie( VC_ALL[, "Total.cannabinoids"]  )
pie( VC_ALL[, "Log.THC.CBD"]  )
pie( VC_ALL[, "THC.CBD"]  )

#Broken unable to download on my version of R
#You can also try variance components with this package
#library(varComp)
#FIT_REML = varComp(y~1,varcov=K)

#Dataset checks; need a numeric dataset without missing data
is.numeric(terp_markers)
anyNA(terp_markers)

#Annotate making and saving
A = K$A
VC_ADD = apply(Y, 2, function(y) reml(y,K=A)$VC)
VC_ADD
write.csv(VC_ADD, file = "VC_ADD.csv")
write.csv(VC_ADD, file = "VC_ADD_cannabinoid.csv")
write.csv(VC_ADD, file = "VC_ADD_ag.csv")

#dataset check data should be numeric
class(Y)

#Making and saving a correlation table 
Y = data.matrix(Y)
GenCor = reml(Y,K=A)
GenCor
save(GenCor, file= "Multivariate_model_results.RData")
save(GenCor, file= "Multivariate_model_results_cannabinoid.RData")
save(GenCor, file= "Multivariate_model_results_ag.RData")

#Saving heritability files as csv format
#write.csv(GenCor$VC$Vg, "Additive_Genetic_Covariance.csv")
#write.csv(GenCor$VC$Ve, "Residual_Covariance.csv")
#write.csv(GenCor$VC$h2, "h2.csv")

write.csv(GenCor$VC$Vg, "Additive_Genetic_Covariance_cannabinoid.csv")
write.csv(GenCor$VC$Ve, "Residual_Covariance_cannabinoid.csv")
write.csv(GenCor$VC$h2, "h2_cannabinoid.csv")

write.csv(GenCor$VC$Vg, "Additive_Genetic_Covariance_ag.csv")
write.csv(GenCor$VC$Ve, "Residual_Covariance_ag.csv")
write.csv(GenCor$VC$h2, "h2_ag.csv")

Vg <- cov2cor(GenCor$VC$Vg)
Vg
write.csv(Vg, "Vg_cannabinoid.csv")
write.csv(Vg, "Vg_ag.csv")

Vg <- read.csv("Additive_Genetic_Covariance.csv", row.names=1, stringsAsFactors=FALSE)

#Dataset check should be zero
mean(is.na(Y))

#Making heat maps of correlation matrices
Vg_ag = read.csv("Vg_ag.csv", header=TRUE, row.names=1)
Vg_ag_matrix = as.matrix(Vg_ag)
heatmap(Vg_ag_matrix)

Vg_cannabinoid = read.csv("Vg_cannabinoid.csv", header=TRUE, row.names=1)
Vg_cannabinoid_matrix = as.matrix(Vg_cannabinoid)
heatmap(Vg_cannabinoid_matrix)

Vg_terpene = read.csv("Vg.csv", header=TRUE, row.names=1)
Vg_terpene_matrix = as.matrix(Vg_terpene)
heatmap(Vg_terpene_matrix)

#Dataset check class should be matrix
class(Vg_ag_matrix)

Vg_ag

