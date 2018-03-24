terp <- read.csv("~/Desktop/f2_terp.csv")
ag <- read.csv("/Users/admin/Desktop/ag_f2.csv")
cannabinoid <- read.csv("~/Desktop/cannabinoid_f2.csv")

data = read.csv("f2_terp.csv")[,2:8]
X = cor(data,use="p") 

data = read.csv("ag_f2.csv")[,2:9]
X = cor(data,use="p") 

data = read.csv("cannabinoid_f2.csv")[,2:10]
X = cor(data,use="p") 

data = read.csv("alltraits_f2.csv")[,2:25]
X = cor(data,use="p") 

X = cor(data,use="p") 

#pca of correlation of cavariance
#data = read.csv("Vg.csv")[,2:8]
#data = read.csv("Vg_cannabinoid.csv")[,2:10]
data = read.csv("Vg_ag.csv")[,2:9]

X = cor(data,use="p") 

#data

# PCA
pca = prcomp(X)
biplot(pca,col=c(0,1))

# Hierarchical clusters
Dst = dist(t(data))
Clust = hclust(Dst)
plot(Clust)

# Graphical model (lambda 0.45)
#install.packages("huge")
require(huge)
GGM = function(X,L) plot(graph.adjacency(as.matrix(huge(X,lambda=L,method="glasso",verbose=F)$path[[1]]),mode="undirected",diag=F),vertex.label = colnames(X),vertex.size=0)
GGM(X,0.45)


#transpose original dataset
cannabinoidsf2 <- read.csv("~/Desktop/hd_markers_sr.csv")
hdmarker <- t(cannabinoidsf2)
write.csv(hdmarker, file = "hdmarkercannabinoids.csv")


cannabinoid_markers <- read.csv("~/Desktop/hdmarkercannabinoids.csv")

#terp_markers <- read.csv("hdmarker.csv", row.names=1, stringsAsFactors=FALSE)
#terp_markers[terp_markers=="a"]=-1
#terp_markers[terp_markers=="h"]=0
#terp_markers[terp_markers=="b"]=1
#terp_markers=data.matrix(terp_markers)
#write.csv(terp_markers, file = "terp_markers2.csv")

cannabinoid_markers <- read.csv("hdmarkercannabinoids.csv", row.names=1, stringsAsFactors=FALSE)
cannabinoid_markers[cannabinoid_markers=="a"]=-1
cannabinoid_markers[cannabinoid_markers=="h"]=0
cannabinoid_markers[cannabinoid_markers=="b"]=1
cannabinoid_markers=data.matrix(cannabinoid_markers)
write.csv(cannabinoid_markers, file = "cannabinoid_markers.csv")

#install.packages("NAM")
library(NAM)

#IMP(X): Imputes missing points from matrix X with the average value of the column. This function does not return anything, rather it modifies X directly.
IMP(terp_markers)

#GAU(X): Created a Gaussian kernel from matrix X.
G = GAU(terp_markers)
G = GAU(cannabinoid_markers)

terp_pheno <- read.csv("f2_terp.csv", row.names=1, stringsAsFactors=FALSE)
cannabinoid_pheno <- read.csv("cannabinoid_f2.csv", row.names=1, stringsAsFactors=FALSE)
ag_pheno <- read.csv("ag_f2.csv", row.names=1, stringsAsFactors=FALSE)
alltraits_pheno <- read.csv("alltraits_f2.csv", row.names=1, stringsAsFactors=FALSE)


#mean(rownames(terp_markers) == rownames(terp_pheno))
#define y with phenotypic trait data
y = terp_pheno$a.pinene_area
#y = terp_pheno$a_guaiene_area
#y = terp_pheno$s_guaiene_area
#y = terp_pheno$Caryophyllene_area
#y = terp_pheno$Humulene_area
#y = terp_pheno$Limonene_area
#y = terp_pheno$myrcene_area

# Get the kernels
library(NAM)
K = G2A_Kernels(terp_markers)
K = G2A_Kernels(cannabinoid_markers)


# Load the package that runs variance components, and reformat the Kernels from above
library(BGLR)
K2 = lapply(K, function(K) list(K=K,model="RKHS") )

# One trait at a time (in this example, y is the phenotypes)
#FIT1 = BGLR(y, ETA = K2, verbose=FALSE)
#VC = unlist(c(lapply(FIT1$ETA, function(x) x$varU ),Ve=FIT2$varE))
#pie(VC)

# All traits at once (in this example, Y is the matrix of phenotypes)
Y = cannabinoid_pheno
Y = ag_pheno
Y = alltraits_pheno

##go back and add back in logthccbd removed it trying to fix!!!!

VC_ALL = apply(Y, 2, function(y){FIT = BGLR(y,ETA=K2,verbose=FALSE);return(unlist(c(lapply(FIT$ETA, function(x) x$varU ),Ve=FIT$varE)))})
VC_ALL
write.csv(VC_ALL, file = "VC_ALL.csv")
write.csv(VC_ALL, file = "VC_ALL_cannabinoid.csv")
write.csv(VC_ALL, file = "VC_ALL_ag.csv")
write.csv(VC_ALL, file = "VC_ALL_alltraits.csv")


VC_ALL_SCALE = apply(VC_ALL, 2, function(x) x/sum(x) )
VC_ALL_SCALE
write.csv(VC_ALL_SCALE, file = "VC_ALL_SCALE_cannabinoid.csv")
write.csv(VC_ALL_SCALE, file = "VC_ALL_SCALE_ag.csv")


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



# You can also try variance components with this package
#library(varComp)
#FIT_REML = varComp(y~1,varcov=K)


#is.numeric(terp_markers)
#anyNA(terp_markers)


A = K$A
VC_ADD = apply(Y, 2, function(y) reml(y,K=A)$VC)
VC_ADD
write.csv(VC_ADD, file = "VC_ADD.csv")
write.csv(VC_ADD, file = "VC_ADD_cannabinoid.csv")
write.csv(VC_ADD, file = "VC_ADD_ag.csv")


class(Y)

Y = data.matrix(Y)
GenCor = reml(Y,K=A)
GenCor
#save(GenCor, file= "Multivariate_model_results.RData")
save(GenCor, file= "Multivariate_model_results_cannabinoid.RData")
save(GenCor, file= "Multivariate_model_results_ag.RData")

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

mean(is.na(Y))

#making heat maps of correlation matrices
Vg_ag = read.csv("Vg_ag.csv", header=TRUE, row.names=1)
Vg_ag_matrix = as.matrix(Vg_ag)
heatmap(Vg_ag_matrix)

Vg_cannabinoid = read.csv("Vg_cannabinoid.csv", header=TRUE, row.names=1)
Vg_cannabinoid_matrix = as.matrix(Vg_cannabinoid)
heatmap(Vg_cannabinoid_matrix)

Vg_terpene = read.csv("Vg.csv", header=TRUE, row.names=1)
Vg_terpene_matrix = as.matrix(Vg_terpene)
heatmap(Vg_terpene_matrix)



class(Vg_ag_matrix)

Vg_ag

