library(corrplot)
library(car)
library(psych)

### PCA on Genes Variables // Dimension Reduction ###
genesData = METABRIC_RNA_Mutation[c(32:520)]
GenesDataClean <- na.omit(genesData)
pCA = prcomp(GenesDataClean, scale = T)
pA = prcomp(GenesDataClean)
plot(pCA)
plot(pA)
abline(10,0,col = "purple")
pCAA = principal(GenesDataClean, nfactors = 6, rotate = "none")

pAA = principal(GenesDataClean, nfactors = 5, rotate = "varimax")
print.psych(pAA, sort=TRUE)


print.psych(pCAA, sort=TRUE)
options(max.print=5000)
print(pCAA, sort = TRUE, PC4)

summary(pCA)

### ========================================================================

# Updated Dataset with added clinical variables 

data.frame(colnames(METABRIC_RNA_Mutation)) #obtain column # of all variables 

cancerDS = METABRIC_RNA_Mutation[c(24,2,3,5,7,9,12,16,20:22,27,29:30,31,52,57:58,60,64:65,71:72,79,88,145,152,174,180,185,197,210,264,276,280,300,341,366,375,483)]
# cancerDS contains clinical variables + gene variables selected from PCA (prior to transformations)

cancerDSS = METABRIC_RNA_Mutation[c(24,2,3,5,7,9,12,16,20:22,27,29:30,31, 32:520)]
cancerDSS$overall_survival_months = (cancerDSS$overall_survival_months)/12
names(cancerDSS)[1]

cancerGeneDS <- na.omit(cancerDSS)
cancerGeneDS <- subset(cancerGeneDSS, cancer_type_detailed!= "Breast")
cancerGeneDS <- subset(cancerGeneDSS, type_of_breast_surgery!= "") 
cancerGeneDS <- subset(cancerGeneDSS, death_from_cancer!= "" )
allGenes = cancerGeneDS[c(16:504)]

cancerDS <- as.factor(cancerDS$cancer_type_detailed)

#== Variable Transformation:
cancerDS$overall_survival_months = (cancerDS$overall_survival_months)/12 #transform survival_months into years
names(cancerDS)[1] <- "survival_years" #renamed variable to 'survival_years'

#== remove NAs, Missing Information, Blank Entries 
cancerGeneDS <- na.omit(cancerDS) # removes NAs from dataset
cancerGeneDS <- subset(cancerGeneDS, cancer_type_detailed!= "Breast") #remove entries with 'breast' only, missing information
cancerGeneDS <- subset(cancerGeneDS, cancer_type_detailed!= "") #remove entries with blanks
cancerGeneDS <- subset(cancerGeneDS, type_of_breast_surgery!= "") #remove entries with "" only, missing information
cancerGeneDS <- subset(cancerGeneDS, death_from_cancer!= "" ) #remove entries with "" only, missing information

cancerGeneDS$type_of_breast_surgery <- as.factor(cancerGeneDS$type_of_breast_surgery) #transform categorical variable to factors (2 levels)
cancerGeneDS$cancer_type_detailed <- as.factor(cancerGeneDS$cancer_type_detailed) # transform categorical variable to factors (5 levels)
cancerGeneDS$death_from_cancer <- as.factor(cancerGeneDS$death_from_cancer) # transform categorical variable to factors (3 levels)

head(cancerGeneDS) #view variables in cancerGeneDS


clinicalDS = cancerGeneDS[c(1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)] 

clinical2DS = cancerGeneDS[c(1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)] 


library(psych)
pairs.panels(clinical2DS)


#### cancerGeneDS == cleaned and updated dataset with clinical and gene variables.

### ==========================================================================



# PCA / PFA Analysis 

genesDS = cancerGeneDS[c(16:40)]
head(genesDS)


corrplot(cor(genesDS), order="AOE", method = "ellipse")
corrplot(cor(genesDS), order="AOE", method = "shade")

#initial PCA:

p = prcomp(genesDS, scale = T) #scaled PCA
summary(p)
plot(p)
abline(1, 0, col = "purple") # Knee and Var=1, shows 4 Factors

print(p)
round(p$rotation, 3)


P = prcomp(genesDS) #unscaled PCA
plot(P) # same as scaled PCA


# Perform Parallel Analysis 

nrow(genesDS)
ncol(genesDS)
m = matrix(rnorm(25*1291, 0, 1), ncol=25)
gDS = data.frame(m)
corrplot(cor(gDS))


pR = prcomp(gDS, scale=T)
plot(pR)
abline(1, 0, col = "green")

library(psych)
pFA = fa.parallel(gDS, n.iter=1000)


# PFA w/ VARIMAX 

pP = principal(genesDS, nfactors = 4, rot = 'varimax') # 4 Factors
print(pP$loadings, cutoff = .4, sort=T)


pQ = principal(genesDS, nfactors = 5, rot = 'varimax')
print(pQ$loadings, cutoff = .4)

pR = principal(genesDS, nfactors = 3, rot = 'varimax')
print(pR$loadings, cutoff = .4, sort=T)


corrplot(cor(genesDS), order = "hclust")

#redo without ctcf variable:

genes2DS = genesDS[-c(23)]
corrplot(cor(genes2DS), order = "hclust")
corrplot(cor(genes2DS), order = "AOE", method = "shade")

PP = principal(genes2DS, nfactors = 3, rot = 'varimax')
print(PP$loadings, cutoff = .4, sort=T)

CF = factanal(genes2DS, 3)
print(CF$loadings, cutoff = .4, sort=T)


library(stats)
#Compare to Common Factor Analysis

cFA = factanal(genesDS, 4)
print(cFA$loadings, cutoff = .4, sort=T)

cFA3 = factanal(genesDS, 3)
print(cFA3$loadings, cutoff = .4, sort=T)


# Goodness of Fit:
KMO(genesDS)


# Factors for Regression:

FacReg = cancerGeneDS[c(1, 16, 22, 26)]
head(FacReg)
fit = lm(survival_years ~., data = FacReg)
summary(fit)


# Multi-Dimensional Scaling (MDS) Using Reduced Dimensionality ///

install.packages("MASS")
library(MASS)

#MDS - Checking Stress
gene.dist = dist(genesDS)
gene.mds = isoMDS(gene.dist)
gene.mds$stress/100


#Plot MDS to Check Clusters 
plot(gene.mds$points)

#Calculate Shepard's Diagram 
gene.sh = Shepard(gene.dist, gene.mds$points)
plot(gene.sh, pch=".")
#lines(gene.sh$x, gene.sh$y, type="S", col="pink")


d = dist(genesDS)
fit = cmdscale(d, eig=TRUE, k=2)
fit

#plot fit
x = fit$points[,1]
y = fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS")
text(x, y, labels = row.names(genesDS), cex=.7)


# Density Cluster
install.packages(c("dbscan", "fpc"))
library(dbscan)

dens = dbscan(genesDS, eps = 1) # no clusters
dens

dens = dbscan(genesDS, eps = 0.5, minPts = 2) 
dens


#Hierarchical Clustering 
clusterH = hclust(d)
plot(clusterH, cex = .6, hang = -1)
clusters = cutree(clusterH, k = 2) 
plot(fit$points, col=clusters)
plot(clusters)
rect.hclust(clusterH, k = 2, border = 2:5)

clusters = cutree(clusterH, k = 10) # k = 2
plot(fit$points, col=clusters)



clusters = cutree(clusterH, k = 5) # k = 
plot(fit$points, col=clusters)


#K-Means Clustering 
geneClust = kmeans(genesDS, 3)



install.packages("mlbench")
install.packages("kernlab")
library(mlbench)
library(kernlab)

install.packages("factoextra")
library(factoextra)

fviz_nbclust(genesDS, pam, method ="silhouette")
fviz_nbclust(genesDS)

install.packages("cluster")

library(cluster)

fviz_nbclust(genesDS, pam, method ="silhouette")
fviz_nbclust(genesDS, hcut, method ="gap")

genePam = pam(genesDS, k = 2)
fviz_cluster(genePam, data = genesDS)
fviz_cluster(genePam, data = genesDS, ellipse.type = "euclid")

kmeansG = kmeans(genesDS, 3)
fviz_cluster(kmeansG, data = genesDS, ellipse.type = "euclid")

geneCut = hcut(genesDS, k = 3, stand=TRUE)
fviz_cluster(geneCut, data = genesDS, ellipse.type = "euclid")

fviz_cluster(geneCut, data = genesDS)

clara = clara(genesDS, 2)
fviz_cluster(clara, data = genesDS)
