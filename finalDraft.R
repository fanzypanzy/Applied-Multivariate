library(doParallel)
library(itertools)
library(missForest)
library(foreach)
library(qgraph)
library(mlr)
library(psych)
library(GPArotation)
library(HapEstXXR)
library(nFactors)
library(cluster)
library(rgl)
# TODO: Check that we're actually using each library
####################### DATA CLEANING #############################

ddata <- read.csv("road.csv")
#change "number NAs" to NA
ddata[ddata== 99] <- NA
ddata[ddata ==98] <- NA
sum(is.na(ddata))

# removing the rows with more than 18 observations
d <- ddata[-which(rowSums(is.na(ddata)) > 18), ]
which(colSums(is.na(ddata)) == 0)
sum(is.na(d))

#remove ID and group number columns
data <- d[, -c(1,2)]
#get rid of unused levels in factors
data$BL_YADB5 <- droplevels(data$BL_YADB5)



#turn BL_YAD85 into a factor
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
data2 <- data
data2$BL_YADB5 <- as.numeric(levels(data2$BL_YADB5))[data2$BL_YADB5]
data2$BL_YADB5 == data$BL_YADB5
str(data2)

# create a complete dataset to compare imputed dataset to 
comp <- data2[-which(rowSums(is.na(data2)) > 0),]
#imput missing values
set.seed(397)
cl <- makeCluster(2)
registerDoParallel(cl)
im.out.2 <- missForest(xmis = data2, maxiter = 10, ntree = 500,
                       variablewise = FALSE,
                       decreasing = FALSE, verbose = FALSE,
                       replace = TRUE,
                       classwt = NULL, cutoff = NULL, strata = NULL,
                       sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                       parallelize = "variables")
im.out.2$OOBerror
cdata <- cbind(im.out.2$ximp)
stopCluster(cl)
sum(is.na(cdata))

write.csv(cdata, "imputed_dat.csv") 
cdata <- read.csv("imputed_dat.csv")
cdata <- cdata[,-1]



########################### ANALYSIS ################################################

####### A. FACTOR ANALYSIS

#1. a) choose rotation (rotate = ), extraction method (fm = ), and method to calculate factor scores (scores = )

# chose method = "pa". Most common. According to help file, "true" minimal resid is probably found using 
# ran FA with both orthogonal and then a oblique rotation and compared the output results. Oblique roation explained a more
# of the variance than orthogonal. However the difference was very small and since the orthogonal rotation is to be carried 
# over to the cluster analysis, we decided to do a orthogonal rotation. 
efa_var <- fa(cdata[,4:31], nfactors = 3, rotate = "varimax", scores = T, fm = "pa")# factor analysis with n selected factors
efa_pro <- fa(cdata[,4:31], nfactors = 3, rotate = "promax", scores = T, fm = "pa") #



#2. Choose number of factors (-> Liangliang code)
ev <- eigen(cor(cdata[,4:31])) # get eigenvalues
ap <- parallel(subject = nrow(cdata[,4:31]),var = ncol(cdata[,4:31]),
               rep=100,cent=.05)
nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea)
plotnScree(nS)


#3. decide how to deal with complex items
colnum <- c(22, 8, 18 , 11, 9, 25, 13, 23, 5, 14, 7) # vector with all of the complex items from the intial FA with "varimax"
sdata <- cdata[,-c(1,2,3)] # data with just the question columns
# all possible subsets
# don't know which to remove or in what order so want to test with all subsets
subs <- powerset(colnum) # creates a list of all the possible subsets

# does a factor analysis with each possible combination of the complex items removed and stores each resulting RMSE in a vector
for (i in 1:length(subs)) {
  delcol <- subs[[i]]
  splits <- fa(sdata[,-delcol], nfactors = 3, rotate = "varimax", scores = T, fm = "pa")
  rmse[i] <- splits$RMSEA[1] 
  print(i)
}

which.min(rmse) # returns the index for the min RMSE in the vector
our_sub <- as.vector(subs[1272]) # uses the index from above to find the complex item subset 
our_sub # the column numbers of the complex items that if removed, produce the best FA. These need to be removed in this order.
         

#FA without complex items 
efa_splits <- fa(sdata[, -c(11,9,25,13,23,7)], nfactors = 3, rotate = "varimax", scores = T, fm = "pa")
efa_splits

loadings(efa_splits)
factor.plot(efa_splits, labels = rownames(efa_splits$loadings)) #?? should be
fa.diagram(efa_splits) # shows contents of each factor


# 4) Defining factors as indices or scales using Crohnbach's alpha (-> Zhenming)

# Cronbach’s alpha is computed by correlating the score for each scale item with the total score for each observation 
# (usually individual survey respondents or test takers), and then comparing that to the variance for all individual item scores:

# The resulting alpha coefficient of reliability ranges from 0 to 1 in providing this overall assessment of a measure’s reliability. 
# the higher the alpha coefficient, the more the items have shared covariance and probably measure the same underlying concept.

# create keys for three factors
keys.list <- list(Reckless=c("BL_YADB16","BL_YADB17","BL_YADB4","BL_YADB18","BL_YADB14","BL_YADB12","BL_YADB15","BL_YADB6","BL_YADB1",
                            "BL_YADB5", "BL_YADB10"),
                  Distracted=c("BL_YADB27","BL_YADB28","BL_YADB2","BL_YADB8","BL_YADB3","BL_YADB20"),
                  Intoxicated=c("BL_YADB24","BL_YADB22","BL_YADB26","BL_YADB21","BL_YADB19")) 
keys <- make.keys(cdata[,4:31],keys.list)

alpha_scores <- scoreItems(keys,cdata[,4:31])




# 5) Find names for factors


# 6) Make pretty plots


### B. CLUSTER ANALYSIS

# 1) Choose method (2 step, hierarchical, k-means,...)

# Hierarchical clusters
dist.e <- dist(efa_splits$scores, method = 'euclidean')
model1 <- hclust(dist.e, method = 'ward')
result <- cutree(model1, k=3)

# k-means clusters
set.seed(31983)
model2 <- kmeans(efa_splits$scores, centers = 3, nstart = 100)

# 2) Choose number of clusters
wss <- (nrow(efa_splits$scores)-1)*sum(apply(efa_splits$scores,2,var))
perc_explained <- c()
for (i in 2:15){
  mycluster <- kmeans(efa_splits$scores, centers=i)
  wss[i] <- sum(mycluster$withinss)
  perc_explained[i] <- mycluster$betweenss/mycluster$totss
  }
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
plot(1:15, perc_explained, type="b", xlab="Number of Clusters",
     ylab="Percent of Variance explained")

# 3) Bing in external variables (gender/age/experience)
with(model2, table(cluster,cdata$GENDER))
with(model2, table(cluster,cdata$DRIVERLICENSE1))
with(model2, table(cluster,cdata$AGE_LIST))
# 4) Pretty plots

# Not so pretty
clusplot(efa_splits$scores, result, color = TRUE, shade = TRUE, labels = 2, lines = 0)
clusplot(efa_splits$scores, model2$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)
plot3d(efa_splits$scores[,1:3], col=model2$cluster, main="k-means clusters")
plot3d(efa_splits$scores[,1:3], col=result, main="Hierarchical clusters")

pairs(efa_splits$scores[,1:3], col=model2$cluster, labels = c("reckless", "intoxicated", "distracted"),
      pch = cdata$GENDER)




