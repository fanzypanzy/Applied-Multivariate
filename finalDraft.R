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
library(plot3D)
library(RColorBrewer)
library(GGally)
library(ggplot2)
library(clues)
library(qgraph)
library(caret)
library(NbClust)
library(Amelia)
library(reshape2)
library(corrplot)
library(MVN)
# TODO: Check that we're actually using each library
####################### DATA CLEANING #############################

ddata <- read.csv("road.csv")
#change "number NAs" to NA
ddata[ddata== 99] <- NA
ddata[ddata ==98] <- NA
sum(is.na(ddata))

newnames <- ddata
names(newnames) <- c("ID", "Group", "Gender", "Age", "Experience", 1:28)
missmap(newnames, col = c("#FFFF33", "#377EB8"), legend = FALSE,
        y.labels = "Obs.", y.at = 300, rank.order = FALSE)


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

# turn column GENDER into factor so they are imputed as whole numbers

data2$GENDER <- as.factor(data2$GENDER)


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

# round driver experience and age columns so they are whole numbers
cdata$AGE_LIST <- round(cdata$AGE_LIST)
cdata$DRIVERLICENSE1 <- round(cdata$DRIVERLICENSE1)

# write csv
write.csv(cdata, "imputed_dat.csv") 
cdata <- read.csv("imputed_dat.csv")
cdata <- cdata[,-1] # remove ID row from CSV



########################### ANALYSIS ################################################

####### Exploratory Analysis
newnames2 <- cdata
names(newnames2) <- c("Gender", "Age", "Experience", 1:28)
melted<- melt(newnames2[,4:31])
ggplot(melted,aes(x=variable, y=value, fill=variable)) + geom_boxplot() +
  xlab("Item") + theme(legend.position="none") + ggtitle("Boxplots for Answers to All Items")
  
ggplot(cdata,aes(x=GENDER)) + geom_bar(aes(fill = GENDER), fill=c("#377EB8","#E41A1C")) +
  scale_x_discrete(labels=c("0" = "Male", "1" = "Female")) +
  ggtitle("Gender")

ggplot(cdata,aes(x=as.factor(AGE_LIST))) + geom_bar(fill=c("#4DAF4A")) +
  ggtitle("Age") + xlab("Age (in years)")

ggplot(cdata,aes(x=as.factor(DRIVERLICENSE1))) + geom_bar(fill=c("#FF7F00")) +
  ggtitle("Experience") + xlab("Driving experience (in months)")

corrplot(cor(cdata[,4:31]), tl.col = "black")

# normality check
hzTest(cdata, qqplot = TRUE)

####### PCA and qgraph
pc <- princomp(cdata[,4:31])
summary(pc)
pc$loadings  #eigenvectors
pc$sdev^2   #eigenvalues
pc$scores
screeplot(pc, col ="#4DAF4A", pch = 16, type = "lines", cex = 2, lwd = 2, main = "Screeplot for PCA")
biplot(pc,xlab = "First principal component", ylab = "Second principal component", main = "Biplot",
       col = c("#377EB8", "#E41A1C"), cex = c(0.5,1), ylabs = 1:28, xlabs = 1:756)

qgraph(cor(cdata[,4:31]), labels = 1:28)
qgraph.pca(cor(cdata[,4:31]),rotation = "varimax",factors=3,factorCors=T, labels = 1:28)




####### A. FACTOR ANALYSIS

#1. a) choose rotation (rotate = ), extraction method (fm = ), and method to calculate factor scores (scores = )

# chose method = "pa". Most common. According to help file, "true" minimal resid is probably found using 
# ran FA with both orthogonal and then a oblique rotation and compared the output results. Oblique roation explained a more
# of the variance than orthogonal. However the difference was very small and since the orthogonal rotation is to be carried 
# over to the cluster analysis, we decided to do a orthogonal rotation. 
efa_var <- fa(cdata[,4:31], nfactors = 3, rotate = "varimax", scores = T, fm = "pa")# factor analysis with n selected factors
efa_pro <- fa(cdata[,4:31], nfactors = 3, rotate = "promax", scores = T, fm = "pa") #

print(fa(cdata[,4:31], rotate = "varimax", nfactors = 3))$Vaccounted # 0.49 explained
print(fa(cdata[,4:31], rotate = "promax", nfactors = 3))$Vaccounted # 0.49 explained

#2. Choose number of factors (-> Liangliang code)
ev <- eigen(cor(cdata[,4:31])) # get eigenvalues
ap <- parallel(subject = nrow(cdata[,4:31]),var = ncol(cdata[,4:31]),
               rep=100,cent=.05)
nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea)
plotnScree(nS)

# same thing in pretty
scree <- data.frame(nofactors = 1:28, eigen = nS$Analysis$Eigenvalues, parallel = nS$Analysis$Par.Analysis)
meltedScree<- melt(scree, id.vars = "nofactors")
ggplot(meltedScree, aes(x = nofactors, y =value, col = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#4DAF4A", "#FF7F00"))




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

# make this pretty

itemToFactor <- apply(loadings(efa_splits)[1:22,c(1,3,2)], 1, which.max)
factors <- data.frame(loadings(efa_splits)[1:22,c(1,3,2)], factor = as.factor(itemToFactor), item = 1:22)

p <- ggpairs(data = factors, columns = 1:3,
        mapping= ggplot2::aes(colour = factor),
        upper = list(continuous = "blank"),
        diag = list(continuous = "densityDiag"),
        title = "Loadings of items on factors",
        columnLabels = c("Factor 1","Factor 2", "Factor 3"))
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_color_brewer(palette =  "Set1") +
      scale_fill_brewer(palette =  "Set1")
  }
}
p




# 4) Defining factors as indices or scales using Crohnbach's alpha (-> Zhenming)

# Cronbach’s alpha is computed by correlating the score for each scale item with the total score for each observation 
# (usually individual survey respondents or test takers), and then comparing that to the variance for all individual item scores.

# The resulting alpha coefficient of reliability ranges from 0 to 1 in providing this overall assessment of a measure’s reliability. 
# the higher the alpha coefficient, the more the items have shared covariance and probably measure the same underlying concept.

# create keys for three factors
keys.list <- list(Reckless=c("BL_YADB16","BL_YADB17","BL_YADB4","BL_YADB18","BL_YADB14","BL_YADB12","BL_YADB15","BL_YADB6","BL_YADB1",
                            "BL_YADB5", "BL_YADB10"),
                  Distracted=c("BL_YADB27","BL_YADB28","BL_YADB2","BL_YADB8","BL_YADB3","BL_YADB20"),
                  Intoxicated=c("BL_YADB24","BL_YADB22","BL_YADB26","BL_YADB21","BL_YADB19")) 
keys <- make.keys(cdata[,4:31],keys.list)

alpha_scores <- scoreItems(keys,cdata[,4:31])

cor(efa_splits$scores)

# Interpreting of scoreItems results
# Cronbach’s alpha is always used to test reliability of questionnairs (If different questions have same underlying concept, then answers
# of these questions are tends to be same). In our case, we can use it to check if it is acceptable to group those driving types in a factor.
# Alpha >= 0.9 means excellent, alpha >= 0.8 means good, alpha >= 0.7 means acceptable, alpha < 0.5 means unacceptable.
# Cronbach’s alpha is calculated by estimating the true variance of each item as the average covariance between items while in Guttman 
# lambda 6 the amount of variance in each item are accounted for by the linear regression of all other items.
# Signal/Noise is the ratio of reliable variance to unreliable variance 

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

# TODO: Discuss again (see lecture today)!!! Results are quite different, apparently k-means is better (Janine)


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

# Ch index
# Ch index
ch_k <- c()
ch_hw <- c()

for(i in 1:10){
  mymodel_k <- kmeans(efa_splits$scores, centers = i, nstart = 100)
  model_hw <- hclust(dist.e, method = 'ward.D')
  mymodel_hw <- cutree(hmodel, k=i)
  ch_k[i] <- get_CH(y = efa_splits$scores, mem = mymodel_k$cluster, disMethod = "Euclidean")
  ch_hw[i] <- get_CH(y = efa_splits$scores, mem = mymodel_hw, disMethod = "Euclidean")
}
plot(ch_k, type = "b", col = "red", ylim = c(min(ch_h, na.rm = TRUE), max(ch_k, na.rm = TRUE)))
points(ch_hw, type = "b", col = "blue")

meltedCh<- melt(data.frame(ch_k, ch_hw, noclusters = 1:10), id.vars = "noclusters")

ggplot(meltedCh, aes(x = noclusters, y =value, col = variable)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#984EA3","#FF7F00"),
                     name = "Method", labels = c("k-means", "hierarchical")) +
  labs(title = "CH for various number of clusters and methods", x = "Clusters", y = "CH index") +
  scale_x_discrete(limits = c(2,4,6,8,10))


# suggests 5

gaps <- clusGap(efa_splits$scores,
        FUN = kmeans, 
        K.max = 20, 
        B = 100)
plot(gaps$Tab[,3])

# Duda and Hart's index (only works with hierarchical, so tried with ward.D method)

dh <- NbClust(efa_splits$scores, method = "ward.D")
# suggests 4 (and definitely not 5...)


# Build model with 5 clusters (and others for comparison)
set.seed(3289)
model5 <- kmeans(efa_splits$scores, centers = 5, nstart = 100)
model3 <- kmeans(efa_splits$scores, centers = 3, nstart = 100)
model4 <- kmeans(efa_splits$scores, centers = 4, nstart = 100)
model6 <- kmeans(efa_splits$scores, centers = 6, nstart = 100)


# 3) Bing in external variables (gender/age/experience)
with(model5, table(cluster,cdata$GENDER))
with(model5, table(cluster,cdata$DRIVERLICENSE1))
with(model5, table(cluster,cdata$AGE_LIST))
# 4) Pretty plots

# Not so pretty
clusplot(efa_splits$scores, result, color = TRUE, shade = TRUE, labels = 2, lines = 0)
clusplot(efa_splits$scores, model5$cluster, color = TRUE, shade = TRUE, labels = 2, lines = 0)
plot3d(efa_splits$scores[,1:3], col=model5$cluster, main="k-means clusters")
plot3d(efa_splits$scores[,1:3], col=result, main="Hierarchical clusters")

pairs(efa_splits$scores[,1:3], col=model5$cluster, labels = c("reckless", "intoxicated", "distracted"),
      pch = cdata$GENDER)

# same thing in pretty:

plotclusters3 <- data.frame(PA1 = efa_splits$scores[,1], PA2 = efa_splits$scores[,3],
                            PA3 = efa_splits$scores[,2], cluster = as.factor(model3$cluster))
p3 <- ggpairs(data=plotclusters3, columns = 1:3,
              mapping=ggplot2::aes(colour = cluster),
              columnLabels = c("reckless","distracted", "intoxicated"),
              upper = list(continuous = "density"))
for(i in 1:p3$nrow) {
  for(j in 1:p3$ncol){
    p3[i,j] <- p3[i,j] + 
      scale_fill_manual(values=brewer.pal(3, "Set1")) +
      scale_color_manual(values=brewer.pal(3, "Set1"))  
  }
}

plotclusters4 <- data.frame(PA1 = efa_splits$scores[,1], PA2 = efa_splits$scores[,3],
                            PA3 = efa_splits$scores[,2], cluster = as.factor(model4$cluster))
p4 <- ggpairs(data=plotclusters4, columns = 1:3,
              mapping=ggplot2::aes(colour = cluster),
              columnLabels = c("reckless","distracted", "intoxicated"),
              upper = list(continuous = "density"))
for(i in 1:p4$nrow) {
  for(j in 1:p4$ncol){
    p4[i,j] <- p4[i,j] + 
      scale_fill_manual(values=brewer.pal(4, "Set1")) +
      scale_color_manual(values=brewer.pal(4, "Set1"))  
  }
}


plotclusters5 <- data.frame(PA1 = efa_splits$scores[,1], PA2 = efa_splits$scores[,3],
                           PA3 = efa_splits$scores[,2], cluster = as.factor(model5$cluster))
p5 <- ggpairs(data=plotclusters5, columns = 1:3,
             mapping=ggplot2::aes(colour = cluster),
             columnLabels = c("reckless","distracted", "intoxicated"),
             upper = list(continuous = "density"))
for(i in 1:p5$nrow) {
  for(j in 1:p5$ncol){
    p5[i,j] <- p5[i,j] + 
      scale_fill_manual(values=brewer.pal(5, "Set1")) +
      scale_color_manual(values=brewer.pal(5, "Set1"))  
  }
}

plotclusters6 <- data.frame(PA1 = efa_splits$scores[,1], PA2 = efa_splits$scores[,3],
                            PA3 = efa_splits$scores[,2], cluster = as.factor(model6$cluster))
p6 <- ggpairs(data=plotclusters6, columns = 1:3,
              mapping=ggplot2::aes(colour = cluster),
              columnLabels = c("reckless","distracted", "intoxicated"),
              upper = list(continuous = "density"))
for(i in 1:p6$nrow) {
  for(j in 1:p6$ncol){
    p6[i,j] <- p6[i,j] + 
      scale_fill_manual(values=brewer.pal(6, "Set1")) +
      scale_color_manual(values=brewer.pal(6, "Set1"))  
  }
}




# Pretty 3D plots (not done yet)
scatter3D(efa_splits$scores[,1], efa_splits$scores[,2], efa_splits$scores[,3], bty = "g", pch = 18, 
          col.var = model5$cluster, 
          col = brewer.pal(5, "Set1"),
          pch = 18, ticktype = "detailed",
          colkey = list(side = 2, length = 1, width = 1, at = seq(from = -1.3, to = 7, length.out = 5),
                        labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")),
          theta = 30, phi = 40,
          xlab = "Reckless", ylab = "Intoxicated", zlab = "Distracted",
          alpha = 1, type = "p")



# Inter-item correlation - maybe moce this further up (kind of exploratory)
r <- lowerCor(cdata)
corPlot(r)
findCorrelation(r)







