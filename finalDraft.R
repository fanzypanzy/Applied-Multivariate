library(doParallel)
library(itertools)
library(missForest)
library(foreach)
library(qgraph)
library(mlr)
library(psych)
library(GPArotation)
library(HapEstXXR)
# TODO: Check that we're actually using each library
heyyyy
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
fa.parallel(cdata[,4:31], fa = "fa", n.iter = 100, show.legend = FALSE) # shows number of factors to use



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
factor.plot(efa, labels = rownames(efa_splits$loadings)) #??
fa.diagram(efa_splits) # shows contents of each factor


# 4) Defining factors as indices or scales using Crohnbach's alpha (-> Zhenming)





# 5) Find names for factors


# 6) Make pretty plots


### B. CLUSTER ANALYSIS

# 1) Choose method (2 step, hierarchical, k-means,...)

# 2) Choose number of clusters

# 3) Bing in external variables (gender/age/experience)

# 4) Pretty plots




