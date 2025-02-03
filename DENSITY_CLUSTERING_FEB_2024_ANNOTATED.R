#  A SCRIPT TO ANALYSE REAL AND ARTIFICIAL CLUSTERING DATA - LATTER USUALLY HAS PRE-ASSIGNED LABELS FOR TESTING
# THE ANALYSIS IS BASED ON A MODIFIED VERSION OF DENSITYCUT ; DING, J., SHAH, S  AND CONDON, A. (2016)   AN EFFICIENT AND VERSATILE TOPOLOGICAL APPROACH TO AUTOMATIC CLUSTERING OF BIOLOGICAL DATA.BIOINFORMATICS, 32, 2567-2576.
# THE MAIN MODIFICATION IS TO REPLACE A SINGLE METHOD OF IDENTIFYING LOCAL DENSITY BASED ON K NEAREST NEIGHBOURS WITH 3 METHODS, BASED ON KNN, SNN(SHARED NEAREST NEIGHBOURS) AND RKNN ( REVERSE K NEAREST NEIGHBOURS)
# THE METHOD FINDS THE OPTIMAL CLUSTER NUMBER AS THAT WITH MAXIMUM FREQUENCY OF OCCURRENCE OVER MEULTIPLE CLUSTERINGS FOUND USING ALL COMBINATIONS OF A RANGE OF ALPHA PARAMETERS FROM DENSITYCUT AND THE THREE METHODS OF DENSITY ESTIMATION
# THE METHOD ALSO COMPUTES THE CLUSTER NUMBERS WITH THE 2ND AND 3RD HIGHEST FREQUENCIES OF OCCURRENCE
# OUPUTS ARE
# 1) - TABLE OF CLUSTER NUMBERS AND FREQUENCIES SHOWING HOW COMMON EACH CLUSTER NUMBER IS
# 2) - A GRAPHIC PDF FILE: ALLCLUSTERSBYFREQUENCY.PDF WHICH SHOWS 1) A PLOT OF FREQUENCIES OF CLUSTER NUMBERS
#                                                           2) CLUSTER MEMBERSHIPS OF EACH OBSERVATION FOR HIGHEST FREQUENCY CLUSTERING
#                                                          3) AS 2) BUT WITH OUTLIERS IN BLACK IF DENOISE ==1
#                                                          4) DISTRIBUTION OF DIFFERENT ALTERNATIVE CONFIGURATIONS WITHIN A GIVEN CLUSTER NUMBER
#                                                          ABOVE REPEATED FOR 2ND AND 3RD MOST COMMON CLUSTER NUMBERS
# 3) CLUSTERING V MODAL CLUSTERING CLUSTER NUMBER FOR EACH SUBJECT. CLUSTERING IS ALLOCATION FOR MOST FREQUENT CONFIGURATION OF MOST COMMON CLUSTER NUMBER. MODAL
#                                                           CLUSTERING IS MOST COMMON CLUSTER ALLOCATION ACROSS ALL CLUSTER CONFIGURATIONS. TABLE IS PRODUCED FOR 3 MOST COMMON CLUSTER NUMBERS.
# 4) CLUSTERING MODESUMS TABLE GIVES PROBABILITY FOR EACH SUBJECT OF BELONGING TO EACH CLUSTER.  REPEATED FOR EACH OF 3 BEST CLUSTER NUMBERS.

library(clusterlab)
library(cluster)
library(pracma)
library(dbscan)
library(densitycut)
library(Spectrum)
library(clusteringdatasets)
library(randomForest)
library(aricode)
library(scales)
library(synthpop)
library(factoextra)     
library(FCPS)
library("doParallel")
library(proxy)
library(GenForImp)
library(reticulate)
library(ggplot2)
library(future.apply)
registerDoParallel(5)
eif<-import("eif")
library(missForest)
library(fossil)
library(dplyr)
library(emstreeR)
library(DescTools)
source("clfunc10c.R")

find_knn_alpha<-function(x){
s<-x%%3
ss<-s
if(s==0){s<-3}
a <- x%/%3
if(ss==0){a<-a-1}
a<-a+1
return(list(metallnames[s],alphas[a]))
}
rowsInTbl <- function(tbl,row){
  sum(apply(tbl, 1, function(x) all(x == row) ))
}
colFrequency <- function(tblall){
  tbl <- unique(tblall)
  results <- matrix(nrow = nrow(tbl),ncol=ncol(tbl)+1)
  results[,1:ncol(tbl)] <- as.matrix(tbl)
  freq <- apply(tbl,1,function(x)rowsInTbl(tblall,x))
  results[,ncol(tbl)+1] <- freq
  return(results)
}

Modefunc3<-function(x,y){
  if(is.matrix(x)){
    dimx<-dim(x)
    lx<-dimx[2]
  }else{
    lx<-length(x)
  }
  outclust<-rep(0,lx)
  if(is.matrix(x)){
  allclustmembs<-matrix(0,lx,length(unique(x[1,])))
  }else{
    allclustmembs<-matrix(0,lx,length(unique(x))) 
  }
  modelengths<-rep(0,lx)
  totfr<-sum(y)
  
  for(i in 1:lx){
    if(is.matrix(x)){
      newmodevec<-rep(x[,i],times=y)
    }else{
      newmodevec<-rep(x[i],times=y) 
    }
    
    z<-Mode(newmodevec)
    if(length(z)>1){
      z<-z[1]
    }
    if(is.matrix(x)){
    clvals<-rep(0,length(unique(x[1,])))
    }else{
    clvals<-rep(0,length(unique(x)))
    }
    mm<-as.matrix(table(newmodevec))
    nn<-as.numeric(rownames(mm))
    index<-1
    for(j in nn){
      clvals[j]<-mm[index,1]
      index<-index+1
    }
    
    allclustmembs[i,]<-clvals
    outclust[i]<-which(clvals==max(clvals))
    modelengths[i]<-length(which(newmodevec==z))
  }
  return(list(modelengths,allclustmembs/sum(freqcss),outclust))
}
Modefunc3_par<-function(x,y){
  if(is.matrix(x)){
    dimx<-dim(x)
    lx<-dimx[2]
  }else{
    lx<-length(x)
  }
  outclust<-rep(0,lx)
  if(is.matrix(x)){
    allclustmembs<-matrix(0,lx,length(unique(x[1,])))
  }else{
    allclustmembs<-matrix(0,lx,length(unique(x))) 
  }
  modelengths<-rep(0,lx)
  totfr<-sum(y)
  
  testaaa<-foreach(i = 1:lx)%dopar%{
    if(is.matrix(x)){
      newmodevec<-rep(x[,i],times=y)
    }else{
      newmodevec<-rep(x[i],times=y) 
    }
    
    z<-Mode(newmodevec)
    if(length(z)>1){
      z<-z[1]
    }
    if(is.matrix(x)){
      clvals<-rep(0,length(unique(x[1,])))
    }else{
      clvals<-rep(0,length(unique(x)))
    }
    mm<-as.matrix(table(newmodevec))
    nn<-as.numeric(rownames(mm))
    index<-1
    for(j in nn){
      clvals[j]<-mm[index,1]
      index<-index+1
    }
    
    list(clvals,z,newmodevec)
    
  }
  for(i in 1:lx){
    outclust[i]<-which(testaaa[[i]][[1]]==max(testaaa[[i]][[1]])) 
    modelengths[i]<-length(which(testaaa[[i]][[3]]==testaaa[[i]][[2]]))
    allclustmembs[i,]<-testaaa[[i]][[1]]
  }
  
  
  return(list(modelengths,allclustmembs/sum(freqcss),outclust))
}
count.duplicates2 <- function(DF){
  df<-as.data.frame(DF)
  return(aggregate(list(numdup=rep(1,nrow(df))), df, length))
  }
changeclustnums<-function(x){
  ux<-unique(x)
  uxx<-length(unique(x))
  newclusts<-rep(0,length(x))
  for(i in 1:uxx){
    wx<-which(x==ux[i])
   newclusts[wx]<-i 
  }
  return(newclusts)
}
changereordclustnums<-function(x){
  ux<-unique(x)
  uxs<-sort(ux)
  uxx<-length(unique(x))
  newclusts<-rep(0,length(x))
  for(i in 1:uxx){
    wx<-which(x==uxs[i])
    newclusts[wx]<-i 
  }
  return(newclusts)
}


findsupermatchnew2<-function(m,x,probe){
  ll<-list()
  listindex<-1
  for(i in 1:x){
    
        if(m[i]==probe){
          ll[[listindex]]<-i
          listindex<-listindex+1
        
    }
  } 
  return(ll)
}

colMins<-function(X){
  dx<-dim(X)
  mins<-rep(0,dx[1])
  for(i in 1:length(mins)){
    mins[i]<-X[i,which(X[i,]==min(X[i,]))]
  }
  return(mins)
}
clustdistprob<-function(prtot,clx,startmodevals,modevals,num){
  startclustids<-rep(0,length(startmodevals))
  for(i in 1:length(startmodevals)){
    for(j in unique(clx)){
      xx<-match(startmodevals[i],which(clx==j))
      if(! is.na(xx)){
        startclustids[i]<-j
      }
    }
  }
  pp<-clprob4(mb1$samples,clx,startmodevals,startclustids)
  relpvals<-pp
  relpvals[which(relpvals==0)]<-1
  relpvals<-apply(relpvals, 1, function(x) {x/sum(x)})
  relpvals<-t(relpvals)
  write(t(relpvals),file=paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"best_",num,"_relative_p_clustmaximp_vals",sep=""),ncolumns=length(unique(clx)))
  for(i in unique(clx)){
    relpvals[clx==i,]<-relpvals[clx==i,]/max(relpvals[clx==i,])
  }
  allcols<-rainbow(max(clx))
  allcolsorig<-allcols
  absmax<-min(6*max(clx),75)
  cexs<-rep(1,dimd[1])
  pchs<-rep(1,dimd[1])
  badpts<-rep(0,dimd[1])
  maxcols<-rep(1,dimd[1])
  for(i in 1:dimd[1]){
    maxcol<-which(relpvals[i,]==max(relpvals[i,]))
    if(length(maxcol)>1){
      maxcol<-maxcol[1]
      badpts[i]<-1
    }
    print(maxcol)
    pchs[i]<-16
    cexs[i]<-2*relpvals[i,maxcol]
    maxcols[i]<-maxcol
  }
  wps<-which(badpts==1)
  pdf(paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"best_",num,"_relative_p_clustmaximp_vals.pdf",sep=""))
  if(length(wps)){
    plot(prtot[-wps,],col=allcols[maxcols[-wps]],pch=pchs[-wps],cex=cexs[-wps],main="BEST CLUSTERING\n COL=HIGHEST PROB CLUSTER, POINT SIZE=PROB")
    points(prtot[startmodevals,],pch=18,col="black",cex=1.6)
    points(prtot[modevals,],pch=17,col="black",cex=2.0)
    points(prtot[wps,],col="black",pch=16)
  }else{
    plot(prtot,col=allcols[maxcols],pch=pchs,cex=cexs,main="BEST CLUSTERING\n COL=HIGHEST PROB CLUSTER, POINT SIZE=PROB")
    points(prtot[startmodevals,],pch=18,col="black",cex=1.6)
    points(prtot[modevals,],pch=17,col="black",cex=2.0)  
  }  
  dev.off()
  return(maxcols)
}

clprob4<-function(x,cla,startmodevals,startclustids){
  ddd<-as.matrix(dist(x)) 
  distprob<-matrix(0,nrow=dim(x)[1],ncol=length(unique(cla)))
  for(i in unique(cla)){
    jj<-which(startclustids==i)
    dists<-matrix(0,nrow=dim(x)[1],ncol=length(jj))
    for(j in 1:length(jj)){
      dists[,j]<-ddd[startmodevals[jj[j]],]
    }
    dist<-colMins(dists)
    dprob<-(max(dist)-dist)/max(dist)
    distprob[,i]<-dprob
  }
  return(distprob)
}
zscore<-function(x){
  xx<-which(! is.na(x))
  out<-(xx-mean(xx))/sd(xx)
  return(list(xx,out))
}
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  require(ggplot2)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of cluster or number of instances do not match")
  }
  
  nC <- length(idsA)
  tupel <- c(1:nA)
  
  # computing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
  
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  
  return(result)
}

remapclusters<-function(clusteringA,clusteringB) {
  matching <- minWeightBipartiteMatching(clusteringA,clusteringB)
  tmp <- sapply(1:length(matching), function(i) {
    clusterA[which(clusteringA == i)] <<- matching[i]
  })
  return(clusterA)
}
# Clustergram function
#' @param clusters data frame giving cluster assignments as produced by 
#'   many_kmeans or all_hclust
#' @param y value to plot on the y-axis.  Should be length
#'   \code{max(clusters$obs)}
clustergram <- function(clusters, y, line.width = NULL) {
  lo<-max(clusters$obs)
  #clusters$y <- y[clusters$obs+((clusters$i-1)*lo)]
  #clusters$out9 <-clusters$out9[clusters$obs+((clusters$i-1)*lo)]
  clusters$y <- y[clusters$obs]
  clusters$out9 <-clusters$out9[clusters$obs]
  # print(clusters$out9)
  clusters$center <- ave(clusters$y, clusters$i, clusters$cluster)
  clusters$out9 <-ave(clusters$out9,clusters$i,clusters$cluster)
  #plot(clusters$out9)
  mp<-dim(clusters)[1]/dim(y)[1]
  clustord<-NULL
  runind<-1
  for(j in 1:mp){
    valnums<-runind:(runind+dim(y)[1]-1)
    ucc<-unique(clusters$center[valnums])
    ords<-order(unique(clusters$center[valnums]))
    cl2new<-rep(0,dim(y)[1])
    for(i in 1:length(ucc)){
      xx<-which(clusters$center[valnums]==ucc[i])
      cl2new[xx]<-which(ords==i)
    }
    clustord<-c(clustord,cl2new)
    runind<-runind+dim(y)[1]
  }
  clusters$ord<-clustord
  if (is.null(line.width)) {
    line.width <- 0.5 * diff(range(clusters$center, na.rm = TRUE)) /
      length(unique(clusters$obs))
  }
  clusters$line.width <- line.width
  # Adjust center positions so that they don't overlap
  clusters <- clusters[with(clusters, order(i, center, y, obs)), ]
  clusters <- ddply(clusters, c("i", "cluster"), transform,
                    adj = center + (line.width * center(seq_along(y)))
  )
  structure(clusters,
            class = c("clustergram", class(clusters)),
            line.width = line.width)
}
plot.clustergram <- function(x) {
  i_pos <- !duplicated(x$i)
  means <- ddply(x, c("cluster", "i","ord"), summarise,
                 min = min(adj), max = max(adj),meanx=mean(out9))
  print(means)
  ggplot(x, aes(i)) +
    geom_ribbon(aes(y = adj, group = obs,ymin = adj - line.width/2, ymax = adj + line.width/2, colour = y)) +
    #geom_errorbar(aes(ymin = min, ymax = max), data = means, width =0.3,lwd=4*abs(means$meanx-1.5),col=round(means$meanx)) +
    geom_errorbar(aes(ymin = min, ymax = max), data = means, width =0.3,lwd=4*abs(means$meanx-1.5),col=means$ord) +
    scale_x_continuous("weights", breaks = x$i[i_pos], labels = x$k[i_pos]) +
    scale_color_gradientn(colors=ccx[1:9])+
    labs(y = "Cluster average", colour = "Obs\nvalue", fill = "Obs\nvalue")
}
#End of clustergram   
rename2<-function(x){
  name<-deparse(substitute(x))
  newname<-paste0(name,"b")
  return(list(newname,as.data.frame(x)))
}
peakweight<-function(x,y){
  if(x[1]==1){
  y[1]<-0  
  }
  rat<-max(y)/sum(y)
  return(rat)
}
maxgaps<-function(x){
  xs<-sort(x)
  ln<-length(xs)
  gaps<-rep(0,(ln-1))
  for(i in 2:ln){
    gaps[i-1]<-xs[i]-xs[i-1]
  }
  sortgaps<-sort(gaps)
  maxgaps2<-c(sortgaps[ln-1],sortgaps[ln-2],sortgaps[ln-3],sortgaps[ln-4])
  mg<-match(maxgaps2,gaps)
  mm1<-match(xs[(mg[1]+1):ln],x)
  mm2<-match(xs[(mg[2]+1):ln],x)
  #mm3<-match(xs[(mg[3]+1):ln],x)
  
  return(list(mm1,mm2))
}
jointclusts<-function(x){
  y<-unique(x)
  out<- match(data.frame(t(x)), data.frame(t(y)))
  return(out)
}
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
getdata<-function(x){
  datin<-x
  dataname<-deparse(substitute(x))
  outdat<-list(data=datin,name=dataname)
  class(outdat)<-"dataobj"
  return(outdat)
}
getdata2<-function(x){
  datin<-x[[2]]
  dataname<-x[[1]]
  outdat<-list(data=datin,name=dataname)
  class(outdat)<-"dataobj"
  return(outdat)
}

DENOISE<-1 # if ==1 ATTEMPTS OUTLIER IDENTIFICATION USING KNN_SUM 

real<-0 #  if ==1 ANALYZE REAL DATA , if ==0 SYNTHETIC CLUSTER WITH PRE_ASSIGNED LABELS

NEUROCOG<-3 # ANALYZE Neurocogntive data
REMONE<-1 # DO NOT ALLOW SINGLE CLUSTER if ==1, if ==0 ALLOW SINGLE CLUSTER
MULTIVIEW<-0 # PERFORM MULTIVIEW CLUSTERING ON REAL DATA
MULTISPLIT<-c(3,5) # SPLITTING OF INPUT DATA BETWEEN VIEWS
PAR<-1 # If ==1, PERFORM PARALLEL PROCESSING FOR SOME ROUTINES
FEAT<-0 # if ==1, DO FEATURE IMPOIRTANCE - QUITE SLOW
# NEXT TWO VARIABLES ARE FOR FEATURE IMPORTANCE
FEATINT<-1
SHOWDEP<-1

#alpha<-0.9
#snn<-0
smooth<-FALSE #if ==1, DO SMOOTHING IN DENSITYCUT ALGORITHM
# NEXT TWO LINES USED IF WANT TO DO WEIGHTING SEARCH, NORMALLY SET AT 1 TO DISABLE SEARCH, SLOW IF IMPLEMEMNTED
#allweights<-seq(0,1,0.05)
allweights<-1
# if ==1 NEXT VARIABLE ONLY IMPORTANT IN PROCESSING CERTAIN TYPES OF SYNTHETIC CLUSTER DATA WHERE x cooords, y coords and labels are on smame line
rename<-1
# If ==1 ALLOWS CONSTRUCTION OF SYNTHETIC DATA - ONLY IMPORTANT FOR TESTING
SYNTH<-0
# NEXT LINE DETERMINES WHICH DENSITY ESTIMATIONS ARE USED IN DENSITYCUT, 1 = KNN, 2 = SNN, 3 = RKNN
methods<-c(1,2,3)


if(exists("K")){
rm(K)
}
metallnames<-c("knn","snn","rknn")
metsused<-paste(metallnames[methods],collapse="_")
source("DensityCutMV5.R")

#####  compute data
moons2<-make_moons(750,0,0.10)
alphas<-c(0.1,0.2,0.3,0.4,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975)
set.seed(42)
mb6f15c_1.25<-make_blobs(1000,n_features=8,centers=15,cluster_std=1.25)
mb2f4c_0.400<-make_blobs(1000,n_features=2,centers=4,cluster_std=0.400)
mb3f6c_1.125<-make_blobs(1000,n_features=3,centers=6,cluster_std=1.125)
mb2f6c_0.25<-make_blobs(1000,n_features=2,centers=6,cluster_std=0.25)
mb30f20c_1.5<-make_blobs(1000,n_features=30,centers=20,cluster_std=1.5)
#mb1<-make_blobs(1000,n_features=10,centers=10,cluster_std=c(1,1,1,1,1,2,2,2,2,2))
mb1<-mb6f15c_1.25
#mb1<-moons2
mb2<-mb1
mb3<-mb1
#mb1_8D_15C<-mb1
#mb1<-make_blobs(n_samples=600,n_features=10,centers=5,cluster_std=c(1,2,1,2,1))
#mb1<-make_blobs(n_samples=600,n_features=5,centers=3,cluster_std=c(1.5,1.5,1.5))
#Spiral_squareB<-as.data.frame(Spiral_square)
#TWOD3Cb<-as.data.frame(TWOD3C)

#Sizes1<-read.csv("Sizes3"))
print("HERE")
if(real==0){
if(rename==1){
  #THE RENAME =1 OPTION IS FOR DATASETS WITH 3 COLUMNS, 1 and 2 ARE DATA, 3 IS LABELS
datain<-getdata2(rename2(Spiral_square))
}else{
  #THIS OPTION IS FOR OTHER DATASETS SUCH AS THOSE CREATED BY MAKE MOONS OR MAKE BLOBS OR DATASETS WITH > 2 DATA COLUMNS AND A SEP LABEL VARIABLE
  #ALSO FOR REAL DATA CONSISTING SOLELY OF VARIABLES WITH NO LABEL
datain<-getdata(Chainlink)
}
}

#datain<-getdata(Sizes5)
if(real==1){
  if(==3){
    wt<-1
    datplusid1<-read.csv("GNG_clust_all_24Sep2023.csv")
    datplusid2<-read.csv("RL_clust_all_24Sep2023.csv")
    datplusid3<-read.csv("SA_clust_all_24Sep2023.csv")
    
    commonids<-Reduce(intersect,list(datplusid1[,1],datplusid2[,1],datplusid3[,1]))
    wn1<-match(commonids,datplusid1[,1])
    wn2<-match(commonids,datplusid2[,1])
    wn3<-match(commonids,datplusid3[,1])
    
    #select the right number of columns depending on your variable selection
    #datain2<-cbind(datplusid1[wn1,1:3],datplusid2[wn2,2:3],datplusid3[wn3,2:3])
    
    datain2<-cbind(datplusid1[wn1,1:5],datplusid2[wn2,2:4],datplusid3[wn3,2:5])
    datain3<-cbind(datplusid1[wn1,7],datplusid2[wn2,6],datplusid3[wn3,7])
    data <- na.omit(datain2)
    
    
    wsn<-complete.cases(datain2)
    diagfinal<-datain3[wsn,3]
    NEUROCOGnewMULTRREAL<-datain2[wsn,]
    write.csv(NEUROCOGnewMULTRREAL, "all_tasks_script_export_22Nov2023.csv", row.names=FALSE)
    
    #write.csv(NEUROCOGnewMULTRREAL, "all_tasks_perf_script_export_22Nov2023.csv", row.names=FALSE)
    NEUROCOGnewMULTRREAL<- NEUROCOGnewMULTRREAL[,-1]
    
    datain<-getdata(NEUROCOGnewMULTRREAL)
    datname<-gsub("\\[.*?\\]", "",datain$name )
    data<-datain$data
    
    origdata<-data
    data2<-apply(data, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    
    data<-data2
    #Ensure that variables are coded in the correct manner
    data[,1] <- data[,1] * -1
    data[,3] <- data[,3] * -1
    data[,4] <- data[,4] * -1
    data[,5] <- data[,5] * -1
    dimd<-dim(data)
    #data2<-apply(data, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    #data2<-scale(data)
    #data<-data2
    data[,3]<-data[,3]*wt
    #data[,4]<-mb1orig[,4]*wt
    data[,5]<-data[,5]*wt
    #data[,6] <- data[,6] * -1
    #data[,7] <- data[,7] * -1
    #data[,11] <- data[,11] * -1
    
    mb1<-list("samples","labels")
    mb1$samples<-data
    mb1orig<-mb1$samples
    mb1$labels<-rep(1,dimd[1])
  }
  if(NEUROCOG==2){
    
    datplusid1<-read.csv("GNG_clust_all_24Sep2023.csv")
    datplusid2<-read.csv("RL_clust_all_24Sep2023b.csv")
    datplusid3<-read.csv("SA_clust_all_24Sep2023.csv")
    commonids<-Reduce(intersect,list(datplusid1[,1],datplusid2[,1],datplusid3[,1]))
    wn1<-match(commonids,datplusid1[,1])
    wn2<-match(commonids,datplusid2[,1])
    wn3<-match(commonids,datplusid3[,1])
    datain2<-cbind(datplusid1[wn1,1:5],datplusid1[wn1,2:5],datplusid1[wn1,2:5])
    wsn<-complete.cases(datain2)
    
    NEUROCOGnewMULTRREAL<-datain2[wsn,-1]
    datain<-getdata(NEUROCOGnewMULTRREAL)
    datname<-gsub("\\[.*?\\]", "",datain$name )
    data<-datain$data
    
    origdata<-data
    #data2<-apply(data, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    
    #data<-data2
    data[,1] <- data[,1] * -1
    data[,3] <- data[,3] * -1
    data[,4] <- data[,4] * -1
    data[,5] <- data[,5] * -1
    #data[,1] <- data[,1] * -1
    #data[,4] <- data[,4] * -1
    #data[,5] <- data[,5] * -1
    #data[,8] <- data[,8] * -1
    #data[,9] <- data[,9] * -1
    #data[,12] <- data[,12] * -1
    data[,3]<-data[,3]*wt
    #data[,4]<-mb1orig[,4]*wt
    data[,5]<-data[,5]*wt
    dimd<-dim(data)
    data2<-apply(data, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
    data<-data2
    for(i in 7:8){
      data[,i]<-data[,i]+runif(dimd[1],-0.1,0.1) 
    }
    for(i in 11:12){
      data[,i]<-data[,i]+runif(dimd[1],-0.2,0.2) 
    }
    mb1<-list("samples","labels")
    mb1$samples<-data
    mb1orig<-mb1$samples
    mb1$labels<-rep(1,dimd[1])
    
    
    
  }  
  if(NEUROCOG==1){
    
    #datplusid<-read.csv("full_ds copy.csv")
    datplusid<-read.csv("GNG_clust_all_18thSep2023.csv")
    #datplusidold<-read.csv("bothtasks_complete.csv")
    
    wsn<-complete.cases(datplusid)
    #wsn2<-complete.cases(datplusidold)
    NEUROCOGnew4<-datplusid[wsn,]
    
    #BTCnew<-datplusidold[wsn2,]
    
    #testmatch<-match(unique(BTCnew[,1]),BTCnew[,1])
    #BTCnew<-BTCnew[testmatch,]
    #uniques<-match(unique(NEUROCOGnew4[,1]),NEUROCOGnew4[,1])
    #NEUROCOGnew4<-NEUROCOGnew4[uniques,]
   
    datain<-getdata(NEUROCOGnew4)
    
    #datainold<-getdata(BTCnew)
   
    datname<-gsub("\\[.*?\\]", "",datain$name )
    
    #dataold<-datainold$data[,2:8]
    #subsold<-datainold$data[,1]
    #subs<-datain$data[,1]
    
    #data<-datain$data[,c(4,5,6,7,8,9,10,11,12)]
    data<-datain$data[,2:5]
    
    #dataold[,c(3,7)]<-(-(dataold[,c(3,7)]))
    
    
    
    datname<-gsub("\\[.*?\\]", "",datain$name )
    
    #matches<-match(subs,subsold)
    #wnnn<-which(! is.na(matches))
    #data[,c(1,2,3,6,7,10)]<-(-(data[,c(1,2,3,6,7,10)]))
    origdata<-data
    data<-scale(data)
    #data[,3] <- data[,3] *0.75
    #data[,4] <- data[,4] *0.75
    data[,1] <- data[,1] * -1
    data[,4] <- data[,4] * -1
    
    
    #dataold<-scale(dataold)
    #data[,c(1,4,8)]<-data[,c(1,4,8)]*10
    dimd<-dim(data)
    #dimd2<-dim(dataold)
      mb1$samples<-data
      mb1orig<-mb1$samples
      mb1$labels<-rep(1,dimd[1])
      if(MULTIVIEW==1){
      mb2<-mb1  
      mb1$samples<-data[,1:MULTISPLIT[1]] 
      mb2$samples<-data[,(MULTISPLIT[1]+1):(MULTISPLIT[2])]
      mblist<-list(mb1$samples,mb2$samples)
      }
    
    
  }
  if(NEUROCOG==0){
  #test2<-read.csv("LEAP_t2_Core clinical variables_12-03-19-withlabels.csv")
  test2[test2==999]<-NA
  test2[test2==777]<-NA
  subnums<-test2[,1]
  dimorig<-dim(test2)
  #test2b<-as.data.frame(test2[,-(c((1:7),32,41,53))])
  test2b<-as.data.frame(test2[,-(c((1:7)))])
  dimt2b<-dim(test2b)
  #First 2 vars, cols 16 and 17 of test2b
  xx1<-which(!is.na(test2b[,16]))
  xx2<-which(!is.na(test2b[,17]))
  yy1<-which(!is.na(test2b[,19]))
  yy2<-which(!is.na(test2b[,20]))
  yy3<-which(!is.na(test2b[,21]))
  zz1<-which(!is.na(test2b[,28]))
  zz2<-which(!is.na(test2b[,29]))
  aa1<-which(!is.na(test2b[,30]))
  aa2<-which(!is.na(test2b[,31]))
  bb1<-which(!is.na(test2b[,32]))
  bb2<-which(!is.na(test2b[,33]))
  cc1<-which(!is.na(test2b[,34]))
  cc2<-which(!is.na(test2b[,36]))
  newdata<-matrix(NA,nrow=764,ncol=6)
  #newdata[xx1,1]<-zscore(test2b[xx1,16])[[2]]
  newdata[xx1,1]<-test2b[xx1,16]
  #newdata[xx2,1]<-zscore(test2b[xx2,17])[[2]]
  newdata[xx2,1]<-test2b[xx2,17]
  #newdata[c(xx1,xx2),1]<-zscore(newdata[c(xx1,xx2),1])[[2]]
  #newdata[yy1,2]<-zscore(test2b[yy1,19])[[2]]
  newdata[yy1,2]<-test2b[yy1,19]
  #newdata[yy2,2]<-zscore(test2b[yy2,20])[[2]]
  newdata[yy2,2]<-test2b[yy2,20]
  #newdata[yy3,2]<-zscore(test2b[yy3,21])[[2]]
  newdata[yy3,2]<-test2b[yy3,21]
  #newdata[c(yy1,yy2,yy3),2]<-zscore(newdata[c(yy1,yy2,yy3),2])[[2]]
  #newdata[zz1,3]<-zscore(test2b[zz1,28])[[2]]
  newdata[zz1,3]<-test2b[zz1,28]
  #newdata[zz2,3]<-zscore(test2b[zz2,29])[[2]]
  newdata[zz2,3]<-test2b[zz2,29]
  #newdata[aa1,4]<-zscore(test2b[aa1,30])[[2]]
  #newdata[c(zz1,zz2),3]<-zscore(newdata[c(zz1,zz2),3])[[2]]
  newdata[aa1,4]<-test2b[aa1,30]
  #newdata[aa2,4]<-zscore(test2b[aa2,31])[[2]]
  newdata[aa2,4]<-test2b[aa2,31]
  #newdata[bb1,5]<-zscore(test2b[bb1,32])[[2]]
  #newdata[c(aa1,aa2),4]<-zscore(newdata[c(aa1,aa2),4])[[2]]
  newdata[bb1,5]<-test2b[bb1,32]
  #newdata[bb2,5]<-zscore(test2b[bb2,33])[[2]]
  newdata[bb2,5]<-test2b[bb2,33]
  #newdata[cc1,6]<-zscore(test2b[cc1,34])[[2]]
  #newdata[c(bb1,bb2),5]<-zscore(newdata[c(bb1,bb2),5])[[2]]
  newdata[cc1,6]<-test2b[cc1,34]
  #newdata[cc2,6]<-zscore(test2b[cc2,36])[[2]]
  newdata[cc2,6]<-test2b[cc2,36]
  #newdata[c(cc1,cc2),6]<-zscore(newdata[c(cc1,cc2),6])[[2]]
  #xn<-missForest(test2b)
  #datain<-getdata(test2b[,16:56])
  #datain<-getdata(newdata[,c(1,2,3,5)])
  #datain<-getdata(newdata[,c(1,2,3,5)])
  #datain<-getdata(newdata[,c(1,2)])
  #datain<-getdata(newdata[,c(1,2,4)])
  newdata<-matrix(NA,nrow=764,ncol=3)
  newdata<-test2b[,c(14,22,23)]
  
  #datain<-getdata(newdata[,c(1,4,6)])
  datain<-getdata(newdata)
  xxb<-which(complete.cases(datain$data)==TRUE)
  datain$data<-datain$data[xxb,]
  #xx<-apply(test2b,2,function(x) length(which(complete.cases(x)==TRUE)))
  #xdat1<-as.matrix(read.table("XDAT1TOM"))
  #xdat1b<-as.data.frame(xdat1)
  dimreal<-dim(datain$data)
  dt<-dim(datain$data)
  mb1$samples<-datain$data
  mb1$labels<-as.numeric(test2[xxb,2])
  }
}else{

dt<-dim(datain$data)
if(rename==1){
  names(datain$data)<-c("x","y","label")  
}
name<-datain$name
#datain$name="blobs"
#mb1$samples<-as.matrix(Chainlink$Data)
#if(rename==0){
#names(datain$data)<-seq(1,dt[2],1)
#}
lx<-length(names(datain$data))
nd<-names(datain$data)
if(names(datain$data[1])!="Data"){
  wl<-which(names(datain$data)=="label")
  wl2<-which(names(datain$data)=="labels") 
  if(sum(wl,wl2)){
    if(length(wl)!=0 && wl>2){
      mb1$samples<-as.matrix(datain$data[,1:(wl-1)])
      mb1$labels<-datain$data[wl]
    }
    if(length(wl2)!=0 && wl2 >2){
      mb1$samples<-as.matrix(datain$data[,1:(wl2-1)])
      mb1$labels<-datain$data[wl]
    }
    if(length(wl2)!=0 && wl2 ==2){
      mb1$samples<-as.matrix(datain$data$samples)
      mb1$labels<-datain$data$labels
    }
    
    
  }else{
    types<-sapply(datain$data,class)
    labs<-which(types=="integer")
    labs2<-which(types=="factor")
    
    
    if(length(labs)!=0){
      mb1$samples<-as.matrix(datain$data[,1:(labs-1)]) 
      mb1$labels<-datain$data[labs]
    }
    if(length(labs2)!=0){
      mb1$samples<-as.matrix(datain$data[,1:(labs2-1)]) 
      mb1$labels<-as.integer(unlist(datain$data[labs2])) 
    }
    if(length(labs2)==0 && length(labs) ==0){
      mb1$samples<-as.matrix(datain$data)
      mb1$labels<-NULL
    }
  }
}else{
  mb1$samples<-as.matrix(datain$data$Data)
  mb1$labels<-datain$data$Cls
}


if(SYNTH == 1 && rename ==1){
  mb1$labels<-datain$data[,3]
  mb1$labels<-as.numeric(mb1$labels)
  if(min(mb1$labels)==0){
    mb1$labels<-mb1$labels+1
  }
} 
if(SYNTH == 0 && rename ==0){
  if(is.list(mb1$labels)){
  mb1$labels<-as.numeric(unlist(mb1$labels))
  }
  if(min(mb1$labels)==0){
    mb1$labels<-mb1$labels+1
  }
} 
}
#mb1$samples<-Ring$ring
#mb1$labels<-Ring$ring.mem
dsorig<-dim(mb1$samples)


if(NEUROCOG==0){
mb1$samples<-scale(mb1$samples)
}
dimd<-dim(mb1$samples)

allprs<-array(0,c(dimd[1],2,length(allweights)))

  
  #mb1$samples[,3]<-mb1orig[,3]*wt
  #mb1$samples[,4]<-mb1orig[,4]*wt

#plot(pra,col=d1$cluster,pch=16,main="VIEW1")
NUMPERM<-length(alphas)*3
nu<-seq(0,1,0.01)
lnu<-length(nu)
prtot<-prcomp(mb1$samples)$x[,1:2]
index<-1
K = ceiling(log2(dimd[1]))
if (K %% 2 != 0) {
  K = K + 1
}
Korig<-K

#ALLKS<-c((K-5),(K-4),(K-3),(K-2),(K-1),K,(K+1),(K+2),(K+3),(K+4),(K+5))
ALLKS<-K

b<-lnu
Kindices<-sapply(ALLKS, function (x) rep(x,b))
Kindices<-as.vector(Kindices)
methods<-c(1,2,3)

ds<-dim(mb1$samples)
prtot<-prcomp(mb1$samples)$x[,1:2]

superclust<-matrix(0,length(ALLKS)*lnu*3*length(alphas),dimd[1])
lengthall<-length(ALLKS)*lnu*3*length(alphas)
  indvalk=which(ALLKS==K)
  print(paste("K = ",indvalk,sep=""))
 #Compute all clusterings using DensityCut, return clustering for all nu values for all methods for all alphas for all Ks
testxxx<-foreach(m = 1:NUMPERM)%dopar%{
  s<-m%%3
  ss<-s
  if(s==0){s<-3}
  a <- m%/%3
  if(ss==0){a<-a-1}
  a<-a+1
  ll<-list()
  foreach(KK = 1: length(ALLKS))%do%{
    K<-ALLKS[KK] 
  
  if(MULTIVIEW==0){
    
    dall<-DensityCut4MV(mb1$samples,alpha=alphas[a],snn=s-1,smooth=smooth,K=K,nu=nu)
     }else{
    dall<-DensityCut4MV(mblist,alpha=alphas[a],snn=s-1,smooth=smooth,K=K,nu=nu) 
     }
  dall$allclusts
  }
  
}
print("COMPLETED K LOOP")
index<-1
#Collect subject cluster allocations for all clusterings in superclust
for(i in 1:NUMPERM){
 nnn<-testxxx[[i]] 
  for(j in 1:length(nnn)){
    nnnn<-as.matrix(nnn[[j]])
    superclust[index:(index+lnu-1),]<-t(nnnn)
    index<-index+lnu
  }
  
  
}
print("COMPLETED CALCULATION OF SUPERCLUST") 
# Find number of clusters in each clustering
superclustnums<-apply(superclust, 1, function(x) length(unique(x)) )
#Calculate frequency distribution of cluster numbers
tt<-table(superclustnums)
tt1<-as.numeric(names(tt))
tt2<-tt

if(REMONE){# Remove single cluster solutions if REMONE = 1
  wuscn1<-0
if(min(tt1)==1){
  
wuscn1<-which(tt1==1)  
}
  if(wuscn1){
    print("HERE")
 tt1<-tt1[-wuscn1] 
 tt2<-tt2[-wuscn1]
  }
} 
clustnumtab<-matrix(0,length(tt1),2)
clustnumtab[,1]<-tt1
clustnumtab[,2]<-tt2
#Save cluster frequency distribution data
write(t(clustnumtab),file=paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"_table_of_clustnums_and_frequencies",sep=""),ncolumns=2)
pdf(paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"_allclustsbyfequency.pdf",sep=""))
# Plot frequency distribution data
plot(tt1,tt2,xlab="cluster number",ylab="frequency",type="lines")

ocww<-order(tt2)
ocw1<-tt1[ocww]
ocw2<-tt2[ocww]
locw<-length(ocw1)
#Compute best three cluster numbers (max frequencies)
maxclusts<-ocw1[locw:(locw-2)]
tables<-list()
wms<-rep(0,3)
listallwms<-rep(list(list()),3)
# Loop through top three cluster solutions
for(i in 1:3){ # START OF LOOPS TO LOOK AT BEST,  2nd BEST and 3rd BEST CLUSTER SOLUTION

# Find index numbers in superclust for clusterings with required number of clusters  
a11<-findsupermatchnew2(superclustnums,lengthall,maxclusts[i] ) 
print(paste("COMPLETED SUPERMATCH LOOP ",i,sep=""))
clustss<-matrix(0,nrow=length(a11),ncol=dimd[1])
vvv<-as.vector(unlist(a11))

for(k in 1:length(a11)){
 vv<-a11[[k]]
 clustss[k,]<-changeclustnums(superclust[vv[1],])
}
print(paste("END OF CHANGECLUSTERNUMS LOOP ",i,sep=""))
#END OF LOOP


cors<-NULL
#Compute matrix of variations on clustering for given cluster number (ucc) with frequency of each variation (freqcss)
ax<-count.duplicates2(clustss)

ucc<-apply(ax[,1:((dimd[1]))],2,function(x) unlist(x))
print(paste("END OF COUNT DUPLICATES AND UNLISTING LOOP ",i,sep=""))
if(is.matrix(ucc) && dim(ucc)[2]>1){
freqcss<-rep(0,dim(ucc)[1])
}else{
  freqcss<-rep(0,length(ucc))
}
if(is.matrix(ucc)&& dim(ucc)[2]>1){
cors<-cor(t(ucc))
corsorig<-cor(t(ucc))
}
freqcss<-ax[,dimd[1]+1]
 wmo<-order(freqcss)
 lf<-length(freqcss)
 wmos<-wmo[lf]  
 wm<-which(freqcss==max(freqcss)) 
 if(length(wm>1)){
   wm<-wm[1]
 }
 frq1<-freqcss[wm]
 listallwms[[i]]<-freqcss
 #f1<-ucc[wmos[1],]
 #f2<-ucc[wmos[2],]
 # f1 is most frequent variation
 f1<-unlist(ax[wm,1:dimd[1]])
 uccorig<-ucc
 # Remap all clusterings onto most frequent variation to produce consistent cluster numbers.
 if(lf>1){
 clusterA<-rep(0,dimd[1]) 
 if(PAR==0){
 for(xxxx in 1:(dim(ucc)[1]))
 {
   #print(xxxx)
   newc<-remapclusters(ucc[xxxx,],f1)
   ucc[xxxx,]<-newc
   
 }
 }else{  
 uccnew<-ucc
 testzzz<-foreach(xxxx = 1:(dim(ucc)[1]))%dopar%{
 
   #print(xxxx)
   newc<-remapclusters(ucc[xxxx,],f1)
   uccnew[xxxx,]<-newc
   uccnew[xxxx,]
 }
 for(ii in 1:length(testzzz)){
 ucc[ii,]<-testzzz[[ii]]  
 }
 }
 print(paste("COMPLETED REMAPCLUSTERS LOOP ",i,sep=""))
 if(is.matrix(ucc)&& dim(ucc)[2]>1){
   cors<-cor(t(ucc))
   
 }
 }
 
 

# Calculate cluster probability data for each subject ( frequencies of membership of each cluster for each subject)
 cexs<-rep(1,dimd[1])
 if(PAR==0){
 cc<-Modefunc3(ucc,freqcss)
 
 }else{
   cc<-Modefunc3_par(ucc,freqcss)
   
 }
 print(paste("COMPLETED MODEFUNC LOOP ",i,sep=""))
 cexs<-cc[[1]]/sum(freqcss)
 probs<-cc[[2]]

 clustxs<-f1
 ux<-unique(clustxs)
 cvxs<-rep(0,dimd[1])
 for(j in 1:length(ux)){
   
   cvxs[clustxs==ux[j]]<-j
   
 }
 clsing<-cc[[3]]
# ff1<-f1[1]
 #ff2<-f2[2]
 
 clustmodclust<-matrix(0,dimd[1],3)
 rownums<-seq(1,dimd[1],1)
 clustmodclust[,1]<-rownums
 clustmodclust[,2]<-clustxs
 clustmodclust[,3]<-cc[[3]]
 # Write frequency data
 write(t(clustmodclust),file=paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"_clustering_v_mod_clustering",i,sep=""),ncolumns=3)
 
allcolx<-rainbow(length(ux))
# Start plotting of clusterings, point size = probability
plot(prtot,col=allcolx[clsing],cex=cexs,pch=16,main=paste("MAX CORRESPONDENCE CLUSTER NUMBER = ",maxclusts[i],"\n","FREQUENCY = ",frq1,sep=""))

if(DENOISE==1){
 # Calculate "noise" points if desired 
  KNOISE<-round(dimd[1]^0.2)
  xxx<-KNN_SUM(mb1$samples,KNOISE)
  numseq<-seq(1,dimd[1],1)
  out1<-numseq[which(xxx>3)]
  allcolsnew<-rainbow(max(unique(clsing))+1)
  allcolsnew[max(unique(clsing))+1]<-"black"
  clanew<-clsing
  clanew[out1]<-max(unique(clsing))+1
  plot(prtot,col=allcolsnew[clanew],pch=16,main=paste("BEST OVERALL CORRESPONDENCE CLUSTER NUMBER = ",maxclusts[i],"\n","FREQUENCY = ",frq1,"\n","OUTLIERS IN BLACK",sep=""),cex.main=0.9)
}

if(lf>1){
  plot(listallwms[[i]],type="lines",main=paste("FREQUENCY OF CORRESPONDENCES TOTAL = ",sum(unlist(listallwms[[i]])),sep=""),ylab="frequency")
  par(new=TRUE)
  plot(cors[wmos[1],],col="red",type="lines", axes = FALSE, xlab = "", ylab = "",ylim=c(0,1))
  axis(side = 4)   
  par(new=TRUE)
  plot(corsorig[wmos[1],],col="green",type="lines", axes = FALSE, xlab = "", ylab = "",ylim=c(0,1))
  
}
if(is.matrix(cors)){
freqcor<-matrix(0,3,length(freqcss))
freqcor[1,]<-freqcss
freqcor[2,]<-cors[wmos[1],]
freqcor[3,]<-corsorig[wmos[1],]
write(t(freqcor),file=paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"_clustering_",i,sep=""),ncolumns=length(freqcss))
}
  if(DENOISE==0){
  clmodev<-matrix(0,dimd[1],dim(probs)[2]+3)
  }else{
    clmodev<-matrix(0,dimd[1],dim(probs)[2]+4) 
  }
  clmodev[,1]<-rownums
  clmodev[,2]<-cvxs
  clmodev[,3]<-cc[[3]]
  clmodev[,4:(dim(probs)[2]+3)]<-probs
  if(DENOISE){
  chars<-rep("O",dimd[1])
  clmodev[,(dim(probs)[2]+4)]<-chars
  clmodev[out1,(dim(probs)[2]+4)]<-"X"
  }
  # Compare feature importance if desired.
  write(t(clmodev),file=paste(datain$name,"_methods_used_",metsused,"_smooth_",smooth,"_clustering_modesums_",i,sep=""),ncolumns=dim(clmodev)[2])
  
  if(FEAT){
    source("PYRFPERM2.R")
    out<-Pypermvarimp2(mb1$samples,clsing,datain$name,0,FEATINT,SHOWDEP,i)
  }

if(real && NEUROCOG==0){
tables[[i]]<-table(test2[xxb,2],cvxs)
}
}
dev.off()

#uniques<-rep(0,dimd[1])
#for(i in 1:dimd[1]){
#  uniques[i]<-length(unique(superclust[,i]))
#}

#newuniques<-changeclustnums(uniques)
#newcols<-rainbow(length(unique(uniques)))