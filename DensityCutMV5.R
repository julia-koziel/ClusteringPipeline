##====================================================================
# densityCut v1.0
#
# Jiarui Ding (jiaruid@cs.ubc.ca)
# Department of Computer Science, The University of British Columbia
# Department of Molecular Oncology, BC Cancer Reseach Centre
# 
# Novmber 25, 2015
# 

##====================================================================
# Fast version of finding local-maxima, Time:O(NK), Space:O(NK)
# Using either out-nodes or in-nodes
#
library(rnndescent)
library(dbscan)
source("SNNCOMBL.R")
library(DDoutlier)
RKNNLIST <-function(X,N,K){
 ll<-length(X) 
 intmat<-matrix(0,nrow=N,ncol=ll)
 for(i in 1:ll){
 intmat[,i]<- KNN_IN(X[[i]],K)
 }
 outmat<-rep(0,N)
 for(i in 1:N){
   outmat[i]<-max(intmat[i,])
 }
 return(outmat) 
}
RKNNLIST2 <-function(X,N,K){
  
  if(is.list(X)){
  ll<-length(X)  
  }else{
    ll<-1
  }
  intmat<-array(0,c(N,K,ll))
  for(i in 1:ll){
    for(m in 1:K){
      out<-rep(0,N)
      if(ll>1){
      test <- KNN_IN(X[[i]],m)
      }else{
       test<- KNN_IN(X,m)
      }
      intmat[1:length(test),m,i]<-test
    }
  }
  outmat<-matrix(0,nrow=N,ncol=K)
  for(i in 1:N){
    for(m in 1:K){
      outmat[i,m]<-max(intmat[i,m,])
    }
  }
  return(outmat) 
}

knnDE2 <-function(X, Grid, k){
  
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k should be a positive integer")
  }
  
  d <- ncol(X)
  n <- nrow(X)
  r.k <- apply(FNN::knnx.dist(X, Grid, k = k, algorithm = "kd_tree"), 1, max)
  v.d <- pi^(d/2) /gamma(d/2+1)
  out <- k / (n * v.d * r.k ^ d)  
  return(out)
}
CheckLocalMax = function(knn, index, V) {
  N = length(knn)
  
  local.maxima = sapply(seq(N), function(z) {
    id = knn[[z]][, 1]
    z1 = index[z]
    
    if (all(V[id] <= V[z1])) {
      return(z1)
    } else {
      return()
    }
  })
  
  return(unlist(local.maxima))
}


##====================================================================
# Fast calculate in-neighbours
GetInNeighbour = function(knn.index.row, knn.index.col, distance) {
  len = length(knn.index.col)
  
  id  = order(distance)
  id1 = order(knn.index.col[id])
  id  = id[id1]
  
  knn.index.col = knn.index.col[id]
  knn.index.row = knn.index.row[id]
  distance = distance[id]
  
  knn.index.col.minus = knn.index.col[-1]
  knn.index.col.minus = c(knn.index.col.minus, knn.index.col[len])
  
  ## because of float-point data, 
  #  and the last one should be an end index
  id.end = which(knn.index.col.minus - knn.index.col >= 0.5)
  id.end = c(id.end, len)
  
  N = length(id.end)
  id.start = c(1, id.end+1)[1:N]
  
  nn = lapply(seq(N), function(z) {
    id = id.start[z]:id.end[z]
    cbind(knn.index.row[id], distance[id])
  })
  index = knn.index.col[id.end]
  
  return(list(knn=nn, index=index))
}


##====================================================================
# Time O(NK), Space O(NK)
# knn[[z]] is sorted ascendently based on distance 
# Either in-neighbours or out-neighbours
# 
NearestNeighbour = function(knn, index, V, id) {
  N = length(knn)
  
  if (missing(id)) {
    id = seq(N)
  }
  
  nearest.neighbour = lapply(id, function(z) {
    x  = knn[[z]]
    id = x[, 1]
    v  = V[id]
    
    z1 = index[z]
    v.diff  = (v - V[z1]) #/ (x[, 2]+.Machine$double.eps)
    id.high = v.diff > 0
    
    id.nn = NA
    if (sum(id.high) > 0) {
      # id.nn.high = which.max(v.diff)
      id.nn.high = which(id.high)[1]
      id.nn = id[id.nn.high]
    }
    
    return(id.nn)
  })
  
  return(nearest.neighbour)
}


##====================================================================
# Given a knn-graph and the estimated densities, initial clustering
# data and detecting valley separating clusters
# 
# Vertices are labeld from 1 to N
# 
# Worst-case memory: O(C^2 + N), Time: O(NK) 
#
AssignCluster = function(nearest.neighbour, knn, index, V, mode, 
                         N, K, adjust=TRUE) {
  M = length(mode)
  
  mode.name = seq(M)
  value = vector("list", M)
  for (i in mode.name) {
    value[[i]] = list()
  }
  cluster = value
  
  # Initialization, assign points to modes
  cluster.assign = rep(0, M)
  V.assign = rep(0, N)
  
  id.mode = order(V[mode], decreasing=TRUE)
  V.assign[mode[id.mode]] = mode.name
  
  id = order(V[index], decreasing=TRUE)
  iter = 0
  for (id.i in id) {
    iter = iter + 1
    
    id.i1 = index[id.i] # convert to the actual coordinate
    label = V.assign[nearest.neighbour[[id.i]]]
    V.assign[id.i1] = label
    cluster.assign[label] = cluster.assign[label] + 1
    if (label == 0) { # the only possibility is outliers
      next
    }
    
    # This is the interesting part. 
    # A new bouhelp(ndary point - neighbours with higher densities 
    # than V[id.i] and labeled differently 
    knn.index = knn[[id.i]][, 1]
    
    id.high   = which(V[knn.index] >= V[id.i1]) # sub index
    if (length(id.high) > 0) {
      id.high.new = knn.index[id.high]
      assign.boundary = V.assign[id.high.new]
      
      # New adjacent boundary cluster. remove outliers
      known.boundary = c(0, label, unlist(cluster[[label]]))
      cluster.index  = !(assign.boundary %in% known.boundary)
      if (sum(cluster.index) <= 0) {
        next
      }
      # Cluster does not include the current label
      cluster.uniq = unique(assign.boundary[cluster.index])
      
      if (adjust == FALSE) {
        # print("no adjust")
        out = sapply(cluster.uniq, function(z) V[id.i1])
      } else {
        # print("adjust")
        valley.adjust = median(V[knn.index])
        
        if (valley.adjust > V[id.i1]) {
          weight = 2 - iter / N
        } else {
          weight = 1
        }
        out = sapply(cluster.uniq, function(z) V[id.i1] * weight)
      }
      
      len = length(value[[label]])
      for (ii in seq(length(cluster.uniq))) {
        id = len + ii
        value[[label]][[id]] = out[[ii]]
        cluster[[label]][[id]] = cluster.uniq[[ii]]
      }
    }
  }
  
  return(list(V.assign=V.assign, 
              cluster.assign=cluster.assign, 
              valley=list(cluster=cluster, value=value),
              id.mode=id.mode)
  )
}


##====================================================================
# Enhance densities based on the transition matrix
EnhanceDensity = function(P, V, smooth=FALSE, debug=TRUE,
                          maxit=50, alpha=0.85, eps=1e-5) {
  iter = 0
  done = FALSE
  
  V0 = V
  while (!done) {    
    iter = iter + 1
    if (smooth == TRUE) {
      V1 =  alpha * P %*% V  + (1-alpha) * V0
    } else {
      V1 =  alpha * V %*% P  + (1-alpha) * V0
    }
    V1 = as.vector(V1)
    V1 = V1 / sum(V1)
    
    v.diff = sum(abs(V1 - V))
    if (debug == TRUE) {
      cat("Iter: ", iter, " V-diff: ", v.diff, "\n")
    }
    
    if (v.diff <= eps) {
      done = TRUE
    }  else if (iter > maxit) {
      cat(paste("WARNING! not converged"), "\n")
      done = TRUE
    }
    V = V1
  }
  
  return(as.vector(V))
}


##====================================================================
MergeCut = function(valley, cluster.assign, V.local.maxima, 
                    nu, show.plot, tip.color, show.tip.label, 
                    text, xlab) {
  M = length(valley$cluster)
  index = which(sapply(valley$cluster, length) > 0) 
  
  cluster = lapply(valley$cluster, function(z) do.call(rbind, z))
  value   = lapply(valley$value, function(z) do.call(rbind, z))
  
  ## Loop through each cluster, find the valley height
  contrast = lapply(index, function(z) {
    valley.a = value[[z]]
    cluster.id = cluster[[z]]
    
    ## Change from min to max -- condition specific?
    valley.height = sapply(seq(length(cluster.id)), function(x) {
      id = cluster.id[x]
      id.overlap = which(cluster[[id]] %in% z)
      
      if (length(id.overlap) > 0) {
        xx = max(value[[id]][id.overlap], valley.a[x])
      } else {
        xx = valley.a[x]
      }
    })
    return(valley.height)
  })
  
  ## Re-label local-maxima and change V, do not change contrast
  MergeClusterThreshold = function(nu, merge.order) {
    for (id.index in merge.order) {
      done  = FALSE
      cluster.id = index[id.index]
      point.boundary = contrast[[id.index]]
      
      ratio = point.boundary / 
        pmin(V[cluster.id], V[cluster[[cluster.id]]])
      
      id = cluster[[cluster.id]][which(ratio > nu)]
      
      ## Use label to find the un-merged clusters
      while (!done) {
        if (length(id) >= 1) {
          id.merge = which(label[id] != label[cluster.id])
          if (length(id.merge) == 0) {
            done = TRUE
          } else {
            id.merge.order = order(ratio[id.merge], decreasing=TRUE)
            id = id[id.merge[id.merge.order][1]]
          }
          
          merged.cluster = c(id, cluster.id)
          label.merge = unique(label[merged.cluster])
          if (length(label.merge) > 1) {
            min.lab = min(label.merge)
            
            id = which(label %in% label.merge)
            label[id] = min.lab
            
            V[id] = max(V[id])   
          }
        } else {
          done = TRUE
        }
        
        ratio = point.boundary / 
          pmin(V[cluster.id], V[cluster[[cluster.id]]])
        
        id = cluster[[cluster.id]][which(ratio > nu)]
      } # End while
    } # End for
    return (list(label, V))
  }
  
  ## Merging by re-setting the heights and labels of local-maxima 
  label = as.numeric(names(V.local.maxima))
  V = V.local.maxima
  
  label.all = vector("list", length(nu))  
  iter = length(nu)
  
  if (length(index) > 0) {
    merge.order = seq(length(index))
  } else {
    merge.order = NULL
  }
  
  CheckNeighbour = function() {
    id = sapply(index[merge.order], function(z) {
      any(label[cluster[[z]][,1]] != label[z])
    })
    if (is.list(id)) {
      id = unlist(id)
    }
    merge.order = merge.order[id]
    return(merge.order)
  }
  
  #====
  for (nu.i in rev(nu)) {
    merge.order = CheckNeighbour()
    out = MergeClusterThreshold(nu.i, merge.order)
    label = out[[1]]
    label.all[[iter]] = label
    V = out[[2]]
    
    iter = iter - 1
  }
  label = label.all
  merged.label = do.call("cbind", label)
  colnames(merged.label) = nu
  
  ## Remove the labels to prevent generating a collaped tree
  if (show.plot == TRUE) {    
    if (length(label[[length(label)]]) > 1) {
      id = sapply(label, function(z) !all(z == z[1]))
      id = which(id)
      merged.label.plot = merged.label[, id, drop=FALSE]
      
      PlotDensitycut(merged.label.plot, 
                     base=nu[1], 
                     xlab=xlab, 
                     tip.color=tip.color,
                     show.tip.label=show.tip.label)
      mtext(text, side=3, line=0.5, adj=-0.08, cex=0.8, col="black")
    } else {
      plot.new()
      text(0.5, 0.5, label="Only one cluster!", cex=0.8)
    }
  }
  
  return(merged.mode=merged.label)
}


SelectCluster = function(label, show.plot=TRUE, xlab=TRUE) {
  ## Determine the number of cluster and centers
  cluster.count = apply(label, 2, function(z) length(unique(z)))
  
  x = table(cluster.count)
  name  = as.numeric(names(x))
  len   = length(x)
  
  SelectClusterNumber = function(x) {
    name  = as.numeric(names(x))
    
    x.max = max(x)
    level = which(x == x.max)
    level = level[length(level)]
    
    level = which(cluster.count == name[level])
    level = names(level)[1]
  }
  
  #==
  if (length(x) > 1) {
    xx = name
    xx.end = xx[-1]
    xx = xx[-len]
    id = xx.end - xx
    
    id.start = which(id <= 0)
    if (length(id.start) > 0) {
      id.end   = id.start + 1
      id.other = setdiff(seq(len), c(id.start, id.end))
      
      len.increase = length(id.start)
      if (len.increase > 1) {
        id.intermediate = 
          id.start[2:len.increase] == id.end[1:(len.increase-1)]
        
        id.intermediate = which(id.intermediate)
        
        if (length(id.intermediate) > 0) {
          id.end   = id.end[-id.intermediate]
          id.start = id.start[-(id.intermediate+1)]   
        }     
      }
      
      id.keep = id.start
      for (z in seq(length(id.start))) {
        id = id.start[z]:id.end[z]
        
        max.x = max(x[id])
        id.max = which(x[id] == max.x)
        id.max = id.max[length(id.max)]
        
        id.keep[z] = id[id.max]
        
        x[id.keep[z]] = sum(x[id])
      }
      
      id = setdiff(seq(len), c(id.keep, id.other))
      if (length(id) > 0) {
        x = x[-id]
      } 
    }
  }
  level = SelectClusterNumber(x)
  
  if (show.plot == TRUE) {
    space = 6 / length(x)
    if (space < 0.5) {
      space = 0.5
    }
    
    if (length(x) == 1) {
      bar = barplot(x/sum(x), xlab=NA, 
                    col=AddTrans("dodgerblue4", 0.45), 
                    ylab=NA,
                    cex.axis=0.8,
                    space=2,
                    width=0.1,
                    xlim=c(0,1),
                    cex=0.8,
                    xaxt="n",
                    yaxt="n")
    } else {
      bar = barplot(x/sum(x), xlab=NA, 
                    col=AddTrans("dodgerblue4", 0.45), 
                    ylab=NA,
                    cex.axis=0.8,
                    space=space,
                    width=0.1,
                    cex=0.8,
                    xaxt="n",
                    yaxt="n")
    }
    
    axis(side=2, tck=-0.015, labels=NA)
    axis(side=2, lwd=0, mgp=c(3,0.5,0), line=-0.4, 
         labels=TRUE, cex.axis=0.8)
    
    axis(side=1, tck=-0.015, labels=NA, at=bar)
    axis(side=1, lwd=0, mgp=c(3,0.5,0), line=-0.4, at=bar,
         labels=names(x), cex.axis=0.8)
    
    if (xlab == TRUE) {
      mtext(side=1, text="# clusters", line=1.0, cex=0.8)
      mtext(side=2, text="Frequency",  line=1.0, cex=0.8) 
    }
  }
  
  return(level)
}


##====================================================================
SelectCluster2 = function(label, show.plot=TRUE, xlab=TRUE,lval=1) {
  ## Determine the number of cluster and centers
  cluster.count = apply(label, 2, function(z) length(unique(z)))
  
  x = table(cluster.count)
  name  = as.numeric(names(x))
  len   = length(x)
  xxx<-order(x)
  SelectClusterNumber = function(x) {
    name  = as.numeric(names(x))
    
    x.max = max(x)
    
    level = which(x == x.max)
    if(level==1 && cluster.count[1]==1 && lval ==1 && x[1] ==max(x)){
      level<-xxx[len-1]  
    }
    if(level==1 && cluster.count[1]==1 && lval ==1 && x[1] !=max(x)){
      level<-xxx[len]  
    }
    if( cluster.count[1]!=1 && lval ==2 && (x[xxx[len-1]])/1001 > 0.10){
      level<-xxx[len-1]  
    }
    if( cluster.count[1]==1 && lval ==2 && (x[xxx[len-2]])/1001 > 0.10 && x[1] ==max(x)){
      level<-xxx[len-2]  
    }
    if( cluster.count[1]==1 && lval ==2 && (x[xxx[len-2]])/1001 > 0.10 && xxx[len-1]==1){
      level<-xxx[len-2]  
    }
    if( cluster.count[1]==1 && lval ==2 && (x[xxx[len-1]])/1001 < 0.10 && xxx[len-2]==1){
      level<-xxx[len-1]  
    }
    if( cluster.count[1]==1 && lval ==2 && (x[xxx[len-2]])/1001 < 0.10 && xxx[len-1]==1){
      level<-xxx[len]  
    }
    if( cluster.count[1]==1 && lval ==2 && (x[xxx[len-1]])/1001 < 0.10 && xxx[len-2]==1){
      level<-xxx[len-1]  
    }
    if( cluster.count[1]!=1 && lval ==2 && (x[xxx[len-1]])/1001 < 0.10){
      level<-xxx[len]  
    }
    
    level = level[length(level)]
    
    level = which(cluster.count == name[level])
    level = names(level)[1]
  }
  
  #==
  if (length(x) > 1) {
    xx = name
    xx.end = xx[-1]
    xx = xx[-len]
    id = xx.end - xx
    
    id.start = which(id <= 0)
    if (length(id.start) > 0) {
      id.end   = id.start + 1
      id.other = setdiff(seq(len), c(id.start, id.end))
      
      len.increase = length(id.start)
      if (len.increase > 1) {
        id.intermediate = 
          id.start[2:len.increase] == id.end[1:(len.increase-1)]
        
        id.intermediate = which(id.intermediate)
        
        if (length(id.intermediate) > 0) {
          id.end   = id.end[-id.intermediate]
          id.start = id.start[-(id.intermediate+1)]   
        }     
      }
      
      id.keep = id.start
      for (z in seq(length(id.start))) {
        id = id.start[z]:id.end[z]
        
        max.x = max(x[id])
        id.max = which(x[id] == max.x)
        id.max = id.max[length(id.max)]
        
        id.keep[z] = id[id.max]
        
        x[id.keep[z]] = sum(x[id])
      }
      
      id = setdiff(seq(len), c(id.keep, id.other))
      if (length(id) > 0) {
        x = x[-id]
      } 
    }
  }
  level = SelectClusterNumber(x)
  
  if (show.plot == TRUE) {
    space = 6 / length(x)
    if (space < 0.5) {
      space = 0.5
    }
    
    if (length(x) == 1) {
      bar = barplot(x/sum(x), xlab=NA, 
                    col=AddTrans("dodgerblue4", 0.45), 
                    ylab=NA,
                    cex.axis=0.8,
                    space=2,
                    width=0.1,
                    xlim=c(0,1),
                    cex=0.8,
                    xaxt="n",
                    yaxt="n")
    } else {
      bar = barplot(x/sum(x), xlab=NA, 
                    col=AddTrans("dodgerblue4", 0.45), 
                    ylab=NA,
                    cex.axis=0.8,
                    space=space,
                    width=0.1,
                    cex=0.8,
                    xaxt="n",
                    yaxt="n")
    }
    
    axis(side=2, tck=-0.015, labels=NA)
    axis(side=2, lwd=0, mgp=c(3,0.5,0), line=-0.4, 
         labels=TRUE, cex.axis=0.8)
    
    axis(side=1, tck=-0.015, labels=NA, at=bar)
    axis(side=1, lwd=0, mgp=c(3,0.5,0), line=-0.4, at=bar,
         labels=names(x), cex.axis=0.8)
    
    if (xlab == TRUE) {
      mtext(side=1, text="# clusters", line=1.0, cex=0.8)
      mtext(side=2, text="Frequency",  line=1.0, cex=0.8) 
    }
  }
  
  return(level)
}


##====================================================================
MergeCluster = function(label, V.assign, 
                        local.maxima, V.local.maxima) {
  AssignLabel = function(label, V.assign) {
    label.uniq = unique(label)
    
    cluster = V.assign
    for (lab in label.uniq) {
      id = label == lab
      id = as.numeric(names(which(id)))
      
      cluster[V.assign %in% id] = lab
    }
    return (cluster)
  }
  cluster = AssignLabel(label, V.assign)
  
  # Cluster center
  mode = local.maxima
  if (length(V.local.maxima) > 1) {
    uniq.label = unique(label)
    id = sapply(uniq.label, function(z) {
      id = names(which(label == z))
      names(which.max(V.local.maxima[id]))
    })
    #id = as.numeric(id)
    mode = local.maxima[id]
  }
  
  return (list(cluster=cluster, mode=mode))
}


##====================================================================
#' The densityCut algorithm
#'
#' @export
#' 
#' @param X A data matrix (columns are features and rows are data points)
#' @param K A integer to specify the number of neighbours in building the Knn graph.
#' Default to \eqn{K=\log_2(N)}, where N is the number of data points
#' @param knn.index An N*K data matrix for the nearest neighbour indices
#' @param knn.dist An N*K data matrix for the nearest neighbour distances
#' 
#' @param threshold A number between 0 and 1 specifying the saliency index to cut the tree. 
#' If not specified, it is selecting by stability analysis of the clustering tree
#' 
#' @param V The initial density vector of length N
#' @param D The dimensionality of data 
#' @param G A sparse Knn graph, reseaved for extension
#' 
#' @param alpha The damping factor between 0 and 1, default to 0.90
#' @param nu The saliency index in merging trees, default to \eqn{seq(0.0, 1.0, by=0.05)}
#' @param adjust Logical, whether to ajdust valley height or not
#' 
#' @param maxit The maximum number of iteration allowed in density refinement, default to 50
#' @param eps The threshold in density refinement, default to 1e-5
#' 
#' @param col A vector of clours
#' @param show.plot Logical, whether to draw clustering results
#' @param show.tip.label Logical, whether to draw the tip labels of trees
#' 
#' @param debug Logical, whether to print debug information
#' @param xlab Logical, whether to show the xlab
#' @param text subplot label
#' @param ... Reserved for extension
#' 
#' @return A list contains the clustering memberships, 
#' the modes of each cluster, and the estimated densities at each point
#' 
#' @importFrom mvtnorm rmvnorm
#' 
#' @examples
#' library(mvtnorm)
#' 
#' data(distinct.col)
#' set.seed(0)
#' 
#' N = 2^12
#' number.cluster = 64
#' N = N / number.cluster
#'   
#' i  = j = seq(-3.5, 3.5, by=1)
#' mu = expand.grid(i, j)
#' mu = as.matrix(mu)
#'
#' sigma = matrix(c(1, 0, 0, 1)*0.05, byrow=TRUE, nrow=2)
#'
#' x = lapply(seq(number.cluster), function(z) rmvnorm(N, mu[z,], sigma))
#' x = do.call(rbind, x)
#' 
#' label = lapply(1:number.cluster, function(z) rep(z, N))
#' col = AssignLabelColor(distinct.col, unlist(label))
#' NeatPlot(x, col=col, pch=4, cex=0.5)
#'
#' K = ceiling(log2(N * number.cluster))
#' a = DensityCut(X=x, K=K, alpha=0.85, nu=seq(0.0, 1.0, by=0.05), 
#'                debug=FALSE, show.plot=TRUE,  
#'                col=distinct.col)
#'
#' col = AssignLabelColor(distinct.col, a$cluster)
#' NeatPlot(x, col=col, pch=4, cex=0.5)
#'
DensityCut4MV = function(X, K, knn.index, knn.dist, V, D, G, threshold, 
                       alpha=0.90, nu=seq(0.0, 1.0, by=0.001), 
                       adjust=TRUE, maxit=500, eps=1e-5, 
                       col, show.plot=TRUE, show.tip.label=FALSE,  
                       debug=FALSE, xlab=TRUE, text=NULL, smooth = FALSE,snn=0,...) {
  
  if (missing(G)) {
    
    if (missing(X) & (missing(knn.index) | missing(knn.dist))) {
      stop("Either X or both knn.index and knn.dist should be provided!")
    }
    
    if (!missing(X)) {
      if(is.list(X)){
        
          llength<-length(X)
          alldat<-NULL
          for(i in 1:llength){
            alldat<-cbind(alldat,X[[i]])
          }
          Xtot<-alldat
        
        N = nrow(Xtot)
        D = ncol(Xtot)
      }else{
        N = nrow(X)
        D = ncol(X) 
      }
    } else if (!missing(knn.dist)) {
      N = nrow(knn.dist)
      K = ncol(knn.dist)
    }
    if (missing(D)) {
      D = 2
    }
    
    if (missing(K)) {
      K = ceiling(log2(N))
      if (K %% 2 != 0) {
        K = K + 1
      }
    }
    
    if (!missing(X) & (missing(knn.index) | missing(knn.dist))) {
      if (D > 20) {
        warning("Dimension D is greater than 20. Considering approximate knn search!")
      }
      
      ## knn search
      Xorig<-X
      if(is.list(X)){
        if(snn==0){
        print("KNN")
        knns<-list()
        KK<-K+1
        for(i in 1:llength){
          knn = FNN::get.knn(X[[i]], k=1*K, algorithm="kd_tree") 
          names(knn)[[1]]<-"idx"
          names(knn)[[2]]<-"dist"
          vec<-rep(0,N)
          vec2<-seq(1,N,1)
          knn[[1]]<-cbind(vec2,knn[[1]])
          knn[[2]]<-cbind(vec,knn[[2]])
          knns[[i]]<-knn
        }
        combknn<-merge_knnl(knns)
        names(combknn)[[1]]<-"nn.index"
        names(combknn)[[2]]<-"nn.dist"
        xa<-combknn[[1]]
        xa<-xa[,-1]
        xb<-combknn[[2]]
        xb<-xb[,-1]
        combknn[[1]]<-xa
        combknn[[2]]<-xb
        knn<-combknn
        knn.index = knn$nn.index[, 1:K]
        knn.dist  = knn$nn.dist[, 1:K]
        X<-Xtot
        }else{
          print("SNN")
          knns<-list()
         knn<-snncombl(X,K,N) 
         
         names(knn)[[1]]<-"nn.index"
         names(knn)[[2]]<-"nn.dist"
        names(knn)[[3]]<-"nn.shared"
         knn[[4]]<-knn[[2]]
         names(knn)[[4]]<-"nn.origdist"
         knn.index = knn$nn.index[, 1:K]
         knn.dist  = knn$nn.shared[, 1:K]
         if(snn==2){
           knns<-list()
           KK<-K+1
           for(i in 1:llength){
             knn = FNN::get.knn(X[[i]], k=1*K, algorithm="kd_tree") 
             names(knn)[[1]]<-"idx"
             names(knn)[[2]]<-"dist"
             vec<-rep(0,N)
             vec2<-seq(1,N,1)
             knn[[1]]<-cbind(vec2,knn[[1]])
             knn[[2]]<-cbind(vec,knn[[2]])
             knns[[i]]<-knn
           }
             combknn<-merge_knnl(knns)
             names(combknn)[[1]]<-"nn.index"
             names(combknn)[[2]]<-"nn.dist"
             xa<-combknn[[1]]
             xa<-xa[,-1]
             xb<-combknn[[2]]
             xb<-xb[,-1]
             combknn[[1]]<-xa
             combknn[[2]]<-xb
             knn<-combknn
             knn.index = knn$nn.index[, 1:K]
             knn.dist  = knn$nn.dist[, 1:K]
             X<-Xtot
             
           mm<-RKNNLIST2(X,N,K) 
           mm<-mm+1
             
         }
         for(n in 1:N){
           distvals<-rep(0,K)
           for(m in 1:K){
             if(snn==1){
             distvals[m]<-(K*m)/sum(knn$nn.shared[n,1:m])
             }
             if(snn==2){
             distvals[m]<-sum(knn.dist[n,1:m])/sum(mm[n,1:m])
             #distvals[m]<-sum(mm[n,1:m])/sum(knn.dist[n,1:m])
             #distvals[m]<-knn$nn.dist[n,m]/hh[n]
             }
             #distvals[m]<-(knn$nn.dist[n,m])/(knn$nn.shared[n,m])
           }
           knn.dist[n,]<-distvals  
         }
         knn[[2]]<-knn.dist
         combknn<-knn
         X<-Xtot
        }
      }else{
        
        print("HERE")
        if(snn==0){
        knn = FNN::get.knn(X, k=1*K, algorithm="kd_tree")
        knnorig<-knn
        knn.index = knn$nn.index[, 1:K]
        knn.dist  = knn$nn.dist[, 1:K]
        }else{
         
          if(snn==2){
          mm<-RKNNLIST2(X,N,K)
          mm<-mm+1
          knn = FNN::get.knn(X, k=1*K, algorithm="kd_tree")
          knn.index = knn$nn.index[, 1:K]
          knn.dist  = knn$nn.dist[, 1:K]
          knnorig<-knn
          }
          if(snn==1){
            knn<- sNN(X,K) 
          knn[[4]]<-knn[[1]]
          knnorig<-knn
         names(knn)[[2]]<-"nn.index"
         names(knn)[[1]]<-"nn.dist"
         names(knn)[[3]]<-"nn.shared"
         names(knn)[[4]]<-"nn.origdist"
         knn.index<-knn[[2]]
         knn.dist<-knn[[1]]
          }
         for(n in 1:N){
           distvals<-rep(0,K)
         for(m in 1:K){
           if(snn==1){
           distvals[m]<-(K*m)/sum(knn$shared[n,1:m])
           }
           if(snn==2){
           distvals[m]<-sum(knn.dist[n,1:m])/sum(mm[n,1:m])
           }
           #distvals[m]<-sum(knn$shared[n,m])/sum(knn$nn.dist[n,m])
         }
         knn.dist[n,]<-distvals  
         }
         
      }
      }
    }
    if (missing(V)) {
      if (D <= 10 & D > 0) {
        V = -D * log(knn.dist[, K])
      } else {
        V = -10 * log(knn.dist[, K])
      }
      id = is.finite(V)      
      V[!id] = max(V[id])
      if(snn==0){
      V = exp(V - LogSumExp(V))
      }
    }
    
    knn.index.col = as.vector(t(knn.index))
    knn.index.row = rep(seq(N), each=1*K)
    G = Matrix::sparseMatrix(knn.index.row, 
                             knn.index.col, 
                             x=1,
                             dims=c(N,N)) # important
    
    knn.out = lapply(seq(N), function(z) cbind(knn.index[z,], knn.dist[z,]))
    index.out = seq(N)
    
    InNeighbour = GetInNeighbour(knn.index.row, 
                                 knn.index.col, 
                                 distance=as.vector(t(knn.dist)))
    knn   = InNeighbour[[1]]
    index = InNeighbour[[2]]
  }
  
  print(G)
  ##--------------------------------------------------
  name.all = seq(N)
  diag(G)  = 1
  P        = G / (K + 1)
  
  id       = V <= .Machine$double.eps
  V[id]    = .Machine$double.eps
  V        = V / sum(V)
  method="m1"
  if(method=="m1"){
    V= EnhanceDensity(P=P, V=V, maxit=maxit, debug=debug,
                      alpha=alpha, eps=eps, smooth=FALSE)
  }
  if(method=="m2"){
    #V=knnDE2(X,X,k=10)
    V=TDA::kde(X,X,h=h)
    V= EnhanceDensity(P=P, V=V, maxit=maxit, debug=debug,
                      alpha=alpha, eps=eps, smooth=smooth)
  }
  if(method=="m3"){
    V<-V
  }
  #=--------------------------------------------------
  # Remove outlier modes before clustering
  #nearest.neighbour = NearestNeighbour(knn, index, V)
  nearest.neighbour = NearestNeighbour(knn, index, V)
  id = which(is.na(nearest.neighbour))
  mode = index[id]
  nearest.neighbour[id] = mode
  
  ## Candidate outliers
  nearest.neighbour.out = 
    NearestNeighbour(knn=knn.out, index=index.out, V=V, id=mode)
  id.outlier = !is.na(nearest.neighbour.out)
  mode.outlier = mode[id.outlier]
  
  mode.outlier.filt = -1
  if (length(mode.outlier) > 0) {
    neighbour.outlier.in = 
      lapply(id[id.outlier], function(z) knn[[z]][,1])
    neighbour.outlier.out = 
      lapply(mode.outlier, function(z) knn.out[[z]][,1])
    
    in.out.intersect = lapply(seq(length(mode.outlier)), function(z) 
      intersect(neighbour.outlier.in[[z]], 
                neighbour.outlier.out[[z]]))
    
    k.top = min(K, 2)
    id.outlier.filt = sapply(in.out.intersect, length) < K/2 | 
      sapply(mode.outlier, function(z) 
        sum(V[knn.out[[z]][1:k.top,1]] > V[z]) > 0)
    
    mode.outlier.filt = mode.outlier[id.outlier.filt]
    mode = setdiff(mode, mode.outlier.filt)
    
    nearest.neighbour[id[id.outlier][id.outlier.filt]] = 
      nearest.neighbour.out[id.outlier][id.outlier.filt]
  }
  
  #=--------------------------------------------------
  # Initial clustering
  id = AssignCluster(nearest.neighbour=nearest.neighbour, index=index, 
                     adjust=adjust, mode=mode, knn=knn, V=V, K=K, N=N)
  mode.index = id$id.mode
  local.maxima = mode[mode.index]
  valley = id$valley
  V.assign = id$V.assign
  cluster.assign = id$cluster.assign  
  
  # Assign outliers without in-neighbours, no boundary points
  id.outlier.in = which(V.assign == 0)
  id.order = order(V[id.outlier.in], decreasing=TRUE)
  id.outlier.in = id.outlier.in[id.order]
  
  V.assign[id.outlier.in] = sapply(id.outlier.in, function(z) {
    id = knn.out[[z]][, 1]
    id.sub = which(!(id %in% mode.outlier.filt))
    if (length(id.sub > 0)) {
      id = id[id.sub]
    } 
    id = id[1]
    
    V.assign[id]
  })
  
  if (is.list(V.assign)) {
    V.assign = unlist(V.assign)
  }
  print("HEREHERE")
 
  #=--------------------------------------------------
  if (!missing(threshold)) {
    if (threshold >= 1) {
      return(list(cluster=V.assign, mode=local.maxima, V=V))
    }
  }
  
  #=--------------------------------------------------
  if (show.plot == TRUE) {
    unique.label = length(unique(V.assign))
    if (missing(col)) {
      col = densitycut::distinct.col
    } 
    
    if (unique.label > length(col)) {
      col = colorRampPalette(col)(unique.label)
    }
  } else {
    col = NULL
  }
  
  #---------------------------------------------------
  # Merge clusters and select the most stable clustering
  V.local.maxima = V[local.maxima]
  mode.name = order(V[local.maxima], decreasing=TRUE)
  names(V.local.maxima) = mode.name
  names(local.maxima)   = mode.name
  
  label = MergeCut(valley=valley, 
                   cluster.assign=cluster.assign,
                   V.local.maxima=V.local.maxima, 
                   nu=nu, 
                   tip.color=col, 
                   show.tip.label=show.tip.label,
                   show.plot=show.plot, 
                   text=text,
                   xlab=xlab)
  
  if (!missing(threshold)) {
    id = which.min(abs(threshold - nu))
    level = colnames(label)[id]
  } else {
    #level = SelectCluster2(label, show.plot=show.plot, xlab=xlab,1)  
    level = SelectCluster(label, show.plot=show.plot, xlab=xlab)  
  }
  allclusts<-matrix(0,nrow=N,ncol=length(nu))
  for(i in 1:length(nu)){
    tmp =label[,i]
    names(tmp)<-mode.name
    clustx<-MergeCluster(tmp, V.assign, local.maxima, V.local.maxima)
    allclusts[,i]<-clustx$cluster
  }
  tmp  = label[, level]
  names(tmp) = mode.name
  cluster = MergeCluster(tmp, V.assign, local.maxima, V.local.maxima)
  mode = cluster$mode
  cluster = cluster$cluster
  cls<-cluster
  clsu<-unique(cls)
  clsnew<-rep(0,length(cls))
  for(i in 1:length(cls)){
    clsnew[i]<-which(clsu==cls[i])
  }
  cluster<-clsnew
  print("HEREHEREHERE")
  if(is.list(Xorig)){
    if(snn>=1){
  return(list(cluster=cluster, mode=mode, V=V,label=label,level=level,knns=knns,combknnsnn=combknn,allclusts=allclusts))
    }else{
      return(list(cluster=cluster, mode=mode, V=V,label=label,level=level,knns=knns,combknn=combknn,allclusts=allclusts))   
    }
  }else{
    if(snn>=1){
    return(list(cluster=cluster, mode=mode, V=V,label=label,level=level,knnssnn=knnorig,allclusts=allclusts)) 
    }else{
      return(list(cluster=cluster, mode=mode, V=V,label=label,level=level,knns=knnorig,allclusts=allclusts))   
    }
  }
}


#=====================================================================
PlotDensitycut = function(label, base, show.tip.label=TRUE, 
                          tip.color, xlab=TRUE) {
  ## For a subtree, start from the roots and traverse to the tips 
  TraverseSubTree = function(x, i, add.internal.node=FALSE, 
                             prob, prev.prob) {
    # Current distance
    if (missing(prob)) {
      prob = nu.max - merge.point[i]
    }
    
    # Previous splitting distance
    if (missing(prev.prob)) {
      prev.prob = single.root
    }
    
    # Internal nodes
    if (i < (ncol(label))) {    
      id = which(as.character(label[,i]) == x)
      next.level.node = as.character(unique(label[id, i+1])) 
      desc = NULL
      
      # No divergence - proceed to the next level
      if (length(next.level.node) == 1) { 
        newickout = TraverseSubTree(x=next.level.node, 
                                    i=i+1, 
                                    add.internal.node, 
                                    prob=prob, 
                                    prev.prob=prev.prob)
      } else {
        prob = merge.point[i+1] - prev.prob
        prev.prob = prob + prev.prob
        
        for (x.i in next.level.node) {
          subtree = TraverseSubTree(x.i, i+1, 
                                    add.internal.node, 
                                    prev.prob=prev.prob)
          desc = c(desc, paste(subtree, sep=""))
        }
        
        internal.node = NULL
        if (add.internal.node == TRUE) internal.node = x
        
        sub.tree  = paste(desc, collapse=",")
        newickout = paste("(", sub.tree, "):", 
                          prob, internal.node, sep="")
      }
    } else {
      newickout = x
      newickout = paste(newickout, prob, sep=":")
    }
  }
  
  merge.point = as.numeric(colnames(label))
  single.root = merge.point[1]
  nu.max = merge.point[length(merge.point)]
  
  ConvertClusterNewick = function(x, add.internal.node=FALSE) {
    newick = NULL
    
    # For each subtree
    level.one = as.character(unique(x[,1]))
    for(x in level.one) {
      subtree.x = TraverseSubTree(x, 1, 
                                  add.internal.node=FALSE, 
                                  prob=nu.max - single.root)
      newick = c(newick, subtree.x)
    }
    newick = paste(newick, collapse=",")
    newick = paste("(", newick, "):", single.root, ";", sep="")
    
    return(newick)
  }
  
  # (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
  out  = ConvertClusterNewick(label, add.internal.node=FALSE)
  tree = ape::read.tree(text=out)
  tree = ape::collapse.singles(tree)
  
  id.col = match(tree$tip.label, unique(label[, ncol(label)]))
  tip.color = tip.color[unique(label[, ncol(label)])][id.col]
  
  ape::plot.phylo(tree, type="phylogram", 
                  label.offset=0.01, 
                  xaxt="n", 
                  show.tip.label=show.tip.label,
                  tip.color=tip.color,
                  ann=FALSE,
                  root.edge=TRUE,
                  no.margin=FALSE, 
                  x.lim=c(base, nu.max))
  ape::tiplabels(bg=tip.color, pch=21, col="ivory4", cex=1) 
  
  line.space = 0.10
  label.tick = seq(0, nu.max, by=nu.max/5)
  axis(side=1, tck=-0.015, labels=NA, at=label.tick, line=line.space)
  axis(side=1, lwd=0, mgp=c(3,0.5,0), line=-0.4,
       at=label.tick, labels=label.tick, cex.axis=0.8)
  
  if (xlab == TRUE) {
    mtext(side=1, text="saliency index", line=1.0, cex=0.8)
  }
}
##====================================================================
#' Efficient Knn-search in high-dimensional search (upto 1000 dimensions)
#' 
#' export
#' 
#' @param X A data matrix (columns are features and rows are data points)
#' @param num.tree The number of trees for random projection
#' @param K A integer to specify the number of neighbours in building the Knn graph.
#' Default to \eqn{K=\log_2(N)}, where N is the number of data points
#' 
#' @return A list containg knn.index, knn.dist, the dimensionality D, 
#' and the projected trees: tree.index
#' 

GetKnnRandomProjection = function(X, num.tree=50, K) {
  X = as.matrix(X)
  D = ncol(X)
  tree.index = new(RcppAnnoy::AnnoyEuclidean, D)
  N = nrow(X)
  
  for (i in seq(N)) {
    tree.index$addItem(i-1, X[i, ])
  }
  tree.index$build(num.tree)
  
  #====================
  if (missing(K)) {
    K  = ceiling(log2(N))
  }
  knn.index = matrix(0, nrow=N, ncol=K)
  knn.dist  = knn.index
  for (i in seq(N)) {
    idx = tree.index$getNNsByItem(i-1, K+1)
    idx = idx[-1]
    
    knn.dist[i, ]  =  sapply(idx, function(z) 
      tree.index$getDistance(i-1, z))
    
    knn.index[i, ] = idx + 1
  }
  
  return(list(knn.index=knn.index, 
              knn.dist =knn.dist, 
              D=D, 
              tree.index=tree.index)
  )
}
#' @useDynLib densitycut sexp_log_sum_exp
LogSumExp = function(x) { 
  # The log-sum-exp trick, x is a vector. 
  # A matrix input will be converted to a vector
  
  x = as.numeric(x)
  x = na.omit(x)
  if (length(x) < 1) {
    stop("Error: x should be a vector of length greater than zero!")
  }
  
  tmp = .Call("sexp_log_sum_exp", x)
  
  return(tmp)
}


#======================================================================
#' Reset the default parameters for the plot function
#' 
#' @export
#' @param ... go to plot
#' @param xlab A title for the x axis 
#' @param ylab A title for the y axis
#' @param xaxt see par 
#' @param yaxt see par
#' @param xtck.label Logical, wheter to draw the x axis tick labels
#' @param ytck.label Logical, wheter to draw the y axis tick labels
#' @param cex.axis cex for axis 
#' 
#' @importFrom graphics plot axis mtext 
#' 
NeatPlot = function(..., xlab, ylab, xaxt="s", yaxt="s", 
                    xtck.label=TRUE, ytck.label=TRUE, cex.axis=0.8) {
  
  plot(xaxt="n", yaxt="n", xlab=NA, ylab=NA, ...) 
  
  if (xaxt == "s") {
    axis(side=1, tck=-0.015, labels=NA)
    axis(side=1, lwd=0, mgp=c(3,0.5,0), line=-0.4,
         labels=xtck.label, cex.axis=cex.axis)
  }
  if (yaxt == "s") {
    axis(side=2, tck=-0.015, labels=NA)
    axis(side=2, lwd=0, mgp=c(3,0.5,0), line=-0.4,
         labels=ytck.label, cex.axis=cex.axis)
  }
  
  if (!missing(xlab)) {
    mtext(side=1, text=xlab, cex=cex.axis, line=1.0)
  }
  if (!missing(ylab)) {
    mtext(side=2, text=ylab, cex=cex.axis, line=1.0)
  }
}


#======================================================================
#' Calcular the normalized mutual information
#' 
#' @export
#' 
#' @param ground.truth The cluster ground truth
#' @param cluster The cluster labels
#' 
ComputeNMI = function(ground.truth, cluster) {
  N = length(ground.truth)
  if (length(cluster) != N) {
    stop("The two vectors should be the same length!")
  }
  
  x = table(ground.truth, cluster)
  x = x / sum(x)
  
  p.x = rowSums(x)
  p.y = colSums(x)
  y = outer(p.x, p.y)
  
  i = which((x > 0) & (y > 0))
  mutual.info = sum(x[i] * log2(x[i] / y[i]))
  entropy.x = sum(p.x * log2(ifelse(p.x > 0, p.x, 1)))
  entropy.y = sum(p.y * log2(ifelse(p.y > 0, p.y, 1)))
  
  nmi = mutual.info / sqrt(entropy.x * entropy.y)
  
  return(nmi)
}


# Calculate the entropy
ComputeEntropy = function(x) {
  N  = length(x)
  prob = table(x) / N
  
  entropy = -sum(prob * log2(prob))
  
  return(entropy)
}


#======================================================================
#' Assign colors to labels
#' 
#' @author Jiarui Ding
#' 
#' @export
#' 
#' @param col A set of colors
#' @param label Cluster labels
#' @param uniq.label The initial labels
#' 
AssignLabelColor = function(col, label, uniq.label) {
  if (is.list(label)) {
    label = unlist(label)
  }
  
  if (missing(uniq.label)) {
    uniq.label = unique(label)
  }
  len = length(uniq.label)
  uniq.label.ind = seq(len)
  
  if (len > length(col)) {
    col = colorRampPalette(col)(len)
  }
  
  col.label = label
  for(z in uniq.label.ind) {
    id = label == uniq.label[z]
    col.label[id] = col[z]
  }
  
  return(col.label)
}


#======================================================================
#' Add transpancy to color
#' 
#' @author Jiarui Ding
#' 
#' @export
#' 
#' @param col A set of colors
#' @param trans A vector or vector of the same length as col, 
#' 0 being fully visable and 1 being fully transpant
#' 
AddTrans = function(col, trans) {
  # Add transcripancy to colours. 
  # Works with either color and trans a vector 
  # of equal length, or one of the two of length 1.
  
  # col  :
  #      colour vector
  # trans:
  #      0 being fully visable and 1 being fully transpant
  # 
  if (any(trans > 1) | any(trans < 0)) {
    stop("Error: trans should in the range of [0, 1]") 
  }
  if (length(col) != length(trans) & 
      !any(c(length(col), length(trans))==1)) 
    stop("Error: vector lengths doesn't match")
  if (length(col)==1 & length(trans)>1)
    color = rep(color,length(trans))
  if (length(trans)==1 & length(col)>1) 
    trans = rep(trans, length(col))
  
  num2hex = function(x) {
    hex = unlist(strsplit("0123456789ABCDEF", split=""))
    return (paste(hex[(x-x%%16)/16 + 1], hex[x%%16 + 1], sep=""))
  }
  rgb = rbind(col2rgb(col), round(255 - trans*255))
  col = apply(apply(rgb, 2, num2hex), 2, paste, collapse="")
  col = paste("#", col, sep="")
  
  return(col)
}

##=====================================================================
