snncombl<-function(X,K,N){
 xlength<-length(X) 
 x<-list()
 for(i in 1:xlength){
  x[[i]]<-sNN(X[[i]],K) 
 }
 outind<-matrix(0,nrow=N,ncol=K)
 outdist<-matrix(0,nrow=N,ncol=K)
 outshared<-matrix(0,nrow=N,ncol=K)

 for(i in 1:N){
   sharedall<-NULL
   indall<-NULL
   distall<-NULL
   for(j in 1:xlength){
   sharedall<-c(sharedall,x[[j]]$shared[i,])  
   indall<-c(indall,x[[j]]$id[i,])
   distall<-c(distall,x[[j]]$dist[i,])
   }
   ord<-rev(order(sharedall))
   #ord<-order(distall)
   outshared[i,]<-sharedall[ord[1:K]]
   outind[i,]<-indall[ord[1:K]]
   outdist[i,]<-distall[ord[1:K]]
   #od<-order(outdist[i,])
   #outdist[i,]<-outdist[i,od]
   #outind[i,]<-outind[i,od]
   
   
 }
 return(list(outind,outdist,outshared))
}