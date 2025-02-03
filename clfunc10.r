clproc<-function(x,y){
  cls<-list()
  nu<-seq(0,1.0,by=0.001) 
  lnu<-length(nu)
  lx<-length(x)
  indval<-list()
  indnum<-list()
  for(j in 1:lx){
    clstr<-x[[j]]
    clallvals<-clstr$allclusts
    cls[[j]]<-clallvals
    lvec<-apply(clstr$label, 2, function(z) length(unique(z)))
    ll<-lvec
    names(ll)<-NULL
    lvu<-unique(ll)
    lvul<-length(lvu)
    stabmat<-matrix(0,nrow=lvul,ncol=2)
    for(i in 1:lvul){
      
      stabmat[i,1]<-lvu[i]
      stabmat[i,2]<-length(lvec[ll==lvu[i]])/(lnu+1)
    }
    #cls2<-clstr$cluster
    #clsu<-unique(cls2)
    #clsnew<-rep(0,length(cls2))
    #for(i in 1:length(cls2)){
    # clsnew[i]<-which(clsu==cls2[i])
    #}
    maxval<-which(stabmat[,2]==max(stabmat[,2]))
    clustnum<-stabmat[maxval,1]
    indval[[j]]<-stabmat[,1]
    indnum[[j]]<-stabmat[,2]
  }
  rr<-Reduce(union,indval)
  totclustnums<-rep(0,length(rr))
  for(i in 1:lx){
    wcl<-match(indval[[i]],rr)
    totclustnums[wcl]<-totclustnums[wcl]+indnum[[i]]
    
  }
  tno<-order(totclustnums)
  tnoo<-rr[order(totclustnums)]
  tn<-rr[totclustnums==max(totclustnums)]
  if(tn==1){
    tn=tnoo[length(tno)-1]
    tn2=tnoo[length(tno)-2]
  }else{
    tn<-tn
    tn2=tnoo[length(tno)-1]
  }
  for(i in 1:lx){
    clall<-cls[[i]]
    clnumsall<-rep(0,lnu)
    for(j in 1:lnu){
      clnumsall[j]<-length(unique(clall[,j]))
    }
  }
  clnumb<-which(clnumsall==tn)
  if(length(clnumb)==0){
    clnumb<-which(clnumsall==tn2) 
  }
  numhits<-length(clnumb)
  clnumb<-clnumb[1]
  cls<-clnumb
  clsu<-unique(clall[,clnumb])
  clsnew<-rep(0,length(clall[,i]))
  allcols<-rainbow(length(clsu))
  for(j in 1:length(clsu)){
    coinc<-which(clall[,clnumb]==clsu[j])
    clsnew[coinc]<-j
    #clsnew[which(clall[,clnumb]==clsu[j])]<-clsu[j]
  }
  #cols<-allcols[clsnew]
  oc<-order(rr)
  return(list(allcols,clsnew,rr[oc],totclustnums[oc]))  
}
clproc2x<-function(x,alphas){
  peakweight<-function(x,y){
    if(x[1]==1){
      y[1]<-0  
    }
    rat<-max(y)/sum(y)
    val<-which(y==max(y))
    return(list(rat,val))
  }
  weightvals<-rep(0,3*length(alphas))
  weightnums<-rep(0,3*length(alphas))
  cls<-list()
  nu<-seq(0,1.0,by=0.001) 
  lnu<-length(nu)
  lx<-length(x)
  indval<-list()
  indnum<-list()
  alphaclusts<-matrix(0,nrow=length(nu),ncol=lx)
  for(j in 1:lx){
    clstr<-x[[j]]
    clallvals<-clstr[[1]]$allclusts
    cls[[j]]<-clallvals
    lvec<-apply(clstr[[1]]$label, 2, function(z) length(unique(z)))
    alphaclusts[,j]<-lvec
    ll<-lvec
    names(ll)<-NULL
    lvu<-unique(ll)
    lvul<-length(lvu)
    stabmat<-matrix(0,nrow=lvul,ncol=2)
    for(i in 1:lvul){
      
      stabmat[i,1]<-lvu[i]
      stabmat[i,2]<-length(lvec[ll==lvu[i]])/(lnu+1)
    }
    #cls2<-clstr$cluster
    #clsu<-unique(cls2)
    #clsnew<-rep(0,length(cls2))
    #for(i in 1:length(cls2)){
    # clsnew[i]<-which(clsu==cls2[i])
    #}
    maxval<-which(stabmat[,2]==max(stabmat[,2]))
    clustnum<-stabmat[maxval,1]
    indval[[j]]<-stabmat[,1]
    indnum[[j]]<-stabmat[,2]
    pw<-peakweight(indval[[j]],indnum[[j]])
    weightvals[j]<-pw[[1]]
    ij<-indval[[j]]
    weightnums[j]<-ij[pw[[2]]]
    
  }
  
  output<-list(weightvals,weightnums)
  return(output)  
}
clproc3x<-function(x,alphas,weights){
  cls<-list()
  nu<-seq(0,1.0,by=0.001) 
  lnu<-length(nu)
  lx<-length(x)
  indval<-list()
  indnum<-list()
  alphaclusts<-matrix(0,nrow=length(nu),ncol=lx)
  for(j in 1:lx){
    clstr<-x[[j]]
    clallvals<-clstr[[1]]$allclusts
    cls[[j]]<-clallvals
    lvec<-apply(clstr[[1]]$label, 2, function(z) length(unique(z)))
    alphaclusts[,j]<-lvec
    ll<-lvec
    names(ll)<-NULL
    lvu<-unique(ll)
    lvul<-length(lvu)
    stabmat<-matrix(0,nrow=lvul,ncol=2)
    for(i in 1:lvul){
      
      stabmat[i,1]<-lvu[i]
      stabmat[i,2]<-length(lvec[ll==lvu[i]])/(lnu+1)
    }
    #cls2<-clstr$cluster
    #clsu<-unique(cls2)
    #clsnew<-rep(0,length(cls2))
    #for(i in 1:length(cls2)){
    # clsnew[i]<-which(clsu==cls2[i])
    #}
    maxval<-which(stabmat[,2]==max(stabmat[,2]))
    clustnum<-stabmat[maxval,1]
    indval[[j]]<-stabmat[,1]
    indnum[[j]]<-stabmat[,2]
  }
  r1<-seq(1,length(alphas)*3,3)
  r2<-seq(2,length(alphas)*3,3)
  r3<-seq(3,length(alphas)*3,3)
  rr<-Reduce(union,indval)
  rr1<-Reduce(union,indval[r1])
  rr2<-Reduce(union,indval[r2])
  rr3<-Reduce(union,indval[r3])
  totclustnums<-rep(0,length(rr))
  for(i in 1:lx){
    wcl<-match(indval[[i]],rr)
    totclustnums[wcl]<-totclustnums[wcl]+indnum[[i]]*weights[i]
    
  }
  totclustnums1<-rep(0,length(rr1))
  for(i in r1){
    wcl<-match(indval[[i]],rr1)
    totclustnums1[wcl]<-totclustnums1[wcl]+indnum[[i]]*weights[i]
    
  }
  totclustnums2<-rep(0,length(rr2))
  for(i in r2){
    wcl<-match(indval[[i]],rr2)
    totclustnums2[wcl]<-totclustnums2[wcl]+indnum[[i]]*weights[i]
    
  }
  totclustnums3<-rep(0,length(rr3))
  for(i in r3){
    wcl<-match(indval[[i]],rr3)
    totclustnums3[wcl]<-totclustnums3[wcl]+indnum[[i]]*weights[i]
    
  }
  tno<-order(totclustnums)
  tnoo<-rr[order(totclustnums)]
  tn<-rr[totclustnums==max(totclustnums)]
  tno1<-order(totclustnums1)
  tnoo1<-rr1[order(totclustnums1)]
  tnr1<-rr1[totclustnums1==max(totclustnums1)]
  tno2<-order(totclustnums2)
  tnoo2<-rr2[order(totclustnums2)]
  tnr2<-rr2[totclustnums2==max(totclustnums2)]
  tno3<-order(totclustnums3)
  tnoo3<-rr3[order(totclustnums3)]
  tnr3<-rr3[totclustnums3==max(totclustnums3)]
  if(tn[1]==1){
    tn=tnoo[length(tno)-1]
    tn2=tnoo[length(tno)-2]
  }else{
    tn<-tn
    tn2=tnoo[length(tno)-1]
  }
  if(tnr1[1]==1){
    tnr1=tnoo1[length(tno1)-1]
    tnr12=tnoo1[length(tno1)-2]
  }else{
    tnr1<-tnr1
    tnr12=tnoo1[length(tno1)-1]
  }
  if(tnr2[1]==1){
    tnr2=tnoo2[length(tno2)-1]
    tnr22=tnoo2[length(tno2)-2]
  }else{
    tnr2<-tnr2
    tnr22=tnoo2[length(tno2)-1]
  }
  if(tnr3[1]==1){
    tnr3=tnoo3[length(tno3)-1]
    tnr32=tnoo3[length(tno3)-2]
  }else{
    tnr2<-tnr2
    tnr32=tnoo2[length(tno2)-1]
  }
  
  clnumsallj<-matrix(0,lnu,lx)
  #CLS IS CLUSTERING (N by NU) FOR EACH ALPHA VALUE
  for(i in 1:lx){
    clall<-cls[[i]]
    clnumsall<-rep(0,lnu)
    for(j in 1:lnu){
      clnumsall[j]<-length(unique(clall[,j]))
      
                           
    }
    clnumsallj[,i]<-clnumsall
  }
  #FINDS NUMBER OF UNIQUE CLUSTERS AT EACH NU VALUE ACROSS ALL ALPHAS
  clustfreq<-apply(clnumsallj,1,function(x)length(unique(x)))
  
  for(i in r1){
    clallr1<-cls[[i]]
    clnumsallr1<-rep(0,lnu)
    for(j in 1:lnu){
      clnumsallr1[j]<-length(unique(clallr1[,j]))
    }
  }
  for(i in r2){
    clallr2<-cls[[i]]
    clnumsallr2<-rep(0,lnu)
    for(j in 1:lnu){
      clnumsallr2[j]<-length(unique(clallr2[,j]))
    }
  }
  for(i in r3){
    clallr3<-cls[[i]]
    clnumsallr3<-rep(0,lnu)
    for(j in 1:lnu){
      clnumsallr3[j]<-length(unique(clallr3[,j]))
    }
  }
  clnumb<-which(clustfreq==tn)
  if(length(clnumb)==0){
    clnumb<-which(clustfreq==tn2) 
  }
  clustmatch<-list()
  # CLUSTMATCHNUMS IS A VECTOR OF LENGTH ALPHAS WHICH CONTAINS THE NUMBER OF NU VALUES AT WHICH THE CLUSTER NUMBER IS THE "BEST" CLUSTER NUMBER (TN)
   clustmatchnums<-apply(clnumsallj,2,function(x) length(which(x==tn)))
   #CLUSTMATCH IS A LIST FOR EACH ALPHA VALUE OF THE NU INDICES WHERE THE CLUSTER NUMBER IS TN
   clustmatch<-apply(clnumsallj,2,function(x) which(x==tn))
   #MCL IS THE ALPHA VALUE WHERE CLUSTMATCH IS MAXIMAL
   mcl<-which(clustmatchnums==max(clustmatchnums))
   if(length(mcl)>1){
     mcl<-mcl[1]
   }
  
  clnumbr1<-which(clnumsallr1==tnr1)
  if(length(clnumbr1)==0){
    clnumbr1<-which(clnumsallr1==tnr12) 
  }
  clnumbr2<-which(clnumsallr2==tnr2)
  if(length(clnumbr2)==0){
    clnumbr2<-which(clnumsallr2==tnr22) 
  }
  clnumbr3<-which(clnumsallr3==tnr3)
  if(length(clnumbr3)==0){
    clnumbr3<-which(clnumsallr3==tnr32) 
  }
  numhits<-length(clnumb)
  #CLNUMB IS THE LST OF NU INDICES WHERE THE CLUSTER NUMBER IS TN AT ALPHA VALUE MCL
  clnumb<-clustmatch[[mcl]]
  #TAKE JUST FIRST VALUE
  clnumb<-clnumb[1]
  #cls<-clnumb
  #CALLL IS THE "ALLCLUST" MATRIX AT ALPHA VALUE MCL
  clall<-cls[[mcl]]
  #CLSU FINGS THE UNIQUE CLUSTER NUMBERS FOR THIS CLUSTERING
  clsu<-unique(clall[,clnumb])
  clsnew<-rep(0,length(clall[,1]))
  allcols<-rainbow(length(clsu))
  #CLSNEW ASSIGNS A CLUSTER COLOUR TO EACH SUBJECT BASED TN CLUSTERS AT ALPHA VALUE MCL 
  #THIS CORRESPONDS TO THE HIGHEST POINT ON THE ALPHA WEIGHTS GRAPH
  for(j in 1:length(clsu)){
    coinc<-which(clall[,clnumb]==clsu[j])
    clsnew[coinc]<-j
    #clsnew[which(clall[,clnumb]==clsu[j])]<-clsu[j]
  }
  numhitsr1<-length(clnumbr1)
  clnumbr1<-clnumbr1[1]
  clsr1<-clnumbr1
  clsur1<-unique(clallr1[,clnumbr1])
  clsnewr1<-rep(0,length(clall[,1]))
  allcols<-rainbow(length(clsu))
  for(j in 1:length(clsur1)){
    coincr1<-which(clallr1[,clnumbr1]==clsur1[j])
    clsnewr1[coincr1]<-j
    #clsnew[which(clall[,clnumb]==clsu[j])]<-clsu[j]
  }
  numhitsr2<-length(clnumbr2)
  clnumbr2<-clnumbr2[1]
  clsr2<-clnumbr2
  clsur2<-unique(clallr2[,clnumbr2])
  clsnewr2<-rep(0,length(clall[,1]))
  allcols<-rainbow(length(clsu))
  for(j in 1:length(clsur2)){
    coincr2<-which(clallr2[,clnumbr2]==clsur2[j])
    clsnewr2[coincr2]<-j
    #clsnew[which(clall[,clnumb]==clsu[j])]<-clsu[j]
  }
  numhitsr3<-length(clnumbr3)
  clnumbr3<-clnumbr3[1]
  clsr3<-clnumbr3
  clsur3<-unique(clallr3[,clnumbr3])
  clsnewr3<-rep(0,length(clall[,1]))
  allcols<-rainbow(length(clsu))
  for(j in 1:length(clsur3)){
    coincr3<-which(clallr3[,clnumbr3]==clsur3[j])
    clsnewr3[coincr3]<-j
    #clsnew[which(clall[,clnumb]==clsu[j])]<-clsu[j]
  }
  #cols<-allcols[clsnew]
  oc<-order(rr)
  ocr1<-order(rr1)
  ocr2<-order(rr2)
  ocr3<-order(rr3)
  output<-list(allcols,clsnew,rr[oc],totclustnums[oc],clsnewr1,rr1[ocr1],totclustnums1[ocr1],clsnewr2,rr2[ocr2],totclustnums2[ocr2],clsnewr3,rr3[ocr3],totclustnums3[ocr3],clustmatchnums)
  return(output)  
}