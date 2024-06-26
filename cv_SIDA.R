cv_SIDA=function(Xdata=Xdata,Y=Y,withCov=FALSE,plotIt=FALSE,Xtestdata=NULL,Ytest=NULL,isParallel=TRUE,ncores=NULL,gridMethod='RandomSearch',AssignClassMethod='Joint',nfolds=5,ngrid=8,standardize=TRUE,maxiteration=20, weight=0.5,thresh=1e-03,seedval=1234,lambda=NULL){

  #check inputs for training data
  dsizes=lapply(Xdata, function(x) dim(x))
  n=dsizes[[1]][1]
  nsizes=lapply(Xdata, function(x) dim(x)[1])

  if(all(nsizes!=nsizes[[1]])){
    stop('The datasets  have different number of observations')
  }


  #check data
  if (is.list(Xdata)) {
    D = length(Xdata)
    if(D==1){
      stop("There should be at least two datasets")
    }
  } else {
    stop("Input data should be a list")
  }


  #set defaults
  #If testing data are not provided, the default is to use training data
  if(is.null(Xtestdata)){
    Xtestdata=Xdata
  }

  #check inputs for testing data
  ntestsizes=lapply(Xtestdata, function(x) dim(x)[1])
  if(all(ntestsizes!=ntestsizes[[1]])){
    stop('The testing datasets  have different number of observations')
  }

  if(is.null(Ytest)){
    Ytest=Y
  }

  if(is.null(withCov)){
    withCov=FALSE
  }

  if(is.null(plotIt)){
    plotIt=FALSE
  }

  if(is.null(standardize)){
    standardize=TRUE
  }


  #standardize if true
  if(standardize==TRUE){
    Xdata=lapply(Xdata,function(x)scale(x,center=TRUE,scale=TRUE))
    Xtestdata=lapply(Xtestdata,function(x)scale(x,center=TRUE,scale=TRUE))
  }

  if(is.null(gridMethod)){
    gridMethod='RandomSearch'
  }

  if(is.null(AssignClassMethod)){
    AssignClassMethod='Joint'
  }

  if(is.null(isParallel)){
    isParallel=TRUE
  }

  if(is.null(nfolds)){
    nfolds=5
  }

  if(is.null(ngrid)){
    ngrid=8
  }


  if(is.null(maxiteration)){
    maxiteration=20
  }

  if(is.null(weight)){
    weight=0.5
  }

  if(is.null(thresh)){
    thresh=1e-03
  }

  set.seed(seedval)
  nK=length(unique(as.vector(Y))) -1

  # nc=length(unique(as.vector(Y)))
  # Nn=mat.or.vec(nc,1)
  # foldid=list()
  # for(i in 1:nc)
  # {
  #   Nn[i]=sum(Y==i)
  #   mod1=Nn[i]%%nfolds
  #   if(mod1==0){
  #     foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds)),Nn[i])
  #   }else if(mod1> 0){
  #     foldid[[i]]=sample(c(rep(1:nfolds,times=floor(Nn[i])/nfolds), 1:(Nn[i]%%nfolds)),Nn[i])
  #   }
  # }
  # foldid=unlist(foldid)
  foldid=createFolds(Y,nfolds,list=F)

  #obtain tuning range common to all K


  if(withCov==TRUE){
    Dnew=D-1
  }else if(withCov==FALSE){
    Dnew=D
  }

  if(Dnew>2){
    ngrid=5
  }

  ##modification scan full lambda range
  ##myTauvec=sidatunerange(Xdata,Y,ngrid,standardize,weight,withCov)
  myTauvec<-list(Tauvec=list())
  myTauvec$Tauvec[1:D]<-list(lambda)
  ###############################

  #define the grid
  mygrid=expand.grid(do.call(cbind,myTauvec))
  #mygrid=data.frame(matrix(rep(lambda,each=D),ncol=D,byrow=T))
  gridcomb=dim(mygrid)[1]
  if(gridMethod=='RandomSearch'){
    if(D==2){
      ntrials=floor(0.2*gridcomb)}
    else if(D>2){
      ntrials=floor(0.15*gridcomb)
    }
    mytune=sample(1:gridcomb, ntrials, replace = FALSE)
    gridValues=mygrid[mytune,]
  }else if(gridMethod=='GridSearch'){
    gridValues=mygrid
  }



  CVOut=matrix(0, nfolds, nrow(gridValues))
  #cross validation
  if(isParallel==TRUE){
      registerDoParallel()
      if(is.null(ncores)){
       ncores=parallel::detectCores()
       ncores=ceiling(ncores/2)}
       cl=makeCluster(ncores)
       registerDoParallel(cl)
    CVOut=matrix(0, nrow(gridValues), nfolds)
    mycv=foreach(i = 1:nrow(gridValues), .combine='rbind',.export=c('sida','sidainner','myfastinner','myfastIDAnonsparse','mysqrtminv','sidaclassify', 'sidatunerange','DiscriminantPlots','CorrelationPlots'),.packages=c('CVXR','RSpectra')) %dopar% {
      mTau=sapply(1:D, function(itau) list(t(gridValues[,itau][i])))
      #cat("Begin CV-fold", i, "\n")
      CVOut[i,]= sapply(1:nfolds, function(j){
        testInd=which(foldid==j)
        testX=lapply(Xdata, function(x) x[testInd,])
        testY=Y[testInd]
        trainX=lapply(Xdata, function(x) x[-testInd,])
        trainY=Y[-testInd]
        mysida=sida(trainX,trainY,mTau,withCov,Xtestdata=testX,testY,AssignClassMethod='Joint',plotIt=FALSE,standardize,maxiteration,weight,thresh)
        return(paste0(min(mysida$sidaerror),"_",paste(mysida$PredictedClass,collapse="")))
      } )
    }
    CVOut=t(mycv)
    stopCluster(cl)
  }else if(isParallel==FALSE){
    CVOut=matrix(0, nfolds, nrow(gridValues))
    for (i in 1:nfolds){
      testInd=which(foldid==i)
      testX=lapply(Xdata, function(x) x[testInd,])
      testY=Y[testInd]
      trainX=lapply(Xdata, function(x) x[-testInd,])
      trainY=Y[-testInd]

      #cat("Begin CV-fold", i, "\n")
      CVOut[i,]= sapply(1:nrow(gridValues), function(itau){
        mTau=sapply(1:D, function(d) list(t(gridValues[itau,][d])))
	sink('/dev/null')
        mysida=sida(trainX,trainY,mTau,withCov,Xtestdata=testX,testY,AssignClassMethod='Joint',plotIt=FALSE,standardize,maxiteration,weight,thresh)
	sink()
        return(paste0(min(mysida$sidaerror),"_",paste(mysida$PredictedClass,collapse="")))
      } )
    }
  }
  CVtot<-apply(CVOut,1,function(x) unlist(strsplit(x,"_")))
  CVOut<-t(apply(CVtot[seq(1,nrow(CVtot),2),],c(1,2),as.numeric))
  PredictedClass<-CVtot[seq(2,nrow(CVtot),2),]

  #compute average classification error
  minEorrInd=max(which(colMeans(CVOut)==min(colMeans(CVOut))))
  optTau=gridValues[ minEorrInd,]

  PredictedClass<-as.numeric(unlist(strsplit(paste(PredictedClass[minEorrInd,],collapse = ""),"")))[order(order(foldid))]

  #Apply on testing data
  moptTau=sapply(1:D, function(i) list(t(gridValues[minEorrInd,][i])))
  sink('/dev/null')
  mysida=sida(Xdata=Xdata,Y=Y,Tau=moptTau,withCov,Xtestdata=Xtestdata,Ytest=Ytest,AssignClassMethod='Joint',plotIt=FALSE,standardize,maxiteration,weight,thresh)
  sink()

  ss=list()
  #sum pairwise RV coefficients
  for(d in 1:D){
    dd=setdiff(seq(1, D, by= 1),d)
    #correlations
    sumCorr2=0;
    for (jj in 1:length(dd)){
      j=dd[jj];
      X1=Xtestdata[[d]]%*%mysida$hatalpha[[d]]
      X2=Xtestdata[[j]]%*%mysida$hatalpha[[j]]
      X1=scale(X1, center=TRUE,scale=FALSE)
      X2=scale(X2, center=TRUE,scale=FALSE)
      X1X2=t(X1)%*%X2/dim(X1)[1]
      X1X1=t(X1)%*%X1/dim(X1)[1]
      X2X2=t(X2)%*%X2/dim(X2)[1]
      sumcorr3=sum(diag(X1X2%*%t(X1X2)))/(sqrt(sum(diag(X1X1%*%X1X1)))*sqrt(sum(diag(X2X2%*%X2X2))))
      sumCorr2=sumCorr2+sumcorr3
    }
    ss[[d]]=sumCorr2/length(dd)
  }

  sidacorrelation=sum(do.call(rbind,ss))/D

  #Produce discriminant and correlation plot if plotIt=T
  if(plotIt==TRUE){
  DiscriminantPlots(Xtestdata,Ytest,mysida$hatalpha)
   CorrelationPlots(Xtestdata,Ytest,mysida$hatalpha)
  }else{
    myDiscPlot=NULL
    myCorrPlot=NULL
  }
  result=list(CVOut=CVOut,sidaerror=mysida$sidaerror,sidacorrelation=sidacorrelation,hatalpha=mysida$hatalpha,PredictedClass=PredictedClass, optTau=moptTau,gridValues=gridValues, AssignClassMethod=AssignClassMethod, gridMethod=gridMethod,foldid=foldid,mysida=mysida)
  return(result)
}
