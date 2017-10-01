getNNG <- function (X, knng) {
  # Purpose: get nearest-neighbour graph (NNG) of X.
  # Input: X, matrix of double -- each row is the feature vector of a data point.
  #        k, int -- k nearest neighbours.
  # Output: nng, igraph -- NNG of X.
  
  library(igraph)
  n <- nrow(X)
  # create adjacency matrix for X: edges weighted by dist(X).
  adjMat <- as.matrix(dist(X)) 
  # find NN.
  adjMat <- adjMat+diag(n)*1.797693e+308 # assign large number to diagonal.
  indkNN <- matrix(apply(adjMat,1,order)[1:knng,],knng,n) # index of kNN: order index of each col from nearest neighbour to furthest, and extract the first k.
  # create adjacency matrix for NN: edge weight 1 if NN, 0 otherwise.
  adjMatNN <- matrix(0,n,n)
  for (i in c(1:n)) {
    adjMatNN[indkNN[,i],i] <- 1
  }
  # when an undirected graph is constructed from a non-symm adjMat, the adjMat of the resulting graph iwill be symm.
  nng <- graph.adjacency(adjMatNN, mode="undirected", weighted=TRUE) # vertices are connected if one is the other's nn.
  
  # plot.
  # install.packages("cccd") # packages in /tmp/RtmpigwNEp/downloaded_packages
  # library(cccd)
  #   nng() doesn't work if X has >3 col., i.e. feature >3 dim. 
  # X <- matrix(runif(10),ncol=2)
  #   old.par <- par(mfrow=c(1,3))
  #   plot(nng(X,k=1))  # directed nn graph where edges point to the nn.
  #   plot(nng(X,k=1,mutual=TRUE))  # vertices are connected if they are each other's nn.
  #   plot(nng)
  #   par(old.par)
  return(nng)
}
getMST <- function (X) {
  # Purpose: get minimum spanning tree (MST) of X.
  # Input: X, matrix of double -- data matrix where each row is the feature vector of a data point.
  # Output: mst, igraph -- MST of X.
  
  # igraph tutorial http://igraph.org/r/doc/aaa-igraph-package.html
  #install.packages("igraph")
  library(igraph)  
  
  # create adjacency matrix for X: edges weighted by dist(X).
  adjMat <- as.matrix(dist(X))   
  g <- graph.adjacency(adjMat, mode="undirected", weighted=TRUE) # create igraph using adjMat.
  mst <- mst(g) # create MST of g.
  return(mst)
}
getUG <- function (g1, g2) {
  # Purpose: get Union Graph (UG) whose edges are the union of two given graphs. 
  # Inputs: g1, g2, igraph -- graphs whose edges are to be combined.
  # Output: ug, igraph -- union of g1 and g2.
  
  #ug <- add_edges(g1, t(get.edgelist(g2, names=TRUE))) # add_edge() gives duplicacted edges.
  n <- length(V(g1))
  if (n != length(V(g2))) {
    stop("g1 and g2 should have the same number of vertices.")
  }
  adjMat1 <- get.adjacency(g1); adjMat2 <- get.adjacency(g2); # adjMat contains either 1 or 0.
  adjMatUG <- matrix(ifelse((adjMat1 | adjMat2),1,0),n,n) # union edges.
  ug <- graph.adjacency(adjMatUG, mode="undirected", weighted=TRUE)
  return(ug)
}
getSSLPred <- function (X,y,graphOpt,knng,edgeWeightOpt,sigma,laplacianOpt,objFunOpt,lambda) {
  # Purpose: get graph-based SSL prediction.
  # Inpupts:   
  #         X, matrix of double -- data matrix, each row is a data point.
  #         y, vector of int -- existing labels of data points, length(y)<nrow(X).
  #         graphOpt, int -- 0 for minimum spanning tree (MST), 1 for k-nearest neighbour graph (kNNG) and MST, 2 for complete graph (CG).
  #         k, int -- k-NN graph, needed when graphOpt==1;
  #         edgeWeightOpt, int -- 0 for binary, 1 for RBF.
  #         sigma, double -- needed when edgeWeightOpt==1.
  #         laplacianOpt, int -- 0 for unnormalised, 1 for normalised.
  #         objFunOpt, int -- 0 ensures exact match of existing labels (hard-margin), 1 allows mismatch (soft-margin).
  #         lambda, double -- needed when objFunOpt==1.
  # Output: yPred, vector of int -- predicted labels.
  library(MASS)
  
  n <- nrow(X) # total number of data points.
  nl <- length(y) # number of labelled data point.
  
  # get graph.
  if (graphOpt==0) { # MST
    g <- getMST(X)
  } else if (graphOpt==1) { # MST and kNNG
    mst <- getMST(X)
    nng <- getNNG(X, knng)
    g <- getUG(nng,mst)
  } else if (graphOpt==2) { # CG
    adjMat <- as.matrix(dist(X))   
    g <- graph.adjacency(adjMat, mode="undirected", weighted=TRUE) # create igraph using adjMat.
  }
  
  # get adjacency matrix.
  adjMatUG <- as.matrix(get.adjacency(g))
  if (objFunOpt == 1) { # soft-margin: allows mismatch of existing labels. 
    # soft-margin solu1: dongle.
    adjMatUG <- cbind(adjMatUG,rbind(diag(lambda,nl),matrix(0,n-nl,nl)))
    adjMatUG <- rbind(adjMatUG,cbind(diag(lambda,nl),matrix(0,nl,n-nl),matrix(0,nl,nl)))
    n <- n+nl
  }
  if (edgeWeightOpt == 0) { # "binary"
    W <- adjMatUG
  } else if (edgeWeightOpt == 1) { # "RBF"
    #sigma <- sd(dist(X))
    Wtmp <- exp(-as.matrix(dist(X))^2/(2*sigma^2))
    W <- Wtmp*adjMatUG
    
    # adjust small numbers to avoid numerical errors.
    #tmp <- W; tmp[which(W==0)] <- max(tmp); 
    #myEps <- ifelse(min(tmp)<.Machine$double.eps*1e2,.Machine$double.eps*1e2,min(tmp))
    myEps <- .Machine$double.eps*1e3 # make myEps a bit larger than default eps s.t. solve() (and may be others) can use the default eps.
    W[which(sign(W)!=adjMatUG)] <- myEps # fix numerical error s.t. any none-zero entries in adjMatUG is assigned none-zero weight.
    W[which(W>0&W<myEps)] <- myEps       # fix numerical error s.t. the min weight is no smaller than myEps.
  }
  D <- diag(apply(W,1,sum))
  
  # get Laplacian matrix.
  if (laplacianOpt==0) { # unnormalised Laplacian.
    L <- D-W
  } else if (laplacianOpt==1) { # normalised Laplacian.    
    D_invSqrt <- diag(1/sqrt(diag(D))) # compute inverse for diagonal matrix.
    L <- D-W
    L <- D_invSqrt %*% L %*% D_invSqrt
  }
  
  # get yPred.
  if (objFunOpt == 0) { # exact match of existing labels. 
    #yPred <- -ginv(t(L[(nl+1):n,(nl+1):n]))%*%L[(nl+1):n,1:nl]%*%y # why need t()??
    yPred <- -solve(L[(nl+1):n,(nl+1):n], L[(nl+1):n,1:nl]%*%y, tol=.Machine$double.eps*1e-3) # solve() can't handle singular mareices.
    
    #    xhat1 faster than xhat2 faster than xhat3
    #     system.time({xhat1 = solve(A,b)})
    #     system.time({xhat2 = solve(A) %*% b})
    #     system.time({xhat3 = ginv(A) %*% b})
  } else if (objFunOpt == 1) { # allow mismatch of existing labels.
    # Note: lambda>1e8 may cause error.
    #yPred <- -ginv(t(L[1:(n-nl),1:(n-nl)]))%*%L[1:(n-nl),(n-nl+1):n]%*%y  # why need t()??
    yPred <- -solve(L[1:(n-nl),1:(n-nl)], L[1:(n-nl),(n-nl+1):n]%*%y)
    yPred <- yPred[(nl+1):(n-nl)]
    # soft-margin solu2: iterative solution (deprecated).
    #     yPred0 <- runif(n-nl, -1, 1)
    #     repeat {
    #       y0 <- ginv(lambda*diag(1,nl)+L[1:nl,1:nl]) %*% (lambda*y-L[1:nl,(nl+1):n]%*%yPred0)
    #       yPred <- -ginv(t(L[(nl+1):n,(nl+1):n]))%*%L[(nl+1):n,1:nl]%*%y0
    #       #print(y0);print(yPred0);print(yPred);print(sqrt(sum((yPred0-yPred)^2)))
    #       if (sqrt(sum((yPred0-yPred)^2))<1e-15) {
    #         break
    #       } else {yPred0 <- yPred}
    #     }
  }
  #return(as.factor(sign(yPred)))
  #return(sign(yPred))
  return(yPred)
}
getSSLPred_multiClass<- function (XTrain, yTrain, XTest, laplacianOpt, edgeWeightOpt, knng, sigma) {
  # Purpose: get multi-class prediction using majority vote of one-vs-one scheme.
  # Input: XTrain -- labeled samples from all classes.
  # Output: yPred -- multi-class predictions.
  
  classNames <- levels(factor(yTrain)); nClass <- length(classNames);
  nLabeled <- length(which(yTrain==classNames[1])) # number of labeled samples per class, use class1 as eg.
  nTest <- nrow(XTest)
  # one-vs-one.
  cs <- matrix(1,nTest,(nClass*(nClass-1))/(2)) # number of labels for binary classification.
  ic <- 0
  for (c1 in c(1:(nClass-1))) {
    for (c2 in c((c1+1):nClass)) {
      print(c(classNames[c1],classNames[c2]))
      ic <- ic+1
      ind <- c(which(yTrain==classNames[c1]), which(yTrain==classNames[c2]))
      X <- rbind(XTrain[ind,], XTest, XTrain[-ind,])
      y <- c(rep(1,nLabeled), rep(-1,nLabeled))
      system.time({yTmp <- getSSLPred(X, y, graphOpt=1, knng, edgeWeightOpt, sigma, laplacianOpt, objFunOpt=0, lambda=0)})
      yTmp <- sign(yTmp[1:nTest])
      cs[which(yTmp==1),ic] <- classNames[c1]
      cs[which(yTmp==-1),ic] <- classNames[c2]
    }
  }
  # get multi-class prediction by getting the majority vote from cs.
  yPred <- 0
  for (ix in c(1: nTest)) {
    resTable <- table(cs[ix,])
    yPred[ix] <- as.numeric(names(which.max(resTable)))
  }
  return(yPred)
}