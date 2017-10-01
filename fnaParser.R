# date: 20150827
##### Auxiliary func to get infoAll and nvAll #####
getNumSegments <- function(file_virus) {
  # Purpose: get number of segments for each virus.
  # Input: file_virus, vector of strings -- list of file names
  # Output: numSeg, a col vector -- elem is the number of segments a virus has.
  m <- length(file_virus)
  numSeg <- rep(0, m)
  for(iVirus in c(1:m)) { # get NC files of each virus.
    file_NC <- list.files(file_virus[iVirus]);
    numSeg[iVirus] <- length(file_NC)
  }
  return(numSeg)
}
getName <- function(s) {  
  # Purpose: get name of virus.
  # Input a string s
  # Output a string of virus name.
  terms <- strsplit(s, "  ")
  return(terms[[1]][length(terms[[1]])]) # get name
}
mergeLines <- function(lines, il) {
  # Purpose: search for period and merge all lines before period.
  # Input: lines, string vector -- each elem has one line.
  #        il, int -- start index at vector lines
  # Ouput: a string with one line.
  i <- 1  # i lines in total.
  posPeriod <- gregexpr("\\.",lines[il+i])
  while (posPeriod < 0) {
    i <- i + 1
    posPeriod <- gregexpr("\\.",lines[il+i])
  }
  # concatenate next i lines
  newLine <- ""
  for (j in c(1:i)) {  
    newLine <- paste(newLine, lines[il+j])
  }
  return(newLine)
}
getBaltimoreClass <- function(s, nucleic) {
  # Purpose: get BaltClass.
  # Input: s, string -- the text to search from.
  #        nucleic, string -- nucleic of the virus, either DNA or RNA.
  # Output: baltimoreClass, int.
  baltimoreClass <- -1
  # unknown BaltClass called "unassigned", unknow ICTV class called "unclassified"
  if (gregexpr("Satellites", s)[[1]][1]>0) {
    baltimoreClass <- 8
  } else if (nucleic == "DNA") {  # || gregexpr("DNA", s)>0) {
    if (gregexpr("Retro-transcribing", s)[[1]][1]>0) {
      baltimoreClass <- 7
    } else if (gregexpr("dsDNA", s)[[1]][1]>0) {
      baltimoreClass <- 1
    } else if (gregexpr("ssDNA", s)[[1]][1]>0) {
      baltimoreClass <- 2
    }
  } else if (nucleic == "RNA"){  # || gregexpr("RNA", s)>0) {
    if (gregexpr("Retro-transcribing", s)[[1]][1]>0) {
      baltimoreClass <- 6
    } else if (gregexpr("negative-strand", s)[[1]][1]>0) {
      baltimoreClass <- 5
    } else if (gregexpr("positive-strand", s)[[1]][1]>0) {
      baltimoreClass <- 4
    } else if (gregexpr("dsRNA", s)[[1]][1]>0) {
      baltimoreClass <- 3
    }
  }# else if (gregexpr("unassigned", s)[[1]][1]>0) { 
  #  baltimoreClass <- 9
  #} 
  return(baltimoreClass)
}
getICTVGroups <- function(s) {
  # Purpose: get ICTV group of order, family, subfamily, genera
  # Input: s, string -- the text to seach from.
  # Output: groups, vector of 4 strings -- elem is the taxa.
  groups <- c("","","","")
  if (gregexpr("Satellites", s)>0) {
    return(groups)
  }
  terms <- strsplit(s, ";")
  gps <- c("virale","viridae","virinae","virus")  # order, family, subfamily, genera
  for (ig in c(1:length(gps))) {
    target <- paste("[A-Za-z0-9]",gps[ig],sep="")
    isExist <- gregexpr(target, terms[[1]])>0
    i <- which(isExist==TRUE)[1]
    if (!is.na(i) && gregexpr("unclassified", terms[[1]][i])<0 && gregexpr("unassigned", terms[[1]][i])<0) {
      tms <- strsplit(terms[[1]][i], " ")
      ii <- which(gregexpr(target, tms[[1]])>0) 
      groups[ig] <- gsub("[ \\.\n\t]", "", tms[[1]][ii])
    }
    # don't do anything since info not available
  }
  return(groups)
}
getNV <- function (s, alphabet, myStat, myKmer) {
  # Purpose: get k-mer NV (natural vector). 
  # Inputs: s, string -- the nucleotide sequence.
  #         alphabet, array -- alphabets, e.g. c('a', 'c', 'g', 't') etc.
  #         myStat, int -- up to (myStat-1)-th order moment, i.e. NV_{4*myStat}.
  #         myKmer, int -- k-mer.
  
  # initialise parameters.
  s=strsplit(s,"")[[1]] # convert string to char array.
  s <- s[which(s %in% alphabet)] # remove none ACGT.
  sLength <- length(s)
  kmerListLength <- sum(length(alphabet)^seq(1,myKmer))
  
  # generate all possible combinations of kmer for a given k.
  kmerList <- alphabet; kmerList0 <- alphabet;
  if (myKmer > 1) {
    for (j in 2:myKmer) {
      kmerList0 <- c(t(outer(kmerList0, alphabet, FUN="paste", sep=""))) # for kmer of a particular k.
      kmerList <- c(kmerList, kmerList0) # for kmer of all k = 1..k.
    }
  }
  
  # create a map whose keys are kmers and values are their positions.
  e <- new.env()
  # create keys.
  for (i in c(1:kmerListLength)) {
    assign(sprintf("e$%s",kmerList[i]),0)
    # get(sprintf("e$%s",kmerList[i]))
  }
  # get positions of each kmer.
  for (j in c(1:myKmer)) { # for eack k of k-mer
    for (i in c(1:(length(s)-j+1))) { # foreach char in the sequence s.
      vn <- sprintf("e$%s", paste(s[i:(i+j-1)],sep="",collapse=""))
      if (get(vn)[1] == 0) { assign(vn, i)
      } else { assign(vn, c(get(vn),i))}
    }
  }
  
  # compute count, mean and variance for each kmer using the env() e, i.e. get NV_{4*(j+1)} for each kmer.
  nv <- matrix(0,kmerListLength,myStat) # each matrix row is a kmer, col is a stat.
  for (i in c(1:kmerListLength)) {
    pos <- get(sprintf("e$%s",kmerList[i]))
    if (pos[1]!=0) {
      nv[i,1] <- length(pos)
      nv[i,2] <- mean(pos)
      for (j in c(2:(myStat-1))) {
        nv[i,j+1] <- sum((((pos-nv[i,2])/(nv[i,1]*(sLength-myKmer+1)))^(j-1))*(pos-nv[i,2]))
      }
    } else {
      nv[i,1] <- 0
      nv[i,2] <- 0
      for (j in c(2:(myStat-1))) {
        nv[i,j+1] <- 0
      }
    }
  }
  return(c(nv)) # make vector.
}
getRTD <- function (s, alphabet, myKmer) {
  # Purpose: get k-mer RTD (return time distribution). 
  # Inputs: s, string -- the nucleotide sequence.
  #         alphabet, array -- alphabets, e.g. c('a', 'c', 'g', 't') etc.
  #         myKmer, int -- k-mer.
  # Output: x, matrix of double -- mean and sd of return time for each kmer.
  
  # initialise parameters.
  s=strsplit(s,"")[[1]] # convert string to char array.
  s <- s[which(s %in% alphabet)] # remove none ACGT.
  sLength <- length(s)
  kmerListLength <- sum(length(alphabet)^seq(1,myKmer))
  
  # generate all possible combinations of kmer for a given k.
  kmerList <- alphabet; kmerList0 <- alphabet;
  if (myKmer > 1) {
    for (j in 2:myKmer) {
      kmerList0 <- c(t(outer(kmerList0, alphabet, FUN="paste", sep=""))) # for kmer of a particular k.
      kmerList <- c(kmerList, kmerList0) # for kmer of all k = 1..k.
    }
  }
  
  # create a map whose keys are kmers and values are their positions.
  e <- new.env()
  # create keys.
  for (i in c(1:kmerListLength)) {
    assign(sprintf("e$%s",kmerList[i]),0)
    # get(sprintf("e$%s",kmerList[i]))
  }
  # get positions of each kmer.
  for (j in c(1:myKmer)) { # for eack k of k-mer
    for (i in c(1:(length(s)-j+1))) { # foreach char in the sequence s.
      vn <- sprintf("e$%s", paste(s[i:(i+j-1)],sep="",collapse=""))
      if (get(vn)[1] == 0) { assign(vn, i)
      } else { assign(vn, c(get(vn),i))}
    }
  }
  
  # compute mean and variance of RTD for each kmer using the env() e.
  X <- matrix(0,kmerListLength,2) # each matrix row is a kmer, col is a stat.
  for (i in c(1:kmerListLength)) {
    pos <- get(sprintf("e$%s",kmerList[i]))
    if (length(pos)>=2) {
      rtd <- pos[2:length(pos)]-pos[1:(length(pos)-1)]-1 # return time distribution.
      X[i,1] <- mean(rtd)
      X[i,2] <- ifelse(length(rtd)==1,0,sd(rtd))
    } else {
      X[i,1] <- 0
      X[i,2] <- 0
    }
  }
  return(c(X)) # make vector.
}
getFeature <- function (s, feature,alphabet, myStat, myKmer) {
  # Purpose: get feature vector.
  # Input: feature, string -- options: "nv", "rtd".
  #        myStat, int -- needed for feature="nv".
  
  if (feature=="nv"){
    X <- getNV(s, alphabet, myStat, myKmer)
  } else if (feature=="rtd"){
    X <- getRTD(s, alphabet, myKmer)
  }
  return (X)
}
##### Get infoAll and nvAll #####
getVirus <- function(fname, needGV, feature,alphabetOpt, myKmer, myStat) {
  # Purpose: get info of a virus from .fna file.
  # Input: fname, string -- name of a .gbk file containing info of a virus.
  #        needGV, bool -- whether need to use Vectorizer() for GV.
  #        mykmer, int -- needed for k-mer NV.
  #        myStat, int -- up to (myStat-1)-th moment: 1 count, 2 mean, 3 var, 4 skew, 5 kurtosis.
  #        alphabetOpt, string -- options: "acgt", "sw", "ry", "mk".
  # Output: list(info, nv), list -- a list of 2 data frames, one for info, another NV.
  
  # read file
  fgbk <- file(fname, "r");
  lines <- readLines(fgbk);
  close(fgbk);
  # parse lines 
  accNum<-""; gi<-""; name<-""; topology<-"";baltimoreClass<-""; groups<-c("","","","");
  lineage<-"";
  for (il in c(1:length(lines))) {
    # get topology: linear or circular
    if (il==1 && gregexpr("LOCUS",lines[il],ignore.case=FALSE) > 0) { 
      if (gregexpr("linear",lines[il],ignore.case=TRUE) > 0) { 
        topology <- "linear"
      } else if (gregexpr("circular",lines[il],ignore.case=TRUE) > 0) {
        topology <- "circular"
      } else {
        stop("neither linear nor circular")
      }
      if (gregexpr("DNA",lines[il],ignore.case=FALSE) > 0) {
        nucleic <- "DNA"
      } else if (gregexpr("RNA",lines[il],ignore.case=FALSE) > 0) {
        nucleic <- "RNA"
      } else {
        stop("neither DNA nor RNA")
      }
    }
    # get accession number
    else if (gregexpr("VERSION",lines[il],ignore.case=FALSE) > 0) { 
      terms <- strsplit(lines[il], "  ")
      for (it in c(2:length(terms[[1]]))) { # it=1 is "VERSION"
        if (terms[[1]][it]!="" && accNum=="") {
          accNum <- gsub(" ", "", terms[[1]][it])
        }
        if (gregexpr("GI:",terms[[1]][it],ignore.case=FALSE) > 0) {
          gi <- strsplit(lines[il], ":")[[1]][2]
        }
        #         if (gregexpr("NC_",terms[[1]][it],ignore.case=FALSE) > 0) {
        #           accNum <- gsub(" ", "", terms[[1]][it])
        #           break
        #         }
      }
    }
    # get taxa
    else if (gregexpr("ORGANISM",lines[il],ignore.case=FALSE) > 0) { 
      name <- getName(lines[il])  # virus name
      newLine <- mergeLines(lines, il)
      baltimoreClass <- getBaltimoreClass(newLine, nucleic)  # Baltimore class
      groups <- getICTVGroups(newLine)  # ICTV groups
      lineage <- gsub("  ", "",newLine)
    }
    # get nucleotide sequence
    else if (gregexpr("ORIGIN",lines[il],ignore.case=FALSE) > 0) { 
      genome <- ""
      for (i in c((il+1):(length(lines)-1))) {
        genome <- paste(genome, gsub("[0-9 ]","",lines[i]), sep="")
      }
      break
    }
  }
  info <- data.frame(gi, accNum, name, topology, baltimoreClass, order=groups[1],
                     family=groups[2],subfamily=groups[3],genera=groups[4],lineage)
  #   if (needGV) {
  #     gv <- Vectorizer(genome, kmer=myKmer, statistic=myStat)
  #     res <- list(info, gv) 
  #   } else {
  if (alphabetOpt=="acgt") {
    X <- getFeature(genome,feature,c('a', 'c', 'g', 't'), myStat, myKmer); # same as gv1mer <- Vectorizer(genome, kmer=1, statistic=3)
  } else if (alphabetOpt=="sw") { # there are some swrymk in the original sequence.
    genome <- gsub("g", "s",genome); genome <- gsub("c", "s",genome);
    genome <- gsub("a", "w",genome); genome <- gsub("t", "w",genome);
    X <- getFeature(genome,feature,c('s','w'), myStat, myKmer);
  } else if (alphabetOpt=="ry") {
    genome <- gsub("a", "r",genome); genome <- gsub("g", "r",genome);
    genome <- gsub("c", "y",genome); genome <- gsub("t", "y",genome);
    X <- getFeature(genome,feature,c('r','y'), myStat, myKmer);
  } else if (alphabetOpt=="mk") {
    genome <- gsub("a", "m",genome); genome <- gsub("c", "m",genome);
    genome <- gsub("g", "k",genome); genome <- gsub("t", "k",genome);
    X <- getFeature(genome,feature,c('m','k'), myStat, myKmer);
  }  
  res <- list(info, X)
  #  }
  return(res)
}
getVirusAll <- function(feature,alphabetOpt,myStat,myKmer) {
  # Purpose: use getVirus() to get info, nv of all viruses.
  # Input: myStat, int -- up to (myStat-1)-th normalised central moment of nucleotide.
  #        input, String -- full path to input folder all.gbk.
  #        output, String -- full path to output file for list(nvAll, infoAll, numSeg).
  # Output: list(nvAll, infoAll, numSeg)
  #         nvAll, matric of double -- each row is nv of a virus.
  #         infoAll, matrix of data frame -- each row is info of a virus.
  #         numSeg, vector of int -- number of segments of viruses.
  #input <- "/cs/research/intelsys/home1/tingtwan/Documents/VirusGenome/Data/20150819all.gbk" # local computer folder.
  #input <- "/cs/research/vision/projects0/tingtwan/VirusGenome/Data/20150918all.gbk" # local project folder.
  #output <- "/cs/research/vision/projects0/tingtwan/phd"
  input <- "/cluster/project2/ADHD/VirusGenome/ClusterInput/20150819all.gbk" # cluster project folder.
  output <- "/cluster/project2/ADHD/VirusGenome/ClusterInput/single_sw_info_nv_kmer=1to6_stat=1to3.RData" # cluster project folder.
  
  segOpt <- "single"  # options: "single","multiple","all".
  needGV <- FALSE # whether use code from GV author.
  if (alphabetOpt=="acgt") {
    xLength <- 4 # length of alphabet. 
  } else if (alphabetOpt=="sw" || alphabetOpt=="ry" || alphabetOpt=="mk") {
    xLength <- 2  
  }
  if (feature=="rtd") {
    myStat <- 2
  }
  nvDim <- sum(xLength^seq(1:myKmer))*myStat # dimension of nv.
  
  fileVirus <- list.files(input,full.names=TRUE)
  m <- length(fileVirus) # total number of viruses (single and multiple segmented)
  numSeg <- getNumSegments(fileVirus) # a m*1 vector
  ind1 <- which(numSeg == 1) # index of single-segmented
  ind2 <- which(numSeg > 1) # index of mult-segmented
  
  if (segOpt == "single") { # single-segmented only
    #myKmer <- 2; myStat <- 3;
    #     if (needGV) {
    #       nvSingle <- matrix(0,m1,ifelse(myKmer==2, 61, 253)); # each row is a virus.
    #     } else {
    #       nvSingle <- matrix(0,m1,(myStat*4)); # each row is a virus.
    #     }
    nvSingle <- matrix(0,length(ind1),nvDim) # each row is a virus.
    for (i in c(1:length(ind1))) { # get NC files of each virus.
      fileNC <- list.files(fileVirus[ind1[i]],full.names=TRUE)
      print(sprintf("i=%g: %s",i,fileNC))
      virus <- getVirus(fileNC, needGV, feature,alphabetOpt, myKmer, myStat)
      nvSingle[i, ] <- virus[[2]]
      if (exists("infoSingle")) {
        infoSingle <- rbind(infoSingle, virus[[1]])
      } else {
        infoSingle <- virus[[1]]
      }
    }
  } else if (segOpt == "multiple") {
    nvMultiple <- matrix(0,sum(numSeg[ind2]),nvDim)
    iSeg <- 1
    for (iVirus in c(1: length(ind2))) {
      fileNC <- list.files(fileVirus[ind2[iVirus]],full.names=TRUE)
      for (iNC in c(1:length(fileNC))) {
        print(c(iVirus,iNC))
        virus <- getVirus(fileNC[iNC], needGV, feature,alphabetOpt,myKmer, myStat)
        nvMultiple[iSeg,] <- virus[[2]]; iSeg <- iSeg+1;
        if (exists("infoMultiple")) {
          infoMultiple <- rbind(infoMultiple, virus[[1]])
        } else {
          infoMultiple <- virus[[1]]
        }
      }
    }
  } else if (segOpt == "all") {
    totalNumSeg <- sum(numSeg)
    nvAll <- matrix(0,totalNumSeg,nvDim);
    i <- 0
    for (iVirus in c(1:m)) {
      fileNC <- list.files(fileVirus[iVirus],full.names=TRUE);
      for (iNC in c(1:length(fileNC))) {
        i <- i+1
        virus <- getVirus(fileNC[iNC], needGV, feature,alphabetOpt,myKmer, myStat)
        nvAll[i, ] <- virus[[2]]
        if (exists("infoAll")) { # info repeats x times when a virus has x segments.
          print(i)
          infoAll <- rbind(infoAll, virus[[1]])
        } else {
          infoAll <- virus[[1]]
        }
      }
    }
  }
  #  write.csv(infoAll$baltimoreClass, file = "~/Documents/VirusGenome/Data/baltAll.csv",row.names=FALSE)
  #  write.csv(nvAll, file = "~/Documents/VirusGenome/Data/nvAll.csv",row.names=FALSE)
  result <- list(nvAll, infoAll, numSeg)
  #  save(result, file=output)
  return(result)
}
##### Get features and labels #####
getXY <- function (fileName, taxon, feature, myKmer1, myKmer2, myStat1, myStat2,xLength,dataKmer,dataStat) {
  # Purpose: get data matrix X and labels y from single_info_nv12mer5.Rdata.
  # Input:taxon, string -- options: balt, unbalt, order, unorder, host3, host4, host.
  #       feature, string -- options: nv, genome length, nucleutide counts, combined nucleutide counts, extended nvr
  #       myKmer12, myStat12, int -- needed when feature=="NV". Start from myKmer1 and end by myKmer2.
  #       fileName, string -- full path to the .RData file containing nv and info.
  #       xLength, int --  4, number of alphabets a,c,g,t.
  #       dataKmer, int -- 6, max. kmer of nv in <fileName> from which XY is extracted.
  #       dataStat, int -- 3, max stat of nv in <fileName> from which XY is extracted.
  
  # Output: list(X, y): X, matrix of double -- the designed data matrix.
  #                         y, vector of factors -- the labels.
  
  # get y.
  if (substring(taxon,1,3)=="all") { # all seq. in the dataset.
    #fileName="/cs/research/vision/projects0/tingtwan/VirusGenome/ClusterInput/single_info_nv_kmer=1to5_stat=1to3.RData"
    load(fileName)
    ind <- c(1:nrow(nv))
    y <- info[ind,]
  } else if (substring(taxon,1,5)=="virus") { # non-satellites.
    #fileName="/cs/research/vision/projects0/tingtwan/VirusGenome/ClusterInput/single_info_nv_kmer=1to5_stat=1to3.RData"
    load(fileName)
    ind <- which(info$baltimoreClass != 8)
    y <- info[ind,]
  } else if (substring(taxon,1,4)=="balt") { # Baltimore Group labelled.
    #fileName="/cs/research/intelsys/home1/tingtwan/Documents/VirusGenome/ClusterInput/nv40_info_single_corrected.RData"
    load(fileName)
    ind <- which(info$baltimoreClass %in% c(1:7))  # only these from BaltClass 1-7
    y <- as.factor(info$baltimoreClass[ind])
    y <- factor(y) # update the number of factor levels
  } else if (substring(taxon,1,6)=="unbalt") { # Baltimore Group unlabelled.
    #fileName="/cs/research/intelsys/home1/tingtwan/Documents/VirusGenome/ClusterInput/nv40_info_single_corrected.RData"
    load(fileName)
    ind <- which(info$baltimoreClass == -1)  
    y <- info[ind,]
  } else if (substring(taxon,1,5)=="order") { # ICTV Order labelled.
    #fileName="/cs/research/intelsys/home1/tingtwan/Documents/VirusGenome/ClusterInput/nv40_info_single_corrected.RData"
    load(fileName)
    ind <- which(info$order != "" & info$baltimoreClass!=8)
    y <- as.factor(info$order[ind])
    y <- factor(y)
  } else if (substring(taxon,1,7)=="unorder") { # ICTV Order unlabelled.
    #fileName="/cs/research/intelsys/home1/tingtwan/Documents/VirusGenome/ClusterInput/nv40_info_single_corrected.RData"
    load(fileName)
    ind <- which(info$order == "" & info$baltimoreClass != 8) # remove Satellites.
    y <- info[ind,]
  } else if (substring(taxon,1,5)=="3host"){ # 3 host labels: Algea, Bactoria, Euc.
    #fileName="/cs/research/vision/projects0/tingtwan/VirusGenome/ClusterInput/single_info_nv_kmer=1to5_stat=1to3.RData"
    load(fileName) # NV12, host and other info.
    ind <- which(info$host!="" & info$host!="environment" & info$baltimoreClass!=8) # single-segmented and host labelled viruses.
    y0 <- info$host[ind]
    y <- rep("E",length(y0))
    y[which(y0=="archaea")] <- "A"; y[which(y0=="bacteria")] <- "B";
    y <- as.factor(y)
    y <- factor(y)
  } else if (substring(taxon,1,5)=="4host"){ # 4 host labels: Algea, Bactoria, EucA, EucB.
    #fileName="/cs/research/vision/projects0/tingtwan/VirusGenome/ClusterInput/single_info_nv_kmer=1to5_stat=1to3.RData"
    load(fileName) # NV12, host and other info.
    ind <- which(info$host!="" & info$host!="environment" & info$baltimoreClass!=8) # single-segmented and host labelled viruses.
    y0 <- info$host[ind]
    y <- rep("E2",length(y0))
    y[which(y0=="archaea")] <- "A"; y[which(y0=="bacteria")] <- "B";
    y[which(y0=="algae" | y0=="fungi" | y0=="protozoa")] <- "E1";
    y <- as.factor(y)
    y <- factor(y)
  } else if (substring(taxon,1,4)=="host"){ # original host labels from NCBI file (should check this after checking substring(taxon,1,5)).
    #fileName="/cs/research/vision/projects0/tingtwan/VirusGenome/ClusterInput/single_info_nv_kmer=1to5_stat=1to3.RData"
    load(fileName) # NV12, host and other info.
    # single-segmented and host labelled viruses. Excluding "invertebrates, vertebrates" b/c only 4 -- too few.
    ind <- which(info$host!="" & info$host!="invertebrates, vertebrates" & info$baltimoreClass!=8) 
    y <- as.factor(info$host[ind])
    y <- factor(y)
    #   else if (substring(taxon,1,5)=="3host"){ # 3 host labels: Algea, Bactoria, Euc.
    #     #fileName="/cs/research/vision/projects0/aoe-collaboration/20151211all.RData"
    #     load(fileName) # NV12, host and other info.
    #     nv <- virus[,1:12]; info <- virus[,13:23];
    #     ind <- which(info$numSeg==1 & info$host!="" & info$host!="environment" & info$baltimoreClass!=8) # single-segmented and host labelled sviruses.
    #     y0 <- info$host[ind]
    #     y <- rep("E",length(y0))
    #     y[which(y0=="archaea")] <- "A"; y[which(y0=="bacteria")] <- "B";
    #     y <- as.factor(y)
    #     y <- factor(y)
    #   } else if (substring(taxon,1,5)=="4host"){ # 4 host labels: Algea, Bactoria, EucA, EucB.
    #     #fileName="/cs/research/vision/projects0/aoe-collaboration/20151211all.RData"
    #     load(fileName) # NV12, host and other info.
    #     nv <- virus[,1:12]; info <- virus[,13:23];
    #     ind <- which(info$numSeg==1 & info$host!="" & info$host!="environment" & info$baltimoreClass!=8) # single-segmented and host labelled sviruses.
    #     y0 <- info$host[ind]
    #     y <- rep("E2",length(y0))
    #     y[which(y0=="archaea")] <- "A"; y[which(y0=="bacteria")] <- "B";
    #     y[which(y0=="algae" | y0=="fungi" | y0=="protozoa")] <- "E1";
    #     y <- as.factor(y)
    #     y <- factor(y)
    #   } else if (substring(taxon,1,4)=="host"){ # original host labels from NCBI file (should check this after checking substring(taxon,1,5)).
    #     #fileName="/cs/research/vision/projects0/aoe-collaboration/20151211all.RData"
    #     load(fileName) # NV12, host and other info.
    #     nv <- virus[,1:12]; info <- virus[,13:23];
    #     ind <- which(info$numSeg==1 & info$host!="" & info$host!="invertebrates, vertebrates" & info$baltimoreClass!=8) # single-segmented and host labelled sviruses.
    #     y <- as.factor(info$host[ind])
    #     y <- factor(y)
  } else {
    stop("unknown taxon")
  }
  
  # get X.
  if (feature=="length") {  # genLength
    X <- apply(nv[ind, 1:xLength], 1, sum)
    # X <- X/max(X)
  } else if (feature=="nv") {  # k-mer nv
    dataKmerListLength <- sum(xLength^seq(1,dataKmer))
    
    # extract the part we need.
    indMyKmer1 <- ifelse(myKmer1==1,1,1+sum(xLength^seq(1,(myKmer1-1))))
    indMyKmer2 <- sum(xLength^seq(1,myKmer2))
    nv1 <- nv[ind,]
    X <- matrix(0,nrow(nv1),(indMyKmer2-indMyKmer1+1)*(myStat2-myStat1+1))
    for(i in c(1:nrow(nv1))) {
      # rearrange vector nv to matrix nv.
      nvTmp <- matrix(nv1[i,],dataKmerListLength,dataStat)
      X[i,] <- c(nvTmp[indMyKmer1:indMyKmer2,myStat1:myStat2])
    }
  } else if (feature=="means") {
    X <- nv[ind, c(5:8)]
  } else if (feature=="gc") {
    na <- nv[ind,1]; nc <- nv[ind,2]; ng <- nv[ind,3]; nt <- nv[ind,4]; 
    X <- (ng+nc)/(na+nc+ng+nt)
  } else if (feature=="combCounts") {
    #na <- nv[ind,1]; nc <- nv[ind,4]; ng <- nv[ind,7]; nt <- nv[ind,10]; # nvSingle_*
    na <- nv[ind,1]; nc <- nv[ind,2]; ng <- nv[ind,3]; nt <- nv[ind,4]; # nv40_info_*
    X <- matrix(c(na+ng, nc+nt, na+nc, ng+nt, na+nt, ng+nc), length(ind), 6)
  } else if (feature=="nvExt1") {  # Mark's
    X1 <- nv[ind,]
    k <- 0
    X2 <- matrix(0,length(ind),(11*12/2))
    for (i in c(1:11)) {
      for (j in c((i+1):12)) {
        k <- k+1
        X2[,k] <- nv[ind,i]*nv[ind,j]
      }
    }
    X <- cbind(X1,X2)
  } else if (feature=="nvExt2") {  # matches NV svmPoly gamma=1,coef=0,degree=2 exactly.
    X1 <- nv[ind,]*nv[ind,]
    k <- 0
    X2 <- matrix(0,length(ind),(11*12/2))
    for (i in c(1:11)) {
      for (j in c((i+1):12)) {
        k <- k+1
        X2[,k] <- nv[ind,i]*nv[ind,j]
      }
    }
    X <- cbind(X1,X2)
  } else {
    stop("unknown feature")
  }
  return(list(as.matrix(X),y))
}