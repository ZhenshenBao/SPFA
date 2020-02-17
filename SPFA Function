SPFA <- function(paths=NULL,mdir=NULL,de=NULL,all=NULL,normal=NULL,cancer=NULL){
  nB <- 2000
  pb <- NULL
  pl <- NULL
  leafscores <- NULL
  Names <- NULL
  pathwayIDs <- NULL
  for(i in 1:length(paths)){
    mapkpathway <- try(parseKGML(paste(mdir,paths[i],sep = "/")),TRUE)
    nodes <- nodes(mapkpathway)
    mapkG <- KEGGpathway2Graph(mapkpathway,genesOnly=FALSE,expandGenes=FALSE)
    mapkG1 <- KEGGpathway2Graph(mapkpathway,expandGenes=TRUE)
    mapkG2 <- KEGGpathway2Graph(mapkpathway,genesOnly=FALSE,expandGenes=TRUE)
    allgenes <- nodes(mapkG1)
    allnodes <- nodes(mapkG2)
    allgeneslist <- strsplit(allgenes,"hsa:")
    allgeneslist <- do.call(rbind,allgeneslist)
    allgenes <- allgeneslist[,2]
    allnodeslist <- strsplit(allnodes,"hsa:")
    allnodeslist <- do.call(rbind,allnodeslist)
    allnodes <- allnodeslist[,2]
    g <- igraph.from.graphNEL(mapkG)
    outdegree <- igraph::degree(g,mode ="out")
    Ue <- names(outdegree[which(outdegree==0)])
    Uenodes <- nodes[Ue]
    types <- sapply(Uenodes, getType)
    endnodes <- names(types[which(types=="gene")])
    nonodes <- setdiff(Ue,endnodes)
    indegree<-igraph::degree(g,mode ="in")
    Us <- names(indegree[which(indegree==0)])
    nonodes <- setdiff(nonodes,Us)
    
    dismatrix <- distances(g,mode="in")
    dismatrix[which(dismatrix=="Inf")] <- 0
    nn <- max(dismatrix)
    dismatrix[which(dismatrix!=0)] <- 1
    
    ni <- 1
    while(length(nonodes)!=0 & ni<=nn){
      adjnodes<-NULL
      for(j in 1:length(nonodes)){
        adjnode <- neighbors(g,nonodes[j],"in")
        for(ij in 1:length(adjnode)){
          adjtype <- nodes[[adjnode[ij]]]
          nonodes <- setdiff(nonodes,names(adjnode[ij]))
          if(adjtype@type == "gene"){
            endnodes <- c(endnodes,names(adjnode[ij]))
          }
          if(adjtype@type != "gene"){
            adjnodes <- c(adjnodes,names(adjnode[ij]))
          }
        }
      }
      adjnodes <- setdiff(adjnodes,Us)
      nonodes <- adjnodes
      ni <- ni+1
    }
    endnodes <- endnodes[!duplicated(endnodes)]

    endtypes <- nodes[endnodes]
    types1 <- sapply(endtypes, getType)
    endnodes <- names(types1[which(types1=="gene")])
    endgenes <- NULL
    for(ii in 1:length(endtypes)){
      if(length(endtypes[[ii]]) != 0){
        endgene <- endtypes[[ii]]@name
        endgenelist <- strsplit(endgene,"hsa:")
        endgenelist <- do.call(rbind,endgenelist)
        endgene <- endgenelist[,2]
        endgenes <- c(endgenes,endgene)
      }
    }
    endgenes <- endgenes[!duplicated(endgenes)]
  
    nde <- intersect(de,allgenes)
    pl[i] <- phyper(q=length(nde)-1,m=length(allgenes),n=length(all)-length(allgenes),k=length(de),lower.tail=FALSE)
    conames <- allgenes
    difnames <- setdiff(conames,colnames(normal))
    conames <- intersect(conames,colnames(normal))
    
    
    nfor <- normal[,conames]
    nfor <- as.data.frame(nfor)
    cfor <- cancer[,conames]
    cfor <- as.data.frame(cfor)
    Inormal <- cor(nfor)
    Icancer <- cor(cfor)
    Inormal[is.na(Inormal)] <- 1
    Icancer[is.na(Icancer)] <- 1
    Iscore <- Icancer-Inormal
    if(length(difnames)!=0){
      ainter <- rep(0,ncol(Iscore))
      ai <- matrix(rep(ainter,length(difnames)),ncol=length(difnames))
      binter <- rep(0,length(allgenes))
      bi <- matrix(rep(binter,length(difnames)),nrow=length(difnames))
      Iscore <- cbind(Iscore,ai)
      Iscore <- rbind(Iscore,bi)
      genenames <- c(conames,difnames)
      colnames(Iscore) <- genenames
      rownames(Iscore) <- genenames
      Iscore <- Iscore[allgenes,allgenes]
    }
    Iscore <- abs(Iscore)
    notends <- setdiff(allnodes,endgenes)
    g1 <- igraph.from.graphNEL(mapkG2)
    dismatrix1 <- distances(g1,mode="in")
    dismatrix1[which(dismatrix1=="Inf")] <- 0
    dism <- dismatrix1
    colnames(dism) <- allnodes
    rownames(dism) <- allnodes
    dism <- dism[allgenes,allgenes]
    dismatrix1[which(dismatrix1!=0)] <- 1
    colnames(dismatrix1) <- allnodes
    rownames(dismatrix1) <- allnodes
    dismatrix1[notends,] <- 0
    allgenes <- as.character(allgenes)
    dismatrix1 <- dismatrix1[allgenes,allgenes] 
    
    edgesL <- Iscore*dismatrix1/dism
    edgesL[is.na(edgesL)] <- 0
    edgescore <- sum(edgesL)
    
    zgenes <- NULL
    zgenes1 <- NULL
    for(iiii in 1:ncol(Iscore)){
      if(sum(Iscore[iiii,]) == 0){
        zgenes[iiii] <- rownames(Iscore)[iiii]
      }
      if(sum(Iscore[,iiii]) == 0){
        zgenes1[iiii] <-colnames(Iscore)[iiii]
      }
    }
    zgenes <- zgenes[!is.na(zgenes)]
    zgenes1 <- zgenes1[!is.na(zgenes1)]
    
    score <- NULL
    for(js in 1:nB){
      esams <- sample(1:ncol(normal),length(allgenes))
      nforsam <- normal[,esams]
      nforsam <- as.data.frame(nforsam)
      cforsam <- cancer[,esams]
      cforsam <- as.data.frame(cforsam)
      Inormalsam <- cor(nforsam)
      Icancersam <- cor(cforsam)
      Inormalsam[is.na(Inormalsam)] <- 1
      Icancersam[is.na(Icancersam)] <- 1
      Iscoresam <- Icancersam-Inormalsam
      colnames(Iscoresam) <- allgenes
      rownames(Iscoresam) <- allgenes
      Iscoresam[zgenes,] <- 0
      Iscoresam[,zgenes1] <- 0
      #Iscoresam <- abs(Iscoresam)*dismatrix1*LLm/dism
      Iscoresam <- abs(Iscoresam)*dismatrix1/dism
      Iscoresam[is.na(Iscoresam)] <- 0
      edgesscore <- sum(Iscoresam)
      score[js] <- edgesscore
    }
    
    
    if(edgescore > 0){
      pb[i] <- sum(score >= edgescore)/length(score)
      if(pb[i] <= 0){pb[i]<-1/nB/100} 
      if(pb[i] > 1){pb[i]<-1} 
    }
    if(edgescore == 0){
      if(all(score == 0)){    #there is nothing to learn from perturbations
        pb[i] <- NA
      }else{
        pb[i] <- 1
      }
    }
    if(edgescore < 0){
      pb[i] <- NA
    }
    Names[i] <- mapkpathway@pathwayInfo@title
    pathwayIDs[i] <- mapkpathway@pathwayInfo@number
    cat("\n");
    cat(paste("Done the ",i," th pathway: ",mapkpathway@pathwayInfo@title,sep=""))
  }
  c <- pb*pl
  pc <- c-c*log(c)
  pcFDR = p.adjust(pc,"fdr")
  pcfwer = p.adjust(pc,"bonferroni") 
  res <- data.frame(Names,ID=pathwayIDs,pl=pl,pb=pb,p=pc,pFdr=pcFDR,pFWER=pcfwer,stringsAsFactors=FALSE)
  
  res <- res[!is.na(res$p),]
  res <- res[order(res$p),]
  rownames(res) <- NULL
  return(res)
}


SPFA_per <- function(mdir = NULL, pathfilename = NULL, normal = NULL, cancer = NULL ){
  mapkpathway<-try(parseKGML(paste(mdir,pathfilename,sep = "/")),TRUE)
  nodes<-nodes(mapkpathway)
  mapkG <- KEGGpathway2Graph(mapkpathway,genesOnly=FALSE,expandGenes=FALSE)
  mapkG1 <- KEGGpathway2Graph(mapkpathway,expandGenes=TRUE)
  mapkG2 <- KEGGpathway2Graph(mapkpathway,genesOnly=FALSE,expandGenes=TRUE)
  allgenes <- nodes(mapkG1)
  allnodes <- nodes(mapkG2)
  allgeneslist <- strsplit(allgenes,"hsa:")
  allgeneslist <- do.call(rbind,allgeneslist)
  allgenes <- allgeneslist[,2]
  allnodeslist <- strsplit(allnodes,"hsa:")
  allnodeslist <- do.call(rbind,allnodeslist)
  allnodes <- allnodeslist[,2]
  g<-igraph.from.graphNEL(mapkG)
  outdegree<-igraph::degree(g,mode ="out")
  Ue<-names(outdegree[which(outdegree==0)])
  Uenodes <- nodes[Ue]
  types <- sapply(Uenodes, getType)
  endnodes <- names(types[which(types=="gene")])
  nonodes <- setdiff(Ue,endnodes)
  indegree<-igraph::degree(g,mode ="in")
  Us<-names(indegree[which(indegree==0)])
  nonodes <- setdiff(nonodes,Us)
  
  g2 <- igraph.from.graphNEL(mapkG1)
  outdegree1 <- igraph::degree(g2,mode ="out")
  Ue1 <- names(outdegree1[which(outdegree1==0)])
  indegree1 <- igraph::degree(g2,mode ="in")
  Us1 <- names(indegree1[which(indegree1==0)])
  onenodes <- intersect(Ue1,Us1)
  onenodeslist <- strsplit(onenodes,"hsa:")
  onenodeslist <- do.call(rbind,onenodeslist)
  onenodes <- onenodeslist[,2]
  
  dismatrix<-distances(g,mode="in")
  dismatrix[which(dismatrix=="Inf")]<-0
  nn <- max(dismatrix)
  dismatrix[which(dismatrix!=0)] <- 1
  
  ni <- 1
  while(length(nonodes)!=0 & ni<=nn){
    adjnodes<-NULL
    for(j in 1:length(nonodes)){
      adjnode <- neighbors(g,nonodes[j],"in")
      for(ij in 1:length(adjnode)){
        adjtype <- nodes[[adjnode[ij]]]
        nonodes <- setdiff(nonodes,names(adjnode[ij]))
        if(adjtype@type=="gene"){
          endnodes <- c(endnodes,names(adjnode[ij]))
        }
        if(adjtype@type!="gene"){
          adjnodes <- c(adjnodes,names(adjnode[ij]))
        }
      }
    }
    adjnodes <- setdiff(adjnodes,Us)
    nonodes <- adjnodes
    ni <- ni+1
  }
  endnodes <- endnodes[!duplicated(endnodes)]
  conames <- allgenes
  difnames <- setdiff(conames,colnames(normal))
  conames <- intersect(conames,colnames(normal))
  
  
  endtypes <- nodes[endnodes]
  types1 <- sapply(endtypes, getType)
  endnodes <- names(types1[which(types1=="gene")])
  endgenes <- NULL
  for(ii in 1:length(endtypes)){
    if(length(endtypes[[ii]])!=0){
      endgene <- endtypes[[ii]]@name
      endgenelist <- strsplit(endgene,"hsa:")
      endgenelist <- do.call(rbind,endgenelist)
      endgene <- endgenelist[,2]
      endgenes <- c(endgenes,endgene)
    }
  }
  
  endgenes <- endgenes[!duplicated(endgenes)]
  
  nfor <- normal[,conames]
  nfor <- as.data.frame(nfor)
  cfor <- cancer[,conames]
  cfor <- as.data.frame(cfor)
  Inormal <- cor(nfor)
  Icancer <- cor(cfor)
  Inormal[is.na(Inormal)] <- 1
  Icancer[is.na(Icancer)] <- 1
  Iscore <- Icancer-Inormal
  if(length(difnames)!=0){
    ainter <- rep(0,ncol(Iscore))
    ai <- matrix(rep(ainter,length(difnames)),ncol=length(difnames))
    binter <- rep(0,length(allgenes))
    bi <- matrix(rep(binter,length(difnames)),nrow=length(difnames))
    Iscore <- cbind(Iscore,ai)
    Iscore <- rbind(Iscore,bi)
    genenames <- c(conames,difnames)
    colnames(Iscore) <- genenames
    rownames(Iscore) <- genenames
    Iscore <- Iscore[allgenes,allgenes]
  }
  Iscore <- abs(Iscore)
  notends <- setdiff(allnodes,endgenes)
  g1<-igraph.from.graphNEL(mapkG2)
  dismatrix1<-distances(g1,mode="in")
  dismatrix1[which(dismatrix1=="Inf")]<-0
  dism <- dismatrix1
  colnames(dism) <- allnodes
  rownames(dism) <- allnodes
  dism <- dism[allgenes,allgenes]
  dismatrix1[which(dismatrix1!=0)] <- 1
  colnames(dismatrix1) <- allnodes
  rownames(dismatrix1) <- allnodes
  dismatrix1[notends,] <- 0
  allgenes <- as.character(allgenes)
  dismatrix1 <- dismatrix1[allgenes,allgenes] 
  
  edgesL <- Iscore*dismatrix1/dism
  edgesL[is.na(edgesL)] <- 0
  
  endscore <- apply(edgesL,1,sum)
  endscore <- endscore[which(endscore!=0)]
  endgenes <- setdiff(names(endscore),onenodes)
  endscore <- endscore[endgenes]
  endscore <- -sort(-endscore)
 
  
  genes <- factor(names(endscore),levels=names(endscore))
  scores <- data.frame(genes,perturbation=as.numeric(endscore))
  p<-ggplot(scores,aes(x=scores$genes,y=scores$perturbation))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+theme(panel.grid=element_blank())+theme(panel.grid.major = element_line(colour = NA))+theme(legend.position="none")
  return(p+xlab("Effector genes")+ylab("Signal variation received by effector genes"))
}

