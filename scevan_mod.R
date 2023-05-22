

pipelineCNA1 <- function(count_mtx, sample="", par_cores = 20, norm_cell = NULL, SUBCLONES = TRUE, beta_vega = 0.5, ClonalCN = TRUE, plotTree = TRUE, AdditionalGeneSets = NULL, SCEVANsignatures = TRUE, organism = "human"){
  
  dir.create(file.path("./output"), showWarnings = FALSE)
  
  start_time <- Sys.time()
  
  normalNotKnown <- length(norm_cell)==0
  
  res_proc <- SCEVAN:::preprocessingMtx(count_mtx,sample, par_cores=par_cores, findConfident = normalNotKnown, AdditionalGeneSets = AdditionalGeneSets, SCEVANsignatures = SCEVANsignatures, organism = organism)
  
  if(normalNotKnown) norm_cell <- names(res_proc$norm_cell)

  res_class <- classifyTumorCells1(res_proc$count_mtx_norm, res_proc$count_mtx_annot, sample, par_cores=par_cores, ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE, beta_vega = beta_vega)
  
  print(paste("found", length(res_class$tum_cells), "tumor cells"))
  classDf <- data.frame(class = rep("filtered", length(colnames(res_proc$count_mtx))), row.names = colnames(res_proc$count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  
  end_time<- Sys.time()
  print(paste("time classify tumor cells: ", end_time -start_time))

  return(classDf)
}


classifyTumorCells1 <- function(count_mtx, annot_mtx, sample = "", distance="euclidean", par_cores=20, ground_truth = NULL, norm_cell_names = NULL, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE, beta_vega = 0.5){
  
  set.seed(1)
  
  if (length(norm_cell_names) < 1){
    print("7) Measuring baselines (pure tumor - synthetic normal cells)")
    count_mtx_relat <- SCEVAN:::removeSyntheticBaseline(count_mtx, par_cores=par_cores)
    
  } else {
    print("7) Measuring baselines (confident normal cells)")
    
    if(length(norm_cell_names) == 1){
      basel <- count_mtx[, which(colnames(count_mtx) %in% norm_cell_names)]
    }else{
      basel <- apply(count_mtx[, which(colnames(count_mtx) %in% norm_cell_names)],1,median)
    }
    
    ##relative expression using pred.normal cells
    count_mtx_relat <- count_mtx-basel
    
  }
  
  
  ##### smooth data ##### 
  if(SMOOTH){
    print("8) Smoothing data")
    
    niters=100
    alpha=0.5
    DeltaT = 0.2
    
    nonLinSmooth <- function(c){
      y <- count_mtx_relat[, c]
      y <- append(0,y)
      for(i in 1:niters){
        DeltaP = (y[2:length(y)]-y[1:(length(y)-1)])/alpha
        DeltaP <- c(DeltaP,0)
        tD <- tanh(DeltaP)
        DeltaM <- tD[2:length(y)]-tD[1:(length(y)-1)]
        DeltaM <- c(0,DeltaM)
        y <- y + DeltaT*DeltaM
      }
      y <- y[2:length(y)]
      return(y)
    } 
    
    #test.mc <- parallel::mclapply(1:ncol(count_mtx_relat), nonLinSmooth, mc.cores = par_cores)
    
    if(Sys.info()["sysname"]=="Windows"){
      cl <- parallel::makeCluster(getOption("cl.cores", par_cores))
      test.mc <- parallel::parLapply(cl, 1:ncol(count_mtx_relat), nonLinSmooth)
      parallel::stopCluster(cl)
    }else{
      test.mc <- parallel::mclapply(1:ncol(count_mtx_relat), nonLinSmooth, mc.cores = par_cores)
    }
    
    count_mtx_smooth <- matrix(unlist(test.mc), ncol = ncol(count_mtx_relat), byrow = FALSE)
    rm(test.mc)
    colnames(count_mtx_smooth) <- colnames(count_mtx_relat)
    count_mtx_relat <- count_mtx_smooth
   }
    
  ##### Segmentation with VegaMC #####
  
  if(SEGMENTATION_CLASS & length(norm_cell_names) > 0){ #){
    
    print("9) Segmentation (VegaMC)")
    
    mtx_vega <- cbind(annot_mtx[,c(4,1,3)], count_mtx_relat)
    colnames(mtx_vega)[1:3] <- c("Name","Chr","Position")
    
    breaks <- getBreaksVegaMC(mtx_vega, annot_mtx[, 3], sample, beta_vega)
    
    subSegm <- read.csv(paste0("./output/ ",sample," vega_output"), sep = "\t")
    
    segmAlt <- abs(subSegm$Mean)>0.05 | (subSegm$G.pv<0.01 | subSegm$L.pv<0.01)
    #segmAlt <- append(segmAlt[1],append(segmAlt, segmAlt[length(segmAlt)]))
    
    CNA_mtx <- SCEVAN:::computeCNAmtx(count_mtx_relat, breaks, par_cores, segmAlt)
    
    SEGM <- TRUE
    
  }else{
    SEGM <- FALSE
    CNA_mtx <- count_mtx_relat
  }
  
  colnames(CNA_mtx) <- colnames(count_mtx_relat)
  CNA_mtx <- apply(CNA_mtx,2, function(x)(x <- x-mean(x)))
  
  print("10) Adjust baseline")
  
  if(length(norm_cell_names) < 1){
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx, method = distance)), method = "ward.D")
    }
    
    #plot heatmap
    print("11) plot heatmap")
    
    #plotCNA(annot_mtx$seqnames, CNA_mtx, hcc, sample)
    #save(CNA_mtx, file = paste0("./output/",sample,"_CNAmtx.RData"))

  } else {
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx, method = distance)), method = "ward.D")
    }
    
    hcc2 <- cutree(hcc,2)
    names(hcc2) <- colnames(CNA_mtx)
    
    cellType_pred <- SCEVAN:::classifyCluster(hcc2, norm_cell_names)
    
    ################removed baseline adjustment
    CNA_mtx_relat <- CNA_mtx-apply(CNA_mtx[,which(cellType_pred=="non malignant")], 1, mean)
    CNA_mtx_relat <- apply(CNA_mtx_relat,2,function(x)(x <- x-mean(x)))
    CNA_mtx_relat <- CNA_mtx_relat/(0.5*(max(CNA_mtx_relat)-min(CNA_mtx_relat)))
    
    #CNA_mtx_relat[abs(CNA_mtx_relat)<0.05] <- 0
    #CNA_mtx_relat <- (CNA_mtx_relat*4)^3
    
    if(SEGM){
      count_mtx_relat <- count_mtx_relat-apply(count_mtx_relat[,which(cellType_pred=="non malignant")], 1, mean)
      count_mtx_relat <- apply(count_mtx_relat,2,function(x)(x <- x-mean(x)))
      count_mtx_relat <- count_mtx_relat/(0.5*(max(count_mtx_relat)-min(count_mtx_relat)))
    }
    
    if(distance=="euclidean"){
      hcc <- hclust(parallelDist::parDist(t(CNA_mtx_relat),threads =par_cores, method = distance), method = "ward.D")
    }else {
      hcc <- hclust(as.dist(1-cor(CNA_mtx_relat, method = distance)), method = "ward.D")
    }
    
    hcc2 <- cutree(hcc,2)
    names(hcc2) <- colnames(CNA_mtx_relat)
    
    cellType_pred <- classifyCluster(hcc2, norm_cell_names)
    
    res <- cbind(names(cellType_pred), cellType_pred)
    colnames(res) <- c("cell.names", "pred")
    
    print("11) plot heatmap")
    
    #plotCNA(annot_mtx$seqnames, CNA_mtx_relat, hcc, sample, cellType_pred, ground_truth)
    
    #save(CNA_mtx_relat, file = paste0("./output/",sample,"_CNAmtx.RData"))
    
  }
  if(length(norm_cell_names) < 1){
    tum_cells <- colnames(CNA_mtx)
    tum_cells <- gsub("\\.","-",tum_cells)
  }else{
    tum_cells <- names(res[,2][res[,2] == "malignant"])
  }
  
  if(SEGM){
    ress <- list(tum_cells, cbind(annot_mtx[,c(4,1,3)], count_mtx_relat), norm_cell_names)
  }else if(length(norm_cell_names) < 1){
    ress <- list(tum_cells, cbind(annot_mtx[,c(4,1,3)], CNA_mtx), norm_cell_names)
  }else{
    ress <- list(tum_cells, cbind(annot_mtx[,c(4,1,3)], CNA_mtx_relat), norm_cell_names)
  }
  
  names(ress) <- c("tum_cells", "CNAmat", "confidentNormal")
  return(ress)
}


classifyCluster <- function(hcc2, norm_cell_names){
  
  perc_norm <- length(intersect(names(hcc2[hcc2==1]), norm_cell_names))/length(hcc2[hcc2==1])
  perc_norm <- c(perc_norm,length(intersect(names(hcc2[hcc2==2]), norm_cell_names))/length(hcc2[hcc2==2]))
  clust_norm <- which(perc_norm==max(perc_norm))
  
  cellType_pred <- names(hcc2)
  cellType_pred[hcc2 == clust_norm] <- "non malignant"
  cellType_pred[!hcc2 == clust_norm] <- "malignant"
  names(cellType_pred) <- names(hcc2)
  
  return(cellType_pred)
}







