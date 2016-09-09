# LINC FUNCTIONS

## 17 08 2016
############################################################
## ISSUES:
## native routines

.onAttach <- function(...) {
  
  packageStartupMessage(paste("This is LINC\n",
                              "(Co-Expression Analysis of lincRNAs)"))
    assign("enrichGO", value = signGO(), envir = globalenv())
}

# set a global binding for variables in ggplot functions
ASSIGNMENT <- NA; CORRELATION <- NA; EXPRESSION <- NA;
NO_SUBJECT <- NA; PC <- NA; SAMPLES <- NA
SUBJECT_1 <- NA; SUBJECT_2 <- NA; SUBJECT_3 <- NA
SUBJECT_4 <- NA; SUBJECT_5 <- NA; TERMS <- NA
VALUE <- NA; VARIANCE <- NA; lincRNA <- NA
protein_coding <- NA; pvalue <- NA; splot_1 <- NA
splot_10 <- NA; splot_2 <- NA; splot_3 <- NA
splot_4 <- NA; splot_5 <- NA; splot_6 <- NA
splot_7 <- NA; splot_8 <- NA; splot_9 <- NA
values <- NA; x <- NA; y <- NA

# HELPING FUNCTION "inlogical"
inlogical <- function(lg_promise, direct){
  if(class(lg_promise) != "logical" ){
    lg_promise <- direct
    warning(paste("argument interpreted as",
                  direct, "; TRUE/FALSE was expected "))
  } else {
    if(!any(lg_promise)) lg_promise <- FALSE
    if(any(lg_promise))  lg_promise <- TRUE 
  }
  return(lg_promise)
}

## HELPING FUNCTION "identifyGenes"
identifyGenes <-  function(gene){
  genepat <- c( "^[0-9]+",
                "ENSG.*",
                "ENST.*",
                "ENSP.*",
                "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",
                "Hs\\.([0-9]).*",
                "(NC|AC|NG|NT|NW|NZ|NM|NR|XM|XR|NP|AP|XP|YP|ZP)_[0-9]+")
  
  genesys <- c("entrezgene",
               "ensemblgene",
               "ensembltranscript",
               "ensemblprotein",
               "uniprot",
               "unigene",
               "refseq")
  
  grepres <- lapply(genepat, function(x, gene) {
    grep(x, gene, perl = TRUE) }, gene)  
  greptfv <- lapply(grepres, function(x){ length(x) > 0 })
  return(genesys[unlist(greptfv)])
}

## HELPING FUNCTION "modZscore"
modZscore <- function(x){
  zscore <- (0.6745*(x - median(x))/mad(x))
  if(any(abs(zscore) > 3.5)){
    x[which.max(abs(zscore))] <- NA
  }
  return(x)
}


## HELPING FUNCTION "svaSolv"
svaSolv <- function(y, mod, svs) {
  x <- cbind(mod, svs)
  intert = solve(t(x) %*% x) %*% t(x)
  beta = (intert %*% t(y)); rm(intert)
  n = seq_len(ncol(mod))
  res <- (y - t(as.matrix(x[, -n]) %*% beta[-n, ]))
  return(res)
}

## HELPING FUNCTION "callCor"

callCor <- function(corMethod, userFun, cor_use){
  if(is.null(userFun)){
    if(corMethod == "spearman"){
      return(Cppspear)
    }  
    if(corMethod != "spearman" &&
       corMethod != "user"){
      useFunction <- stats::cor
      method <- corMethod
      if(cor_use == "pairwise"){
        use <- "pairwise.complete.obs"
      } else {
        use <- "evertyhing"
      } 
      nonCppf <- function(useFunction, method,
                          use, ...){
        function(pcmatrix, ncmatrix, ...){
          x <- seq_len(nrow(pcmatrix));  lx <- length(x)
          y <- seq_len(nrow(ncmatrix)); ly <- length(y)          
          MM <- matrix(ncol = lx * ly , nrow = 2)
          MM[1,] <- rep(x, times = ly )
          MM[2,] <- as.vector(unlist(lapply(y,
                                            function(y, lx) {rep(y, times = lx)}, lx)))           
          cormatrix <- matrix(ncol = ly , nrow = lx )           
          for(i in seq_len(ncol(MM))){
            z1 <- MM[1,i]
            z2 <- MM[2,i]
            cormatrix[z1,z2] <- useFunction(
              pcmatrix[z1, ], ncmatrix[z2, ],
              method = method, use = use, ...)  
          } 
          return(cormatrix)
        } # function end inner
      } # function end    
      to_apply <- nonCppf(useFunction, method, use)   
      return(to_apply) 
    } 
  } else {
    # else: apply user function
    if(!is.function(match.fun(userFun))){
      stop("'userFun' is not a valid function")
    }
    useFunction <- match.fun(userFun)
    if(length(formals(useFunction)) != 2){
      warning(paste("'userFun' should have",
                    "two arguments 'x' and 'y'"))  
    }
    nonCppf2 <- function(useFunction){
      function(pcmatrix, ncmatrix, ...){
        x <- seq_len(nrow(pcmatrix));  lx <- length(x)
        y <- seq_len(nrow(ncmatrix)); ly <- length(y)                     
        MM <- matrix(ncol = lx * ly , nrow = 2)           
        MM[1,] <- rep(x, times = ly )
        MM[2,] <- as.vector(unlist(lapply(y,
                                          function(y, lx) {rep(y, times = lx)},lx)))           
        cormatrix <- matrix(ncol = ly , nrow = lx)           
        for(i in seq_len(ncol(MM))){
          z1 <- MM[1,i]
          z2 <- MM[2,i]
          cormatrix[z1,z2] <- useFunction(
            pcmatrix[z1, ],ncmatrix[z2, ])  
        } 
        return(cormatrix)
      } # function end inner
    } # function end          
    to_apply <- nonCppf2(useFunction)
    return(to_apply) 
  } 
}

## HELPING FUNCTION "correctESD"
correctESD <- function(x, alpha, rmax){
  canpos <- doesd(query = x, wquery = x, rmax, 
                  querypos = rep(0, rmax),
                  rmaxvec = rep(0, rmax),
                  mdiff = rep(0, length(x)),
                  alpha = alpha)
  maxtorm <- tail(canpos, 1)
  if(maxtorm > 0) x[canpos[seq_len(maxtorm)]] <- NA
  return(x)
}

## HELPING FUNCTION "detectesd"
detectesd <- function(x, alpha, rmax){
  canpos <- doesd(query = x, wquery = x, rmax, 
                  querypos = rep(0, rmax),
                  rmaxvec = rep(0, rmax),
                  mdiff = rep(0, length(x)),
                  alpha = alpha)
  maxtorm <- tail(canpos, 1)
  x <- FALSE
  if(maxtorm > 0) x <- TRUE
  return(x)
}

## HELPING FUNCTION "doCorTest"
doCorTest <- function(corMethod, cor_use){
  use <- cor_use; method <- corMethod
  useFunction <- stats::cor.test
  nonCppf3 <- function(useFunction, method,
                       use, ...){
    function(pcmatrix, ncmatrix, ...){
      x <- seq_len(nrow(pcmatrix));  lx <- length(x)
      y <- seq_len(nrow(ncmatrix)); ly <- length(y)          
      MM <- matrix(ncol = lx * ly , nrow = 2)
      MM[1,] <- rep(x, times = ly )
      MM[2,] <- as.vector(unlist(lapply(y,
                                        function(y, lx) {rep(y, times = lx)}, lx)))           
      cormatrix <- matrix(ncol = ly , nrow = lx )           
      for(i in seq_len(ncol(MM))){
        z1 <- MM[1,i]
        z2 <- MM[2,i]
        cormatrix[z1,z2] <- useFunction(
          pcmatrix[z1, ], ncmatrix[z2, ],
          method = method, use = use, 
          alternative = "greater", ...)$p.value  
      } 
      return(cormatrix)
    } # function end inner
  } # function end    
  to_apply <- nonCppf3(useFunction, method, use)   
  return(to_apply)
}  


# HELPING FUNCTION 'deriveCluster'
deriveCluster <- function(inmatrix, alpha){
  pc_genes <- rownames(inmatrix)
  tf_list <- apply(inmatrix, 2, function(x, pc_genes){
    sign_pc <- (x < alpha)
    if(any(sign_pc)){
      return(pc_genes[sign_pc])
    } else {
      return(NA)  
    }
  }, pc_genes)
  if(all(is.na(tf_list))){
    stop("No neighbour genes for this significance level")
  }
  # convert to a list
  if(is.matrix(tf_list)){
    tf_sp  <- split(tf_list, rep(seq_len(ncol(tf_list)),
                                 each = nrow(tf_list)))
    names(tf_sp) <- colnames(tf_list)
  } else {
    tf_sp <- tf_list
  }
  return(tf_sp)
}

# HELPING FUNCTION 'deriveCluster2'
deriveCluster2 <- function(inmatrix, n){
  pc_genes <- rownames(inmatrix)
  tf_list <- apply(inmatrix, 2, function(x, pc_genes){
    sign_pc <- sort(x, decreasing = FALSE)[seq_len(n)] 
    return(sign_pc)
  }, pc_genes)
  if(all(is.na(tf_list))){
    stop("No neighbour genes for this significance level")
  }
  # convert to a list
  if(is.matrix(tf_list)){
    tf_sp  <- split(tf_list, rep(seq_len(ncol(tf_list)),
                                 each = nrow(tf_list)))
    names(tf_sp) <- colnames(tf_list)
  }  else {
    tf_sp <- tf_list
  }
  return(tf_sp)
}

# Helping function 'cDiceDistance'
cDiceDistance <- function(inmatrix){
  size <- ncol(inmatrix)
  distmatrix <- matrix(data = rep(-1, times = size^2),
                       ncol = size, nrow = size)
  rownames(distmatrix) <- colnames(inmatrix)
  colnames(distmatrix) <- colnames(inmatrix)
  
  for(out in seq_len(size)){
    query <- inmatrix[, out]
    for(inn in seq_len(size)){
      subject <- inmatrix[, inn]
      
      INT <- length(which(is.element(
        (query), (subject))))
      ALL <- (length((query)) +
                length((subject)) - INT)
      
      CDD <- (ALL - INT)/(ALL + INT)
      if(!is.finite(CDD)) CDD <- 0
      if(any(is.na(query)) | any(is.na(subject))) CDD <- NA
      distmatrix[out, inn] <- CDD                   
    }    
  }
  return(distmatrix)
}

querycluster <- function(query = NULL,
                         queryTitle = NULL,
                         traits = NULL,
                         method = "spearman",
                         returnDat = FALSE,
                         mo_promise = NULL,
                         ...){
  errorm01 <- "Usage of argument 'query' is required"
  errorm02 <- "Input for 'traits' is invalid"
  errorm03 <- "Argument 'method' was 'NULL'" 
  errorm04 <- "'method' has to be 'spearman' or 'dicedist'"
  errorm05 <- "method 'dicedist' requires 'traits' > 2"
  errorm06 <- "min. two 'LINCcluster' objects are required"
  errorm07 <- "Not enough queries with valid traits"
  errorm08 <- paste("Computation for method 'spearman'",
                    "failed; method 'dicedist' may work instead")
  warnim01 <- "duplicated names of datasets"
  
  ## SECTION0: INPUT CONTROL
  if(is.null(query)) stop(errorm01)
  arg_lang <- match.call()
  arg_list <- lapply(arg_lang, as.character)
  if(!is.character(query)) stop(errorm01)
  if(!is.numeric(traits) && !is.null(traits)) traits <- NULL
  if(!is.character(method)) method <- "spearman"
  if(!is.character(queryTitle)) queryTitle <- query
  if(!is.logical(returnDat)) returnDat <- FALSE
  if(is.null(queryTitle)) query <- queryTitle
  
  # interpretation of 'traits'
  tr_promise <- traits
  if(length(tr_promise) != 1 &&
     !is.null(tr_promise)) stop(errorm02)
  if(!is.numeric(tr_promise) &&
     !is.null(tr_promise)) stop(errorm02)
  if(tr_promise < 3 &&
     !is.null(tr_promise)) stop(errorm02)
  
  # interpretation of 'method'
  if(is.null(method)) stop(errorm03)
  mh_promise <- try(match.arg(method,
                              c("spearman", "dicedist")), silent = TRUE)
  if(class(mh_promise) == "try-error") stop(errorm04)  
  if((mh_promise == "dicedist") && is.null(tr_promise))
    stop(errorm05) 
  
  # interpretation of 'returnDat'
  returnDat <- inlogical(returnDat, direct = FALSE) 
  
  ## SECTION1: OBJECTS FROM CALL
  # capture all 'LINCcluster' objects
  if(class(mo_promise) != 'list'){
    arg_list <- as.character(match.call())    
    mo_promise <- mget(unlist(arg_list)[-1], envir =
                         globalenv(), ifnotfound = NA)
  }
  
  class_promise <- vapply(mo_promise, function(x){
    if(class(x) == "LINCcluster" ){
      TRUE
    } else {
      FALSE
    }
  }, FALSE)
  if(length(which(class_promise)) < 2) stop(errorm06)
  linc_class <- mo_promise[class_promise]   
  
  ## SECTION2: TRAITS
  # derive the traits
  trait_list <- lapply(linc_class, function(x){
    qn_promise  <- names(
      results(linCenvir(x)$cluster)[[1]]$neighbours)
    
    tt_promise <- try(index <- is.element(qn_promise, query),
                      silent = TRUE)
    if(class(tt_promise) == "try-error") tt_promise <- NULL
    if(!is.null(tt_promise)){
      y <- results(linCenvir(x)$cluster)[[1]]$neighbours[index] 
      if(is.list(y) && length(y) > 1) y <- NULL
    } else {
      y <- NULL
    }
    if(is.null(y) | anyNA(y)) {
      return(NULL) 
    } else{
      if(is.null(tr_promise)){
        return(unlist(y))
      } else {
        return(unlist(y)[seq_len(tr_promise)])
      }
    }
  })
  
  # missing traits, index
  na_in <- !unlist(lapply(trait_list, is.null))
  if(length(which(na_in)) < 2) stop(errorm07)
  
  # get custom IDs and colours
  c_name <- vector(); c_color <- vector()
  for(n in seq_len(length(linc_class))){
    if(exists("customID", history(linc_class[[n]]))){
      c_name[n] <- paste(n, history(linc_class[[n]])$customID,
                         sep = "_") 
      c_color[n] <- history(linc_class[[n]])$customCol
    } else {
      c_name[n] <- paste(n, "COND", sep = "_")
      c_color[n] <- "black"
    }
  }
  if(!identical(c_name, unique(c_name))){
    warning(warnim01)   
  } 
  
  # correct for missing traits
  trait_list <- trait_list[na_in]
  c_name <- c_name[na_in]
  c_color <- c_color[na_in]
  
  ## SECTION3: DISTANCE METRIC
  # do the "spearman" method
  if(mh_promise == "spearman"){
    cp_union <- unique(unlist(trait_list)) 
    union_matrix <- matrix(NA, nrow = length(cp_union),
                           ncol = length(trait_list))   
    rownames(union_matrix) <- cp_union; n = 1
    for(m in seq_len(length(linc_class))[na_in]){
      x <- correlation(linCenvir(linc_class[[m]])$linc)$cormatrix
      rest_union <- intersect(cp_union, rownames(x))
      spear_value <- x[rest_union, as.character(query)]
      union_matrix[rest_union, n] <- spear_value
      n = n + 1
    }
    union_cor <- try(callCor("spearman", NULL,"pairwise")(
      t(union_matrix), t(union_matrix)), silent = TRUE)
    if(class(union_cor) == "try-error" |
       base::anyNA(union_cor)) stop(errorm08)
    crude_dist <- (1 - union_cor)
  }
  
  # do the 'dicedist' method
  if(mh_promise == "dicedist"){
    c_matrix <- matrix(NA, ncol = length(trait_list),
                       nrow = tr_promise)
    for(k in seq_len(length(trait_list))){
      c_matrix[,k] <- trait_list[[k]]
    }
    crude_dist <- cDiceDistance(c_matrix)
  }
  
  # distmatrix and dendrogram
  colnames(crude_dist) <- c_name
  rownames(crude_dist) <- c_name
  distmat <- as.dist(crude_dist)
  dist_data <- as.data.frame(crude_dist)
  colnames(dist_data) <- paste((seq_len(length(trait_list))),
                               "_", sep = "")
  rownames(dist_data) <- c_name
  dist_clust <- hclust(distmat, method = "average")
  #suppressPackageStartupMessages(require("ape")) 
  
  dist_phylo <- as.phylo(dist_clust)
  dist_phylo$tip.label <- colnames(crude_dist)
  
  #shared interaction partners
  m_len <- length(trait_list)
  si_promise <- matrix(0, ncol = m_len, nrow = m_len)
  for(s in seq_len(m_len)){
    for(i in seq_len(m_len)){
      if(si_promise[i,s] == 0){
        partners <- length(intersect(trait_list[[i]],
                                     trait_list[[s]]))
        si_promise[s,i] <- partners
      }
    }
  }
  n_partner <- as.data.frame(si_promise)
  rownames(n_partner) <- colnames(crude_dist)
  colnames(n_partner) <- colnames(crude_dist)
  
  ## SECTION4: PLOTTING
  # plot the cluster
  tree <- ggtree(dist_phylo, colour = "dodgerblue4",
                 layout= "rectangular", alpha = 0.5, size= 1 ) +
    coord_flip() + scale_x_reverse() +
    geom_tiplab(size = 3.5, angle = -90,
                colour = c_color, hjust = 0)
  
  clust_heat <- gheatmap(tree, n_partner, offset = 0.3,
                         width = 1.2, colnames = TRUE,
                         colnames_position = "top",
                         low = "white", high = "black")
  
  plot_tree <- (tree +  ggtitle(queryTitle) + theme(
    plot.title = element_text(face = "bold", size = 25)))
  
  querycluster_img <- readPNG(system.file("extdata", "querycluster_img.png",
                                          package ="LINC"))
  querycluster <- rasterGrob(querycluster_img, interpolate = TRUE)
  
  plot_it <- grid.arrange(querycluster, clust_heat, ncol = 2)
  # return
  if(returnDat){
    return(list(cluster = dist_phylo,
                distmatrix = distmat))
  } else {
    return(invisible(plot_it))
  }
} # function end  


# Helping function listQuery
setGeneric(name = "listQuery",
           def = function(
             input){
             standardGeneric("listQuery")
           })
setMethod(f   = "listQuery",
          signature = c("LINCmatrix"),  
          def = function(
            input){
            return(colnames(correlation(linCenvir(input)$linc))[[1]])
          })
setMethod(f   = "listQuery",
          signature = c("LINCcluster"),  
          def = function(
            input){
            return( names(results(linCenvir(input)$cluster)[[1]]$neighbours)  ) 
          })
setMethod(f   = "listQuery",
          signature = c("LINCbio"),  
          def = function(
            input){
            return( names(results(linCenvir(input)$bio)[[1]]) ) 
          })
setMethod(f   = "listQuery",
          signature = c("LINCsingle"),  
          def = function(
            input){
            return(results(linCenvir(input)$single)$query ) 
          })

# function getcoexpr
getcoexpr<- function(input, query = NULL){
  
  # check the class          
  if(!any(is.element(class(input), 
                     c("LINCmatrix", "LINCcluster",
                       "LINCsingle", "LINCbio")))){
    stop("input is not of a supported 'LINC' class")
  } 
  
  if(class(input) == "LINCcluster"){
    ip_promise  <- try(results(input)$cluster$neighbours[
      names(results(input)$cluster$neighbours) == query],
      silent = TRUE)
    if(class(ip_promise) != "try-error"){
      ip_promise <- unlist(ip_promise)
      names(ip_promise) <- NULL
      return(ip_promise)
    }
  }
  
  if(class(input) == "LINCsingle"){
    ip_promise  <- try(names(results(input)$cor),
                       silent = TRUE)
    if(class(ip_promise) != "try-error"){
      return(ip_promise)
    }
  }
}


