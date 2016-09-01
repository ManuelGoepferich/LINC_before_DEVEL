## LINC_METHODS

## getter methods
setGeneric("results", function(x) standardGeneric("results"))
setMethod("results", "LINCmatrix", function(x) x@results)

setGeneric("assignment", function(x) standardGeneric("assignment"))
setMethod("assignment", "LINCmatrix", function(x) x@assignment)

setGeneric("correlation", function(x) standardGeneric("correlation"))
setMethod("correlation", "LINCmatrix", function(x) x@correlation)

setGeneric("express", function(x) standardGeneric("express"))
setMethod("express", "LINCmatrix", function(x) x@expression)

setGeneric("history", function(x) standardGeneric("history"))
setMethod("history", "LINCmatrix", function(x) x@history)

setGeneric("linCenvir", function(x) standardGeneric("linCenvir"))
setMethod("linCenvir", "LINCmatrix", function(x) x@linCenvir)

## setter methods
setGeneric("results<-", function(x, value) standardGeneric("results<-"))
setReplaceMethod("results", "LINCmatrix",
          function(x, value) {x@results <- value; validObject(x); x})

setGeneric("assignment<-", function(x, value) standardGeneric("assignment<-"))
setReplaceMethod("assignment", "LINCmatrix",
                 function(x, value) {x@assignment <- value; validObject(x); x})

setGeneric("correlation<-", function(x, value) standardGeneric("correlation<-"))
setReplaceMethod("correlation", "LINCmatrix",
                 function(x, value) {x@correlation <- value; validObject(x); x})

setGeneric("express<-", function(x, value) standardGeneric("express<-"))
setReplaceMethod("express", "LINCmatrix",
                 function(x, value) {x@expression <- value; validObject(x); x})

setGeneric("history<-", function(x, value) standardGeneric("history<-"))
setReplaceMethod("history", "LINCmatrix",
                 function(x, value) {x@history <- value; validObject(x); x})

setGeneric("linCenvir<-", function(x, value) standardGeneric("linCenvir<-"))
setReplaceMethod("linCenvir", "LINCmatrix",
                 function(x, value) {x@linCenvir <- value; validObject(x); x})


## GENERIC FUNCTION "justlinc"
setGeneric(name = "justlinc",
           def  = function(object,
                           targetGenes = "lincRNA",
                           rmPC        = TRUE
           ){
             standardGeneric("justlinc")
           })    
setMethod(f = "justlinc",
          signature = c("matrix", "ANY", "ANY"),   
          definition  =   function(object,
                                   targetGenes = "lincRNA",
                                   rmPC   = TRUE
          ){
            ## method for a matrix as input
            ## PREDEFINITIONs
            #data(ENSG_BIO, package = "LINC")
            #data(ENTREZ_BIO, package = "LINC")          
            #data(ENSG_PC, package = "LINC")
            
            ENSG_BIO_DIR <- system.file("extdata", "ENSG_BIO.RData", package = "LINC")
            load(ENSG_BIO_DIR)
            ENTREZ_BIO_DIR <- system.file("extdata", "ENTREZ_BIO.RData", package = "LINC")
            load(ENTREZ_BIO_DIR)
            ENSG_PC_DIR <- system.file("extdata", "ENSG_PC.RData", package = "LINC") 
            load(ENSG_PC_DIR)
            
            
            errorm00 <- "'justlinc' failed, please use 'linc'"
            errorm01 <- "no match for 'targetGenes', valid targets:"
            
            errorm03 <- "Only one gene biotype is allowed"
            errorm05 <- "Gene names as 'rownames(object)' required"
            errorm06 <- "No Ensembl or Entrez ids, method failed"
            errorm07 <- "'targetGenes' supplied as gene ids not found"
            
            warnim00 <- paste("correlation matrix contains infinite",
                              "or missing values; converted to 0")
            warnim01 <- "This method is intended for Ensembl gene ids"
            warnim02 <- paste("long computation times expected for",
                              "> 100 'targetGenes'")
            
            exit_print <- names(table(ENSG_BIO$gene_biotype))
            on.exit(print(exit_print))
            
            ## SECTION0: INPUT CHECK
            # targetGenes
            tg_promise <- try(is.element(targetGenes,
                                         c(ENSG_BIO$gene_biotype, rownames(object))),
                              silent = TRUE)
            
            if(class(tg_promise) == "try-error" |
               !all(tg_promise)) stop(errorm01) 
            
            # suggest gene systems
            gN_promise <- rownames(object)
            if(is.null(gN_promise)) stop(errorm05)
            
            gD_promise <- try(identifyGenes(gN_promise),
                              silent = TRUE)
            if(class(gD_promise) == "try-error" |
               length(gD_promise) == 0) stop(errorm00)
            if(!any(is.element(gD_promise, "ensemblgene")))
              warning(warnim01)
            if(!any(is.element(gD_promise, c("ensemblgene",
                                             "entrezgene")))) stop(errorm06)
            
            # removal of PC 
            rm_promise <- inlogical(rmPC, direct = TRUE)
            
            ## SECTION1: MATRIX SEPARATION    
            # remove PCs (> 30% of main variance)
            if(rm_promise){
              pca_object <- prcomp(object, center = FALSE,
                                   scale. = FALSE)  
              mvar30 <- sum(pca_object$sdev) * 0.3
              
              n = 1; mvar_sum = 0
              while(mvar_sum < mvar30){
                mvar_sum <- mvar_sum + pca_object$sdev[n]
                n = n + 1
              }
              object <- pca_object$x[,n:ncol(object)] %*% t(
                pca_object$rotation[,n:ncol(object)])  
            }
            
            # separation
            if(any(is.element(gD_promise, "ensemblgene"))){
              in_match <- match(ENSG_PC$ENSG, rownames(object))
              pc_matrix  <- object[in_match[!is.na(in_match)], ]
              rownames(pc_matrix) <- ENSG_PC$ENTREZ
            } else {
              pc_matrix  <- object[is.element(
                rownames(object), ENSG_PC$ENTREZ), ]
            }
            
            # reduce samples and compute the correlation matrix
            if(ncol(object) > 30){
              samples <- seq_len(30)
            } else {
              samples <- seq_len(ncol(object))  
            }
            
            # no queries are supplied 
            if(any(is.element(targetGenes, ENSG_BIO$gene_biotype))){
              if(length(targetGenes) != 1) stop(errorm03) 
              
              if(any(is.element(gD_promise, "ensemblgene"))){
                nc_genes  <- ENSG_BIO$ensembl_gene_id[is.element(
                  ENSG_BIO$gene_biotype, targetGenes)]
              } else {
                e_index <- is.element(ENTREZ_BIO$entrez_gene_id,
                                      rownames(object))
                t_index <- is.element(ENTREZ_BIO$gene_biotype,
                                      targetGenes)
                et_index <- ((e_index + t_index) == 2)
                nc_genes <- as.character(ENTREZ_BIO$entrez_gene_id[
                  et_index]) 
              }
              nc_matrix  <-  object[nc_genes, ]
              if(nrow(pc_matrix) < 5000) stop(errorm00)
              if(nrow(nc_matrix) < 500) stop(errorm00)
              
              # select for median expression and variance
              pc_median <- median(pc_matrix, na.rm = TRUE)
              nc_median <- median(nc_matrix, na.rm = TRUE)
              pc_rw_median <- rowMedians(pc_matrix, na.rm = TRUE)
              nc_rw_median <- rowMedians(nc_matrix, na.rm = TRUE)
              pc_matrix <- pc_matrix[(pc_rw_median > 0.25*pc_median), ]
              
              if(nrow(pc_matrix) < 5000) stop(errorm00)
              pc_var <- apply(pc_matrix, 1, var)
              pc_matrix <- pc_matrix[order(pc_var,
                                           decreasing = TRUE)[1:5000], ]
              
              if(nrow(nc_matrix) < 10) stop(errorm00)
              nc_matrix <- nc_matrix[(nc_rw_median > 0.25*nc_median), ]        
              nc_var <- apply(nc_matrix, 1, var)
              nc_matrix <- nc_matrix[order(nc_var,
                                           decreasing = TRUE)[1:100], ]
              
              cormatrix <- try(Cppspear(pc_matrix[, samples ],
                                        nc_matrix[, samples ]), silent = TRUE)
              if(class(cormatrix) == "try-error") stop(errorm00)
              rownames(cormatrix) <- rownames(pc_matrix)
              colnames(cormatrix) <- rownames(nc_matrix)
              if(!all(is.finite(cormatrix))){
                warning(warnim00)
                cormatrix[is.infinite(cormatrix) |
                            is.na(cormatrix) ] <- 0
              }
              
              # select 10 genes with best correlations
              max_cor <- apply(cormatrix, 2, max)
              q10_genes <- order(max_cor, decreasing = TRUE)[1:10]
              th_matrix <- (cormatrix[,q10_genes] > 0.75)
              
              pc_list <- list()
              q10_list <- list()
              for(q in 1:10){
                i_partners <- cormatrix[th_matrix[,q], q10_genes[q]]
                i_partners <- i_partners[order(i_partners,
                                               decreasing = TRUE)][1:50]
                i_partners <- i_partners[!is.na(i_partners)] 
                q10_list[[q]] <- i_partners
                if(length(i_partners) > 25){
                  pc_list[[q]] <- names(i_partners)
                } else {
                  pc_list[[q]] <- NULL
                }
              }
              names(pc_list) <- colnames(th_matrix)
              names(q10_list) <- colnames(th_matrix)
              null_index <- vapply(pc_list, is.null, TRUE)
              
              if(length(which(!null_index)) < 10) stop(errorm00)
              pc_list <- pc_list[!null_index]
              pc_path <- try(compareCluster(pc_list, fun =
                                              "enrichPathway"), silent = TRUE)
              if(class(pc_path) == "try-error") stop(errorm00)
              
              
              # a <- compareCluster(pc_list)
              
              ## SECTION 3: PLOTTING
              
              # plot expression and correlation of best
              for(pp in 1:10){
                subject <- pc_matrix[pc_list[[pp]][1], ]
                query <- nc_matrix[names(pc_list)[pp], ]
                s_df <- data.frame(lincRNA = query, protein_coding
                                   = subject, SAMPLES = seq_len(ncol(pc_matrix)))
                s_cor <- cor(subject, query, method = "spearman")
                splot <- ggplot(s_df) + geom_point(aes(x =
                                                         protein_coding, y = lincRNA), colour = "firebrick4",
                                                   size = 2) + ggtitle(paste(names(pc_list)[pp], "vs", 
                                                                             pc_list[[pp]][1],"| CORRELATION:", round(s_cor, 2))) +
                  theme(title = element_text(size = 12, face = "bold")) 
                assign(paste("splot", pp, sep = '_'), splot )
              }
              grid.arrange( splot_1, splot_6,
                            splot_2, splot_7,
                            splot_3, splot_8,
                            splot_4, splot_9,
                            splot_5, splot_10, 
                            ncol =2, nrow = 5, top = textGrob(
                              "PLOT 1: EXPRESSION & CORRELATION (DETAILS)",
                              gp = gpar(fontsize = 30, font = 2, face = "bold")))
              
              # Plot for enriched pathways
              cplot <- plot(pc_path, showCategory = 10)
              grid.arrange(cplot, top = textGrob(
                "PLOT 2: ENRICHED PATHWAYS",
                gp = gpar(fontsize = 30, font = 2, face = "bold")))
              
              final_list <- list(REACTOMEPA = pc_path,
                                 INTERACTION_PARTNERS = q10_list)
            } else {
              
              sf_targets <- is.element(rownames(object), targetGenes)        
              if(!any(sf_targets)) stop(errorm07)
              if(length(targetGenes) > 100) warning(warnim02) 
              
              if(length(targetGenes) > 1){
                nc_matrix <- object[targetGenes, ]
              } else {
                nc_matrix <- rbind(object[targetGenes, , drop = FALSE],
                                   object[targetGenes, , drop = FALSE])
                rownames(nc_matrix) <- c(targetGenes, paste(targetGenes,
                                                            "1", sep = '_'))
              }
              
              # select for median expression and variance
              pc_median <- median(pc_matrix, na.rm = TRUE)
              pc_rw_median <- rowMedians(pc_matrix, na.rm = TRUE)
              pc_matrix <- pc_matrix[(pc_rw_median > 0.25*pc_median), ]
              
              if(nrow(pc_matrix) < 5000) stop(errorm00)
              pc_var <- apply(pc_matrix, 1, var)
              pc_matrix <- pc_matrix[order(pc_var,
                                           decreasing = TRUE)[1:5000], ]
              
              c_matrix <- rbind(pc_matrix[ ,samples],
                                nc_matrix[ ,samples])
              l_matrix <- linc(object = c_matrix, codingGenes =
                                 is.element(rownames(c_matrix),
                                            rownames(pc_matrix)), verbose = FALSE)
              if(length(targetGenes) > 1){
                l_cluster <- clusterlinc(l_matrix,
                                         pvalCutOff = 0.0005,
                                         verbose = FALSE)
                final_list <- getbio(l_cluster, annotateFrom =
                                       "enrichPathway",   
                                     translate = "none", verbose = FALSE)
                #to_return <- l_bio
                plotlinc(final_list)
              } else {
                final_list <- singlelinc(l_matrix, query = targetGenes,
                                         threshold = 0.00005, annotateFrom =
                                           "enrichPathway", verbose = FALSE)
                plotlinc(final_list)
              }
            }
            print <- function(x = final_list){ return(x) }
          })         

## GENERIC FUNCTION "linc"
## create a LINCmatrix instance
setGeneric(name = "linc",
           def  = function(object,
                           codingGenes = NULL,
                           corMethod   = "spearman",
                           batchGroups = NULL,
                           nsv         = 1,
                           rmPC        = NULL,
                           outlier     = NULL,
                           userFun     = NULL,
                           verbose     = TRUE){
             standardGeneric("linc")
           })    
setMethod(f = "linc",
          signature = c("matrix"),   
          definition  = function(object,
                                 codingGenes = NULL,
                                 corMethod   = "spearman",
                                 batchGroups = NULL,
                                 nsv         = 1,
                                 rmPC        = NULL,
                                 outlier     = NULL,
                                 userFun     = NULL,
                                 verbose     = TRUE){
            ## method for a matrix as input
            ## PREDEFINITIONs
            
            # errors, warnings and messages  
            errorm00 <- paste("Assignment of protein-coding genes",
                              "in 'codingGenes' is required")  
            errorm01 <- paste("'codingGenes' must have the same",
                              "length as 'nrow(object)'")  
            errorm02 <- paste("'corMethod' needs to be 'pearson',",
                              "'kendall' or 'spearman'")  
            errorm03 <- "A numeric matrix is required as input"  
            errorm04 <- "Size or content of matrix insufficient"  
            errorm05 <- "Gene names as 'rownames(object)' required" 
            errorm06 <- paste("'batchGroups' need to be of the",
                              "same length as the columns")  
            errorm07 <- paste("Not allowed to use the same name", 
                              "for every entry in 'batchGroups'")  
            errorm08 <- paste("unable to use 'rmPC' as an index", 
                              "vector for the removal of pcs")
            errorm09 <- paste("'outlier' needs to be 'zscore',",
                              "or 'esd'")  
            errorm10 <- paste("'codingGenes' needs to be a gene",
                              "annotation or a logical vector") 
            errorm11 <- paste("Error in argument 'codingGenes',",
                              "not enough protein-coding genes")
            errorm12 <- paste("unable to compute correlation matrix:",
                              "1. check input for infinite values / NAs",
                              "2. check user-defined correlation function", sep = '\n') 
            errorm13 <- "computation of cor. matrix lnc vs lnc failed"
            warnim01 <- "Input 'object' contains infinite values"
            warnim02 <- "'linc' was unable to identify a gene system"
            warnim03 <- paste(
              "single outliers and high sample variance were detected",
              "by ESD and ANOVA; statistical correction is recommended",
              sep = '\n')  
            warnim04 <- paste("Subsequent use of sva and removal of",
                              "principle components is not intended")
            warnim05 <- paste("correlation matrix contains infinite",
                              "or missing values; converted to 0")
            inform01 <- quote(paste("linc: removed ", infrm, 
                                    "rows contaning only infinite values"))
            inform02 <- quote(paste("removed", length(obvar[obvar
                                                            == 0]), "zero variance genes from input"))
            inform22 <- "removed genes with duplicated names"
            inform03 <- "linc: gene system(s) assumed:"
            inform04 <- "linc: correction by sva was called"  
            inform05 <- "linc: remove principle components"
            inform06 <- quote(paste("linc: The outlier method '",
                                    ol_promise, "' was called"))  
            inform07 <- quote(paste("linc: Correlation function",
                                    " with '", cor_use,  "' called", sep = ''))
            inform08 <- paste("linc: Computation of correlation",
                              "matrix started")
            
            # environments and object
            store <- new.env(parent = emptyenv())
            out_history <- new.env(parent = emptyenv())
            
            ## SECTION0: INPUT CONTROL
            # use of the 'codingGenes' argument  
            if(is.null(codingGenes)) stop(errorm00)
            if(length(codingGenes) != nrow(object)) stop(errorm01)
            pc_promise <- codingGenes
            
            # interpretation of'verbose'
            if(class(verbose) != "logical" ){
              verbose <- TRUE
            } else {
              if(!any(verbose)) verbose <- FALSE
              if(any(verbose))  verbose <- TRUE
            }
            if(!verbose) message <- function(x) x
            
            # interpretation of the 'corMethod' argument
            cM_promise <- try(match.arg(corMethod,
                                        c("pearson", "kendall", "spearman")),
                              silent = TRUE)
            if(class(cM_promise) == "try-error") stop(errorm02)
            if(!is.null(userFun)) cor_Method <- "user-defined"
            
            # evaluation of 'object' and its gene ids
            # numeric matrix
            if(!is.numeric(object)) stop(errorm03)
            
            # only infinite values
            if(!all(is.finite(object))){ 
              warning(warnim01)
              mobject <- object[apply(object, 1, function(x){
                any(is.finite(x)) }), ]
              pcobject <- object; rownames(pcobject) <- pc_promise
              pcobject <- pcobject[apply(pcobject, 1, function(x){
                any(is.finite(x)) }), ]
              infrm  <- nrow(object) - nrow(mobject)
              if(infrm != 0){ message(inform01); object <- mobject  
              pc_promise <- rownames(pcobject)
              }
            }
            
            # 0 variance rows
            obvar <- apply(object, 1, var)
            if(is.element(0, obvar)){
              object <- object[obvar != 0, ]
              pc_promise <- pc_promise[obvar != 0]
              message(eval(inform02))
            }
            
            # rows duplicated
            if(any(duplicated(rownames(object)))){
              pc_promise <- pc_promise[!duplicated(rownames(object))]
              object <- object[(!duplicated(rownames(object))), ]
              message(inform22)
            }
            out_object <- object
            
            object <- object[!is.na(rownames(object)),]
            pc_promise <- pc_promise[!is.na(pc_promise)]
            
            # is the matrix to small
            if(!all(dim(object) > 5)) stop(errorm04)
            colnum <- ncol(object)
            
            # suggest gene systems
            gN_promise <- rownames(object)
            if(is.null(gN_promise)) stop(errorm05)
            
            gD_promise <- try(identifyGenes(gN_promise),
                              silent = TRUE)
            if(class(gD_promise) == "try-error" |
               length(gD_promise) == 0){
              warning(warnim02)
              out_history$gene_system <- NA
            }else{
              out_history$gene_system <- gD_promise
              message(inform03); sapply(gD_promise,
                                        function(x) message(x))
            }
            
            # if 'batchGroups' should be used
            if(!is.null(batchGroups)){
              if(length(batchGroups) != colnum)    stop(errorm06)
              if(1 == length(unique(batchGroups))) stop(errorm07)
              store$SVA <- TRUE; message(inform04)
              if(length(nsv) == 1 && is.numeric(nsv) &&
                 is.vector(nsv)){
                bn_promise <- nsv
              } else {
                bn_promise <- 1
              }
            }
            
            # if 'rmPC' should be used
            if(!is.null(rmPC)){
              col_sel <- try(seq_len(colnum)[-rmPC], silent = TRUE)
              if(class(col_sel) == "try-error" ) stop(errorm08)
              if(length(col_sel) == 0 |
                 anyNA(col_sel)) stop(errorm08)
              rm_promise <- seq_len(colnum)[-rmPC]
              store$PCA <- TRUE; message(inform05)
            }
            
            # interpretation of the 'outlier' argument
            if(!is.null(outlier)){
              ol_promise <- try(match.arg(outlier,
                                          c("zscore", "esd")), silent = TRUE)
              if(class(ol_promise) == "try-error") stop(errorm09)  
              store$outlier <- TRUE; message(eval(inform06))
            }
            
            ## SECTION1: STATISTICS  
            
            # test samples using ANOVA
            av_promise <- suppressMessages(reshape2::melt(data.frame(object)))
            colnames(av_promise) <- c("group", "y")
            anova_test <- anova( lm(y~group, data = av_promise ))
            f_sample <- anova_test$`F value`[1]; f_df <- anova_test$Df
            f_critical <- df(0.95, df1 = f_df[1] , df2 = f_df[2])
            anova_passed <- (f_sample <= f_critical)
            out_history$F_critical <- round(f_critical, 2)
            out_history$F_sample   <- round(f_sample, 2)
            out_history$F_anova    <- anova_passed
            
            # test genes using esd
            out_genes <- apply(object, 1, detectesd,
                               alpha = 0.05, rmax = 4)
            outlier_det <- (100 * sum(out_genes,
                                      na.rm = TRUE))/nrow(object)
            out_history$outlier_detected <- round(outlier_det, 1)
            
            # does the input fail for both tests
            stats_fail <- all((outlier_det > 10) && !anova_passed)
            
            # give a warning for no correction and failed tests
            if(!exists("SVA", store) &
               !exists("PCA", store) &
               !exists("outlier", store)){
              out_sobject <- object; sobject <- out_sobject
              stats_applied <- "none"
              if(stats_fail) warning(warnim03)  
            } else {
              stats_applied <- paste(ls(store), collapse = ",")
            }
            
            if(exists("SVA", store) &
               exists("PCA", store)) warning(warnim04) 
            
            if(exists("outlier", store)){
              if(ol_promise == "esd"){
                sobject <- t(apply(object, 1, correctESD,
                                   alpha = 0.05, rmax = 4))
              }
              if(ol_promise == "zscore"){
                sobject <- t(apply(object, 1, modZscore))
              }
              out_sobject <- sobject  
            } else {
              sobject <- object; out_sobject <- object
            }
            
            if(exists("PCA", store)){
              pca_object <- prcomp(sobject, center = FALSE,
                                   scale. = FALSE)  
              out_sobject <- pca_object$x[,rm_promise] %*% t(
                pca_object$rotation[,rm_promise])  
              sobject <- out_sobject
            }
            
            if(exists("SVA", store)){
              #suppressPackageStartupMessages(require(sva))    
              exbatch <- as.factor(batchGroups)
              mod1 <- model.matrix(~exbatch)
              mod0 <- cbind(mod1[,1])
              svse <- svaseq(sobject, mod1, mod0,
                             n.sv = bn_promise)$sv
              out_sobject <- svaSolv(sobject, mod1, svse)
              sobject <- out_sobject
              # detach(pos = 2L, force = TRUE) 
              #  detach(pos = 2L, force = TRUE) 
              # detach(pos = 2L, force = TRUE) 
            }
            
            ## pairwise for correlation
            if(anyNA(sobject)){
              cor_use <- "pairwise"
            } else {
              cor_use <- "everything"
            }
            
            ## SECTION2: MATRIX SEPARATION AND CORRELATION 
            # index for coding genes
            if(is.vector(pc_promise) && is.logical(pc_promise)){
              store$pc_index <- pc_promise
              out_assignment <- gN_promise[store$pc_index]
            }
            if(is.vector(pc_promise) && is.character(pc_promise)){
              store$pc_index <- is.element(pc_promise,
                                           c('protein_coding',
                                             'coding',
                                             'protein',
                                             'protein-coding',
                                             'protein coding'))
              out_assignment <- gN_promise[store$pc_index]
            }
            if(!exists("pc_index", store)) stop(errorm10)
            if(length(which(store$pc_index)) < 5) stop(errorm11)
            
            pc_matrix <- sobject[store$pc_index, ]
            nc_matrix <- sobject[!store$pc_index, ]
            
            # to calculate the correlation
            message(eval(inform07)); message(inform08)
            out_cormatrix <- try(callCor(corMethod, userFun, cor_use)(
              pc_matrix, nc_matrix), silent = TRUE)
            if(class(out_cormatrix) == "try-error") stop(errorm12)
            rownames(out_cormatrix) <- rownames(pc_matrix)
            colnames(out_cormatrix) <- rownames(nc_matrix)
            if(!all(is.finite(out_cormatrix))){
              warning(warnim05)
              out_cormatrix[is.infinite(out_cormatrix) |
                              is.na(out_cormatrix) ] <- 0
            }
            
            out_ltlmatrix <- try(callCor(corMethod, userFun, cor_use)(
              nc_matrix, nc_matrix), silent = TRUE)
            if(class(out_ltlmatrix) == "try-error") stop(errorm13)
            rownames(out_ltlmatrix) <- rownames(nc_matrix)
            colnames(out_ltlmatrix) <- rownames(nc_matrix)
            if(!all(is.finite(out_ltlmatrix))){
              out_ltlmatrix[is.infinite(out_ltlmatrix) |
                              is.na(out_ltlmatrix) ] <- 0
            }
            
            #out_history$outlier_detected <- 
            #c("corMethod", "cor_use", "cor_max", "F_critical", "F_sample",  "",  "outlier_detected")
            
            ## SECTION2: PREPARATION OF OUTPUT
            out_history$cor_max    <- round(max(out_cormatrix,
                                                na.rm = TRUE), 2)
            out_history$corMethod  <- corMethod
            out_history$cor_use    <- cor_use
            out_history$pc_matrix  <- pc_matrix
            out_history$nc_matrix  <- nc_matrix
            out_history$stats_use  <- stats_applied
            
            #out_linc             <- LINCmatrix()
            out_linc             <- new("LINCmatrix")
            results(out_linc)     <- list(statscorr = out_sobject) 
            assignment(out_linc)  <- out_assignment
            correlation(out_linc) <- list(cormatrix = out_cormatrix,
                                         lnctolnc  = out_ltlmatrix)
            express(out_linc)  <- out_object
            history(out_linc)     <- out_history
            out_linCenvir        <- NULL
            out_linCenvir        <- new.env(parent = emptyenv())
            out_linCenvir$linc   <- out_linc
            #lockEnvironment(out_linCenvir, bindings = TRUE)
            linCenvir(out_linc)   <- out_linCenvir
            
            return(out_linc)
          }) # method end

setMethod(f = "linc",
          signature = c("data.frame"),   
          definition  = function(object,
                                 codingGenes = NULL,
                                 corMethod   = "spearman",
                                 batchGroups = NULL,
                                 nsv         = 1,
                                 rmPC        = NULL,
                                 outlier     = NULL,
                                 userFun     = NULL,
                                 verbose     = TRUE){
            ## method for a data.frame as input
            ## PREDEFINITIONs
            object <- as.matrix(object)
            linc(object,
                 codingGenes,
                 corMethod,
                 batchGroups,
                 rmPC,
                 outlier,
                 userFun,
                 verbose)
          }) # method end

setMethod(f = "linc",
          signature = c("ExpressionSet"),   
          definition  = function(object,
                                 codingGenes = NULL,
                                 corMethod   = "spearman",
                                 batchGroups = NULL,
                                 nsv         = 1,
                                 rmPC        = NULL,
                                 outlier     = NULL,
                                 userFun     = NULL,
                                 verbose     = TRUE){
            ## method for an ExpressionSet as input
            ## PREDEFINITIONs
            #require(Biobase)          
            object <- Biobase::exprs(object)
            linc(object,
                 codingGenes,
                 corMethod,
                 batchGroups,
                 rmPC,
                 outlier,
                 userFun,
                 verbose)
          }) # method end

setMethod(f = "linc",
          signature = c("LINCmatrix", "missing"),   
          definition  = function(object,
                                 codingGenes = NULL,
                                 corMethod   = "spearman",
                                 batchGroups = NULL,
                                 nsv         = 1,
                                 rmPC        = NULL,
                                 outlier     = NULL,
                                 userFun     = NULL,
                                 verbose     = TRUE){
            ## method for a LINCmatrix as input
            ## PREDEFINITIONs
            
            linc(object = linCenvir(object)$expression,
                 codingGenes = linCenvir(object)$assignment,
                 corMethod,
                 batchGroups,
                 rmPC,
                 outlier,
                 userFun,
                 verbose)
          }) # method end


## getter methods
setMethod("results", "LINCcluster", function(x) x@results)
setMethod("assignment", "LINCcluster", function(x) x@assignment)
setMethod("correlation", "LINCcluster", function(x) x@correlation)
setMethod("express", "LINCcluster", function(x) x@expression)
setMethod("history", "LINCcluster", function(x) x@history)
setMethod("linCenvir", "LINCcluster", function(x) x@linCenvir)

## setter methods
setReplaceMethod("results", "LINCcluster",
                 function(x, value) {x@results <- value; validObject(x); x})

setReplaceMethod("assignment", "LINCcluster",
                 function(x, value) {x@assignment <- value; validObject(x); x})

setReplaceMethod("correlation", "LINCcluster",
                 function(x, value) {x@correlation <- value; validObject(x); x})

setReplaceMethod("express", "LINCcluster",
                 function(x, value) {x@expression <- value; validObject(x); x})

setReplaceMethod("history", "LINCcluster",
                 function(x, value) {x@history <- value; validObject(x); x})

setReplaceMethod("linCenvir", "LINCcluster",
                 function(x, value) {x@linCenvir <- value; validObject(x); x})

## create a 'LINCcluster' instance
setGeneric(name = "clusterlinc",
           def           = function(linc,
                                    distMethod    = "dicedist",
                                    clustMethod   = "average",
                                    pvalCutOff    = 0.05,
                                    padjust = "none",
                                    coExprCut     = NULL,
                                    cddCutOff     = 0.05,
                                    verbose       = TRUE){
             standardGeneric("clusterlinc")
           })
setMethod(f     = "clusterlinc",
          signature     = c("LINCmatrix"),  
          def           = function(linc,
                                   distMethod    = "dicedist",
                                   clustMethod   = "average",
                                   pvalCutOff    = 0.05,
                                   padjust = "none",
                                   coExprCut     = NULL,
                                   cddCutOff     = 0.05,
                                   verbose       = TRUE){
            ## method for a LINCmatrix
            ## PREDEFINITIONs
            
            validObject(linc)
            out_history <- history(linc)
            out_history$type <- "cluster"
            
            # errors, warnings and messages  
            errorm00 <- "LINCmatrix is empty, input corrupted"
            errorm01 <- paste("'distMethod' needs to be ",
                              "'correlation', pvalue' or 'dicedist'")  
            errorm02 <- paste("'ward.D', 'ward.D2', 'single',",
                              "'complete', 'average', 'mcquitty',",
                              "'median', 'centroid'")
            errorm03 <- "incorrect 'pvalCutOff' supplied"
            errorm04 <- "incorrect 'coExprCut' supplied"
            errorm05 <- "incorrect 'cddCutOff' supplied"
            errorm06 <- paste("clustering failed, use 'none'",
                              "in 'padjust' or 'dicedist",
                              "in 'distMethod'")
            warnim00 <- "call of hidden arguments"
            warnim01 <- paste("cor. test matrix contains infinite",
                              "or missing values; converted to 0")
            warnim02 <- paste("'corMethod' changed to 'spearman';",
                              "use 'singlelinc' to apply a user-defined",
                              "correlation test function")
            warnim03 <- paste("unable to convert 'hclust' object",
                              "output will be corrupted")
            inform01 <- paste("clusterlinc: computation for the",
                              "correlation test started")
            inform02 <- quote(paste("clusterlinc: distance matrix",
                                    "called with the method", dM_promise))
            inform03 <- paste("clusterlinc: co-expressed",
                              "genes selected based on 'coExprCut'") 
            inform04 <- paste("clusterlinc: co-expressed",
                              "genes selected based on 'pvalCutOff'") 
            
            ## SECTION0: INPUT CONTROL      
            
            if(!exists("linc", linCenvir(linc))) stop(errorm00)
            
            # interpretation of'verbose'
            if(class(verbose) != "logical" ){
              verbose <- TRUE
            } else {
              if(!any(verbose)) verbose <- FALSE
              if(any(verbose))  verbose <- TRUE
            }
            if(!verbose) message <- function(x) x
            
            # interpretation of the 'distMethod' argument
            dM_promise <- try(match.arg(distMethod,
                                        c("correlation", "pvalue", "dicedist")),
                              silent = TRUE)
            if(class(dM_promise) == "try-error") stop(errorm01)
            out_history$distMethod <- dM_promise
            
            # interpretation of the 'clustMethod' argument
            cM_promise <- try(match.arg(clustMethod, c("ward.D",
                                                       "ward.D2", "single", "complete", "average",
                                                       "mcquitty", "median", "centroid")),
                              silent = TRUE)
            if(class(cM_promise) == "try-error") stop(errorm02)
            out_history$clustMethod <- cM_promise
            
            # interpretation of the 'pvalCutOff' argument
            co_promise <- pvalCutOff
            if(length(co_promise) != 1) stop(errorm03)
            if(!is.numeric(co_promise)) stop(errorm03)
            if(co_promise >= 1 | co_promise < 0) stop(errorm03)
            out_history$pvalCutOff <- co_promise
            
            #interpretation of 'coExprCut'
            if(!is.null(coExprCut)){
              ct_promise <- try(as.integer(coExprCut), silent = TRUE)
              if(class(ct_promise) == "try-error") stop(errorm04)
              if(length(ct_promise) != 1) stop(errorm04)
              if(!is.element(ct_promise, seq_len(nrow(
                linCenvir(linc)$correlation(linc)[[1]])))) stop(errorm04)
              out_history$pvalCutOff <- NULL
            }
            
            # interpretation of 'cddCutOff'
            if(length(cddCutOff) != 1 | !is.numeric(cddCutOff) |
               !is.vector(cddCutOff)) stop(errorm05)
            if(cddCutOff > 1 | cddCutOff < 0) stop(errorm05)
            
            ## SECTION1: DISTANCE MATRIX AND CLUSTERING
            # do the correlation test
            message(inform01)
            pc_matrix <- history(linc)$pc_matrix
            nc_matrix <- history(linc)$nc_matrix
            cor_use   <- history(linc)$cor_use
            corMethod <- history(linc)$corMethod
            if(corMethod == "user"){
              corMethod <- "spearman"; warning(warnim02)
            }
            
            cortest_matrix <- suppressWarnings(doCorTest(
              corMethod, cor_use)(
                pc_matrix, nc_matrix))
            m_adjust <- p.adjust(cortest_matrix, method =
                                   padjust)
            cortest_matrix <- matrix(m_adjust, ncol = ncol(
              cortest_matrix))
            colnames(cortest_matrix) <- rownames(nc_matrix)
            rownames(cortest_matrix) <- rownames(pc_matrix)
            if(!all(is.finite(cortest_matrix))){
              warning(warnim01)
              cortest_matrix[is.infinite(cortest_matrix) |
                               is.na(cortest_matrix) ] <- 1
            }
            cortest_matrix[cortest_matrix < 1e-26] <- 1e-26
            
            out_history$pval_min <- min(cortest_matrix, na.rm = TRUE) 
            if(min(cortest_matrix, na.rm = TRUE) < 1e-26){
              out_history$pval_min <- 1e-26
            }
            
            # compute a distance matrix
            message(eval(inform02))
            
            if(dM_promise == "correlation"){
              cormat <- as.dist(correlation(linCenvir(linc)$linc)[[2]])
              distmat <- (1 - cormat)
            }
            if(dM_promise == "dicedist"){
              msize <- nrow(nc_matrix) 
              for_cdd_test <- -log10(cortest_matrix)
              blanc_matrix <- matrix(rep(-1, msize^2), ncol = msize)
              cdd_result <- docdd(for_cdd_test, blanc_matrix,
                                  (-log10(cddCutOff)))
              distmat <- as.dist(cdd_result)
            }
            if(dM_promise == "pvalue"){
              cortest_lnctolnc <- suppressWarnings(doCorTest(
                corMethod = corMethod, cor_use)(
                  nc_matrix, nc_matrix))
              
              l_adjust <- p.adjust(cortest_lnctolnc, method =
                                     padjust)
              cortest_lnctolnc  <- matrix(l_adjust, ncol = ncol(
                cortest_lnctolnc ))
              cormat <- as.dist(cortest_lnctolnc)
              distmat <- (1 - cormat)
            }
            
            # perform clustering
            out_clust <- try(hclust(distmat, method = cM_promise), silent = TRUE)
            if(class(out_clust) == "try-error") stop(errorm06)
            
            #suppressPackageStartupMessages(require("ape")) 
            if(is.function(as.phylo)){
              linc_phylo <- as.phylo(out_clust)
              linc_phylo$tip.label <- rownames(nc_matrix)
              out_clust <- linc_phylo
            } else {
              warning(warnim03)
            }
            
            # select neighbour genes
            if(!is.null(coExprCut)){
              neighbours <- deriveCluster2(cortest_matrix,
                                            n = ct_promise)
              message(inform03)
            } else {
              neighbours <- deriveCluster(cortest_matrix,
                                           alpha = co_promise) 
              message(inform04)
            }
            out_clust$neighbours <- neighbours
            
            ## SECTION4: PREPARATION OF OUTPUT
            out_linc              <- new("LINCcluster")
            results(out_linc)      <- list(cluster = out_clust) 
            assignment(out_linc)   <- assignment(linc)
            correlation(out_linc)  <- list(correlation(linc), 
                                          cortest = cortest_matrix) 
            express(out_linc)   <- express(linc)
            history(out_linc)      <- out_history
            out_linCenvir         <- new.env(parent = emptyenv())
            out_linCenvir$linc    <- linc
            out_linCenvir$cluster <- out_linc
            #lockEnvironment(out_linCenvir, bindings = TRUE)
            linCenvir(out_linc)    <- out_linCenvir
            return(out_linc)
          }) # method end

setMethod(f    = "clusterlinc",
          signature     = c("LINCcluster"),  
          def           = function(linc,
                                   distMethod    = "dicedist",
                                   clustMethod   = "average",
                                   pvalCutOff    = 0.05,
                                   coExprCut     = NULL,
                                   cddCutOff     = 0.05,
                                   verbose       = TRUE){
            ## method for a LINCcluster
            ## PREDEFINITIONs
            clusterlinc(linc          = linCenvir(linc)$linc,
                        distMethod    = distMethod,
                        clustMethod   = clustMethod,
                        pvalCutOff    = pvalCutOff,
                        coExprCut     = coExprCut,
                        cddCutOff     = cddCutOff,
                        verbose       = verbose)
          }) # method end


## getter methods
setMethod("results", "LINCbio", function(x) x@results)
setMethod("assignment", "LINCbio", function(x) x@assignment)
setMethod("correlation", "LINCbio", function(x) x@correlation)
setMethod("express", "LINCbio", function(x) x@expression)
setMethod("history", "LINCbio", function(x) x@history)
setMethod("linCenvir", "LINCbio", function(x) x@linCenvir)

## setter methods
setReplaceMethod("results", "LINCbio",
                 function(x, value) {x@results <- value; validObject(x); x})

setReplaceMethod("assignment", "LINCbio",
                 function(x, value) {x@assignment <- value; validObject(x); x})

setReplaceMethod("correlation", "LINCbio",
                 function(x, value) {x@correlation <- value; validObject(x); x})

setReplaceMethod("express", "LINCbio",
                 function(x, value) {x@expression <- value; validObject(x); x})

setReplaceMethod("history", "LINCbio",
                 function(x, value) {x@history <- value; validObject(x); x})

setReplaceMethod("linCenvir", "LINCbio",
                 function(x, value) {x@linCenvir <- value; validObject(x); x})

## create a LINCbio instance
setGeneric(name = "getbio",
           def            = function(cluster,
                                     translate      = "mygene",
                                     annotateFrom   = 'enrichGO',
                                     verbose        = TRUE,
                                     ...){
             standardGeneric("getbio")
           })
setMethod(f     = "getbio",
          signature       = c("LINCcluster"),  
          def            = function(cluster,
                                    translate      = "mygene",
                                    annotateFrom   = 'enrichGO',
                                    verbose        = TRUE,
                                    ...){                                
            ## method for a LINCcluster
            ## PREDEFINITIONs
            
            # environments and object
            validObject(cluster)
            store <- new.env(parent = emptyenv())
            out_history <- history(cluster)
            
            # errors, warnings and messages  
            errorm00 <- "LINCmatrix is empty, input corrupted"
            errorm01 <- paste("'translate' needs to be 'mygene",
                              " or 'none'")  
            errorm02 <- "interpretation of 'annotateFrom' failed"
            errorm03 <- "incorrect 'pvalCutOff' supplied"
            warnim00 <- "call of hidden arguments"
            warnim01 <- paste("unintended call to 'mygene' in a case",
                              "where gene ids do not belong to the Ensembl",
                              "gene id system")
            warnim02 <- "translation by 'mygene' time-consuming "
            inform01 <- paste("getbio: Package 'clusterProfiler' was",
                              "called for gene annotation")
            inform02 <- paste("getbio: Package 'DOSE' was",
                              "called for gene annotation")
            inform03 <- paste("getbio: Package 'ReactomePA' was",
                              "called for gene annotation")
            inform04 <- paste("getbio: Package 'mygene' was",
                              "called for gene id translation")
            
            ## SECTION0: INPUT CONTROL      
            
            # interpretation of'verbose'
            if(class(verbose) != "logical" ){
              verbose <- TRUE
            } else {
              if(!any(verbose)) verbose <- FALSE
              if(any(verbose))  verbose <- TRUE
            }
            if(!verbose) message <- function(x) x
            
            # interpretation of 'translate'
            tl_promise <- try(match.arg(translate,
                                        c("mygene", "xxx", "none")),
                              silent = TRUE)
            if(class(tl_promise) == "try-error") stop(errorm01)
            out_history$translate <- tl_promise
            
            # interpretation of annotateFrom
            cP_promise <- try(match.arg(annotateFrom,
                                        c("enrichGO",
                                          "enrichKEGG",
                                          "enrichPathway",
                                          "enrichDO")),
                              silent = TRUE)
            if(class(cP_promise) == "try-error") stop(errorm02)
            if(any(is.element(cP_promise, c("enrichGO",
                                            "enrichKEGG",
                                            "enrichDO",
                                            "enrichPathway")))) {
              #suppressPackageStartupMessages(require(clusterProfiler))  
              #suppressPackageStartupMessages(require(org.Hs.eg.db))
              message(inform01); store$cP_promise <- cP_promise
              if(cP_promise == "enrichDO"){
                #suppressPackageStartupMessages(require(DOSE))
                message(inform02)
              }
              if(cP_promise == "enrichPathway"){
                #suppressPackageStartupMessages(require(ReactomePA))  
                message(inform03)
              }
            }
            
            ## !! FROM THE ENVIRONMENT
            qy_promise <- results(cluster)[[1]]$neighbours
            
            ## SECTION1: GENE TRANSLATION
            gn_promise <- identifyGenes(unlist(qy_promise))
            
            if(tl_promise == "mygene"){
              if(!any(is.element(gn_promise, "ensemblgene"))){
                warning(warnim01)
              }
              # suppressPackageStartupMessages(require(mygene))  
              message(inform04); warning(warnim02)
              transThis <- function(x){
                query <- mygene::getGenes(x, fields = "entrezgene")
                entrez_query <- query@listData$entrezgene
                return(entrez_query)
              }
            }
            if(tl_promise == "none"){
              transThis <- function(x) x
            }
            
            ## SECTION2: CALL TO GENE ANNOTATION
            if(exists("cP_promise", store)){
              annotateFun <- match.fun(cP_promise)
              
              OrgDb <- 'org.Hs.eg.db'
              if(cP_promise == "enrichPathway"){ OrgDb <- "human"}
              
              out_bio_terms <- lapply(qy_promise,
                                      function(x){ entrez_query <- transThis(x)
                                      cP_result   <- annotateFun(gene = entrez_query, OrgDb, ...)
                                      if(class(cP_result) == "enrichResult"){
                                        qvalues     <- slot(cP_result, "result")$qvalue 
                                        bio_terms   <- slot(cP_result, "result")$Description
                                        query_entry <- list(qvalues, bio_terms)
                                      } else {
                                        query_entry <- list(NA, NA)
                                      }
                                      return(query_entry)
                                      })  
            } # clusterProfiler 
            
            ## SECTION3: PREPARATION OF OUTPUT
            #out_linc             <- LINCbio()
            out_linc             <- new("LINCbio")
            results(out_linc)     <- list(bio = out_bio_terms)  
            assignment(out_linc)  <- assignment(cluster)
            correlation(out_linc) <- correlation(cluster)
            expression(out_linc)  <- express(cluster)
            history(out_linc)     <- out_history
            out_linCenvir        <- NULL
            out_linCenvir        <- linCenvir(cluster)
            out_linCenvir$bio    <- out_linc
            #lockEnvironment(out_linCenvir, bindings = TRUE)
            linCenvir(out_linc)   <- out_linCenvir
            return(out_linc)
          }) # method end


## PLOTTING METHOD
setGeneric(name = "plotlinc",
           def = function(
             input,
             showCor,
             returnDat = FALSE){
             standardGeneric("plotlinc") 
           })
setMethod(f   = "plotlinc",
          signature = c("LINCmatrix",
                        "missing"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            
            ## method for a LINCbio object            
            ## PREDEFINITIONs
            # on.exit(options(stringsAsFactors = TRUE))         
            
            validObject(input)       
            em_promise  <- results(linCenvir(input)$linc)[[1]]
            ac_promise  <- correlation(linCenvir(input)$linc)[[1]]
            hs_promise <-  history(linCenvir(input)$linc)
            ep_promise  <- express(linCenvir(input)$linc)
            # suppressPackageStartupMessages(require(reshape))
            
            # create a box plot
            df_boxplot  <- suppressMessages(reshape2::melt(em_promise))
            names(df_boxplot) <- c(NA, "SAMPLES", "VALUE")
            gg_box <- ggplot(df_boxplot,
                             aes(SAMPLES, VALUE)) + geom_boxplot(outlier.color =
                                                                   "firebrick", colour = "dodgerblue3") +  theme(
                                                                     panel.border = element_rect(color = "grey", fill = NA),   
                                                                     panel.background = element_blank()) +
              ggtitle("BOXPLOT OF EXPRESSION VALUES") +
              theme(plot.title = element_text(face = "bold",
                                              color = "steelblue4"))
            
            # create a frequency plot
            gg_freq <- ggplot(df_boxplot, aes(VALUE)) +
              geom_histogram(bins = 30, fill = "gray95" ) +
              geom_freqpoly(colour = "firebrick", linetype = "dashed", size = 0.7) +
              theme(
                panel.border = element_rect(color = "grey", fill = NA),   
                panel.background = element_blank()) +
              ggtitle("FREQUENCY OF EXPRESSION VALUES") +
              theme(plot.title = element_text(face = "bold",
                                              color = "gray47"))
            
            
            # create a histogram of correlation values
            df_cor <- data.frame(CORRELATION = as.vector(ac_promise))
            gg_cor <- ggplot(df_cor, aes(CORRELATION)) + geom_histogram(binwidth = 0.02,  #bins = 100,
                                                                        size = 1, fill = c(rep("grey", 95), rep("dodgerblue3", 6))) +  theme(
                                                                          panel.border = element_rect(color = "grey", fill = NA),   
                                                                          panel.background = element_blank()) +
              xlim(-1,1) + #scale_x_continuous(breaks = 0.01 ) +
              geom_vline(xintercept = 0, colour = "firebrick", linetype = 'dashed') +
              geom_vline(xintercept = 0.9, colour = "dodgerblue3") + 
              ggtitle("HISTOGRAM OF CORRELATION VALUES") +
              theme(plot.title = element_text(face = "bold",
                                              color = "violetred4"))
            
            # plot PCA analysis
            pca_object <- prcomp(ep_promise, center = FALSE,
                                 scale. = FALSE)  
            df_pca <- data.frame(PC = seq_len(length(pca_object$sdev)),
                                 VARIANCE = (pca_object$sdev/sum(pca_object$sdev) * 100))
            gg_pca <- ggplot(df_pca, aes(PC, VARIANCE)) + geom_point(
              size = 4, colour = "dodgerblue3") +  theme(
                panel.border = element_rect(color = "grey", fill = NA),   
                panel.background = element_blank()) +
              ylim(0,100) + scale_x_continuous(breaks=seq(1, 
                                                          length(pca_object$sdev),1)) +
              ggtitle("PRINCIPLE COMPONENT ANALYSIS") +
              theme(plot.title = element_text(face = "bold",
                                              color = "grey25"))   
            
            # get and plot the statistical parameters
            get_this <- c("corMethod", "cor_use", "cor_max",
                          "F_critical", "F_sample",  "F_anova", 
                          "outlier_detected", "stats_use")
            par_description <- c("Method for correlation:  ",
                                 "Usage of observations:   ",
                                 "Maximum cor. observed:   ",
                                 "Critical F-value:        ",
                                 "F-value of sample:       ",
                                 "ANOVA passed:            ",
                                 "Genes with outliers [%]: ",
                                 "Correction applied:      ")
            stat_par <- mget(get_this, envir = hs_promise) 
            par_described <- mapply(paste, par_description, stat_par)
            df_stat <- data.frame(value = c(" ", par_described, " "),
                                  y = -c(1:10), x = rep(0,10), group =
                                    c(rep("cor", 4), rep("samples", 4), rep("expr", 2)))
            
            pty_pl <- (ggplot(df_stat, aes(x,y)) +
                         geom_point(color = "white") + xlim(0, 1) +
                         theme(axis.line = element_line(colour =
                                                          "white"), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),   
                               panel.background = element_blank(),
                               axis.title.y = element_text(color ="white"),
                               axis.title.x = element_text(color ="white"),
                               axis.text.x = element_text(color = "white"),
                               axis.text.y = element_text(color = "white")))
            
            linc_stats_plot <- pty_pl + geom_text(aes(label = 
                                                        df_stat$value, colour = df_stat$group) ,hjust = 0, vjust = 0,
                                                  size = 5, fontface = "bold") + ggtitle("STATISTICAL PARAMETERS") +
              scale_colour_manual(values = c(
                "violetred4", "gray47", "steelblue4"), guide = "none") + 
              theme(plot.title = element_text(face = "bold"))
            
            #suppressPackageStartupMessages(require(grid))  
            #  suppressPackageStartupMessages(require(png))
            stats_img <- readPNG(system.file("extdata", "statslinc_img.png",
                                             package ="LINC"))
            #stats_img <- readPNG("statslinc_img.png")
            stats_plot <- rasterGrob(stats_img, interpolate = TRUE)
            
            #  suppressPackageStartupMessages(require(gridExtra)) 
            
            customid <- ""
            if(exists("customID", envir = history(input))){
              customid <- history(input)$customID
            }
            
            left_side <- suppressMessages(suppressWarnings(arrangeGrob(
              gg_cor, gg_box , gg_pca, ncol = 1)))
            
            right_side <- arrangeGrob(stats_plot, linc_stats_plot, ncol = 1, bottom = customid)
            
            grid.arrange(left_side, right_side,  nrow = 1 )
            
          })


## plot a cluster without bio_terms
setMethod(f   = "plotlinc",
          signature = c("LINCcluster",
                        "missing"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            
            ## method for a LINCcluster object            
            ## PREDEFINITIONs
            validObject(input)         
            #  suppressWarnings(suppressPackageStartupMessages(
            # require(ggtree)))          
            
            cluster  <- results(linCenvir(input)$cluster)[[1]]
            
            ## SECTION0: INPUT CONTROL  
            # check for a cluster
            
            #+ to be added          
            
            # plot the cluster
            tree <- ggtree(cluster, colour = "firebrick") +
              coord_flip() + scale_x_reverse()
            tree <- tree + geom_tiplab(size = 3.5, angle = 90,
                                       colour = "firebrick", hjust = 1)
            
            # derive neighbour genes for different CutOffs
            cortest_matrix  <- correlation(linCenvir(input)$cluster)[[2]]
            pval_default <- c(1e-25, 1e-10, 1e-5, 1e-4, 1e-3, 1e-2)
            nrow <- length(pval_default)
            msize <- dim(cortest_matrix); mextn <- nrow  * msize[2]
            fill_matrix <- matrix(rep(0, mextn), ncol = msize[2],
                                  nrow = nrow )
            
            for(i in seq_len(nrow) ){
              ng_derived <- deriveCluster(cortest_matrix,
                                           alpha = pval_default[i])
              fill_matrix[i,] <- unlist(lapply(ng_derived, length))
            }
            
            
            # convert this matrix for plotting
            plot_df <- as.data.frame(t(fill_matrix))
            rownames(plot_df) <- colnames(cortest_matrix )
            colnames(plot_df) <-  as.character(pval_default)
            
            clust_heat <- gheatmap(tree, plot_df, offset = 0.1,
                                   width = 0.5, colnames = TRUE,
                                   low = "white", high = "violetred4",
                                   colnames_position = "top")

            hs_promise <- history(linCencir(input)$cluster)

            # get and plot the statistical parameters
            get_this <- c("distMethod", "clustMethod", 
                          "corMethod",  "cor_use", "pvalCutOff",
                          "pval_min" , "stats_use", "gene_system")
            par_description <- c("Distance measure:       ",
                                 "Clustering algorithm:   ",
                                 "Method for correlation: ",
                                 "Usage of observations:  ",
                                 "p-value cut-off:        ",
                                 "Best p-value observed:  ",
                                 "Statistical correction: ",
                                 "Gene system identified: ")
            
            stat_par <- mget(get_this, envir = hs_promise) 
            par_described <- mapply(paste, par_description, stat_par)
            df_stat <- data.frame(value = c(" ", par_described, " "),
                                  y = -c(1:10), x = rep(0,10))
            
            pty_pl <- (ggplot(df_stat, aes(x,y)) +
                         geom_point(color = "white") + xlim(0, 1) +
                         theme(axis.line = element_line(colour =
                                                          "white"), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),   
                               panel.background = element_blank(),
                               axis.title.y = element_text(color ="white"),
                               axis.title.x = element_text(color ="white"),
                               axis.text.x = element_text(color = "white"),
                               axis.text.y = element_text(color = "white")))
            
            clust_para_plot <- pty_pl + geom_text(aes(label = 
                                                        df_stat$value) ,hjust = 0, vjust = 0,
                                                  size = 5, fontface = "bold",
                                                  colour = "gray47") +
              ggtitle(paste("PARAMETERS OF CLUSTER",
                            "AND COEXPRESSED GENES")) +
              theme(plot.title = element_text(face = "bold"))
            
            #  suppressPackageStartupMessages(require(grid))  
            #  suppressPackageStartupMessages(require(png))
            stclust_img <- readPNG(system.file("extdata", "stclust_img.png",
                                               package ="LINC"))
            # stclust_img <- readPNG("stclust_img.png")
            stclust_plot <- rasterGrob(stclust_img, interpolate = TRUE)
            
            # suppressPackageStartupMessages(require(gridExtra))     
            # left_side <- suppressMessages(suppressWarnings(arrangeGrob(
            #  gg_cor, gg_box , gg_freq, ncol = 1)))
            customid <- ""
            if(exists("customID", envir = history(input))){
              customid <- history(input)$customID
            } 
            
            right_side <- arrangeGrob(stclust_plot, clust_para_plot, ncol = 1, bottom = customid)
            
            grid.arrange(clust_heat, right_side,  nrow = 1 )
            
          }) # end of method  

setMethod(f   = "plotlinc",
          signature = c("LINCsingle",
                        "missing"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            
            ## method for a LINCsingle object            
            ## PREDEFINITIONs
            # on.exit(options(stringsAsFactors = TRUE)) 
            validObject(input)  
            #require(RColorBrewer)
            
            linCenvir(input)$single
            
            query  <- results(linCenvir(input)$single)$query
            bio    <- results(linCenvir(input)$single)$bio
            corl   <- results(linCenvir(input)$single)$cor            
            pval   <- results(linCenvir(input)$single)$pval 
            pval10 <- -log10(unlist(pval))
            if(all(is.na(pval))) pval[1] <- 0 # ggplot rescue
            
            # limit the terms to plot
            size <- length(bio[[2]])
            priority <- paste("[", 1:9, "]", sep = '')
            if(size >= 9){
              bio[[2]] <- bio[[2]][1:9]
              bio[[1]] <- bio[[1]][1:9]
            } else {
              bio[[2]][(size + 1):9] <- "NA"
              bio[[1]][(size + 1):9] <- 1
            }
            
            # make the data.frame for plotting
            priority <- paste("[", 1:9, "]", sep = '')
            term_assign <- mapply(function(x, y) { paste(x, y) }, x = priority, y = bio[[2]])
            annotation_df <- data.frame(priority, -log10(bio[[1]]), term_assign)
            names(annotation_df) <- c("TERMS", "pvalue", "ASSIGNMENT")
            ############################################################  
            
            # plot the annotation
            annotation_plot <- ggplot(annotation_df, aes(TERMS,
                                                         y = pvalue, fill = ASSIGNMENT)) + geom_bar( stat =
                                                                                                       "identity", width = 0.8) +
              theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank()) +
              theme(axis.text.x = element_text(size = 15,
                                               hjust = 0, vjust = 0, colour = "violetred4")) +
              scale_fill_manual(values= c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6"),
                                name="ANNOTATION",
                                breaks= term_assign,
                                labels= term_assign) +
              theme(legend.text = element_text(size = 15, colour = "violetred4"),
                    legend.title = element_text(size = 17, colour = "#023858")) +
              ggtitle(paste("QUERY:", query, ";", length(corl), "CO-EXPRESSED GENES" )) + theme(plot.title = element_text(size = 17, face = "bold", vjust = -3, hjust = 1)) +
              theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 16)) 
            
            #'-log10(p-value)'
            
            # plot histograms of correlation and p-values
            ############################################################    
            # create a histogram of correlation values
            df_cor <- data.frame(CORRELATION = unlist(corl))
            gg_cor <- ggplot(df_cor, aes(CORRELATION)) +
              geom_histogram(binwidth = 0.02, size = 1, fill = "grey") +
              theme(panel.border = element_rect(color = "grey",#
                                                fill = NA), panel.background = element_blank()) +
              xlim(-1,1) + geom_vline(xintercept = 0, colour =
                                        "firebrick", linetype = 'dashed') + geom_vline(
                                          xintercept = 0.75, colour = "dodgerblue3") + 
              ggtitle("CO-EXPRESSED GENES: CORRELATION VALUES") +
              theme(plot.title = element_text(face = "bold",
                                              color = "ivory4"))
            
            df_pval <- data.frame(pvalue = pval10)
            gg_pval <- ggplot(df_pval, aes(pvalue)) +
              geom_histogram(binwidth = 1, size = 1, fill = "grey") +
              theme(panel.border = element_rect(color = "grey",
                                                fill = NA), panel.background = element_blank()) +
              xlim(0,100) + geom_vline(xintercept = -log10(0.05),
                                       colour = "firebrick", linetype = 'dashed') + geom_vline(
                                         xintercept = 16, colour = "dodgerblue3") + 
              ggtitle("CO-EXPRESSED GENES: P-VALUES (COR.TEST)") +
              theme(plot.title = element_text(face = "bold",
                                              color = "lightpink4"))
            
            # suppressPackageStartupMessages(require(grid))  
            #suppressPackageStartupMessages(require(png))
            single_img <- readPNG(system.file("extdata", "singlelinc_img.png",
                                              package ="LINC"))
            #single_img <- readPNG("singlelinc_img.png")
            single_plot <- rasterGrob(single_img, interpolate = TRUE)
            
            customid <- ""
            if(exists("customID", envir = history(input))){
              customid <- history(input)$customID
            } 
            
            # suppressPackageStartupMessages(require(gridExtra))     
            right_bottom <- suppressMessages(suppressWarnings(arrangeGrob(
              gg_cor, gg_pval, ncol = 1)))
            right_side <- arrangeGrob(single_plot, right_bottom, ncol = 1, bottom = customid)
            
            grid.arrange( annotation_plot, right_side, nrow = 1 )
            
          })


setMethod(f   = "plotlinc",
          signature = c("LINCbio",
                        "missing"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            
            ## method for a LINCbio object            
            ## PREDEFINITIONs
            validObject(input)        
            # suppressWarnings(suppressPackageStartupMessages(
            # require(ggtree)))          
            
            cluster  <- results(linCenvir(input)$cluster)[[1]]
            bio_list <- results(linCenvir(input)$bio)[[1]]
            
            ## SECTION0: INPUT CONTROL  
            # check for a cluster
            
            #+ to be added          
            
            # plot the cluster
            tree <- ggtree(cluster, colour = "firebrick") +
              coord_flip() + scale_x_reverse()
            tree <- tree + geom_tiplab(size = 3.5, angle = 90,
                                       colour = "firebrick", hjust = 1)
            
            # prepare biological terms for plotting
            # preparation of a matrix
            
            # MAKE THIS ROBUST
            term_ext <- 1
            while(term_ext < 20){
              term_ext <- (term_ext + 1)
              raw_names <- names(bio_list)
              term_crude <- lapply(bio_list, function(x, term_ext){
                x[[2]][seq_len(term_ext)] }, term_ext)
              term_unique <- unique(unlist(term_crude))
              if(length(term_unique) > 20) break
            } 
            term_unique[is.na(term_unique)] <- "NA"
            m <- length(raw_names); n <- length(term_unique)
            plot_matrix <- matrix(rep(0, (m*n)), ncol = n, nrow = m )
            colnames(plot_matrix) <- term_unique
            rownames(plot_matrix) <- raw_names
            
            # now fill matrix with biological terms
            for(query in seq_len(m)){
              if(length(bio_list[[query]][[2]]) > 0 &&
                 is.character(bio_list[[query]][[2]]) &&
                 is.numeric(bio_list[[query]][[1]])){
                
                bio_query <- bio_list[[query]][[2]]
                pvalues <- (-log10(bio_list[[query]][[1]]))
                row_entry <- vapply(colnames(plot_matrix),
                                    function(x, pvalues){
                                      if(any(x == bio_query)){
                                        pvalues[x == bio_query][1]
                                      } else {
                                        return(0)   
                                      }
                                    }, 0, pvalues)
                row_entry[row_entry > 16] <- 16
                plot_matrix[query,] <- row_entry
              }
            }
            
            # melt this matrix for plotting

            # convert this matrix for plotting
            term_assign <- seq_len(length(term_unique))
            colnames(plot_matrix) <- paste(term_assign, "          ")
            #plot_matrix > 16 <- 16 
            plot_df <- as.data.frame(plot_matrix)
            
            # return data or return the plot
            if(returnDat){
              colnames(plot_df) <- term_unique
              final_dat <- list(cluster = cluster, annotation = plot_df)
              return(final_dat) 
            } else { 
              clust_heat <- gheatmap(tree, plot_df, offset = 0.1,
                                     width = 0.5, colnames = TRUE,
                                     low = "white", high = "dodgerblue3",
                                     colnames_position = "top")
              
              #grid.arrange(clust_heat, clust_heat, ncol= 0, nrow=0)
              
              
              # additional plot for term assignments
              bio_assignment <- mapply(function(x,y){ paste(x,y) },
                                       x = term_assign, y = term_unique)
              df_assign <- data.frame(bio_assignment, y = -term_assign,
                                      x = rep(0, length(term_unique)))
              
              pty_pl <- (ggplot(df_assign, aes(x,y)) +
                           geom_point(color = "white") + xlim(0, 1) +
                           theme(axis.line = element_line(colour =
                                                            "white"), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),   
                                 panel.background = element_blank()) +
                           theme(axis.title.y = element_text(color =
                                                               "white")) + theme(axis.title.x = element_text(
                                                                 color = "white")) + theme(axis.text.x = 
                                                                                             element_text(color = "white")) + theme(
                                                                                               axis.text.y = element_text(color = "white")))
              
              bio_assign_plot <- pty_pl + geom_text(aes(label = 
                                                          bio_assignment),hjust=0, vjust=0,
                                                    size = 4, colour = "grey18") 
              
              # additional plot for title and explanations  
              # suppressPackageStartupMessages(require(grid))  
              # suppressPackageStartupMessages(require(png))
              
              clust_img <- readPNG(system.file("extdata", "clust_img.png",
                                               package ="LINC"))
              #clust_img <- readPNG("clust_img.png")
              clust_plot <- rasterGrob(clust_img, interpolate = TRUE)
              
              
              # assembly of the final plot
              # suppressPackageStartupMessages(require(gridExtra))    
              customid <- ""
              if(exists("customID", envir = history(input))){
                customid <- history(input)$customID
              } 
              
              right_side <- arrangeGrob(clust_plot, bio_assign_plot,
                                        ncol = 1, bottom = customid)
              final_plot <- grid.arrange(clust_heat, right_side,
                                         nrow = 1)
              return(invisible(final_plot))  
            }
          }) # method end  


# method for showCor
setMethod(f   = "plotlinc",
          signature = c("LINCmatrix",
                        "character"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            
            validObject(input)  
            input <- (input + feature(setLevel = "LINCmatrix"))
            returnDat <- inlogical(returnDat, direct = FALSE)
            if(returnDat) warning("'returnDat' unused in this method")
            
            errorm01 <- paste("'showCor' must be a vector",
                              "supplying 2 up to 6 gene ids")            
            errorm02 <- "unable to find all gene ids in input"       
            
            # check showCor
            spg_promise <- showCor
            if(length(spg_promise) < 2 | length(spg_promise) > 6 |
               is.numeric(spg_promise))
              stop(errorm01)
            
            select_try <- try(is.element(showCor,
                                         rownames(express(input))), silent = TRUE)             
            if(class(select_try) == "try-error") stop(errorm01)      
            if(!all(select_try)) stop(errorm02)
            
            # predefine empty vectors
            for(n in 2:6){
              assign(paste("SUBJECT", n, sep = '_'),
                     rep(NA, ncol(express(input))))
            }
            
            # define input
            QUERY <- results(input)[[1]][spg_promise[1], ]
            for(n in (c(2:length(spg_promise))) - 1){
              assign(paste("SUBJECT", n, sep = '_'),
                     results(input)[[1]][spg_promise[n + 1], ])
            }
            
            
            expr_df <- data.frame(EXPRESSION = QUERY,
                                  SAMPLES = seq_len(ncol(express(input))), SUBJECT_1,
                                  SUBJECT_2,  SUBJECT_3, SUBJECT_4, SUBJECT_5,
                                  NO_SUBJECT = rep(0, ncol(express(input))))
            
            qy_gg <- ggplot(expr_df) 
            
            qy_gg_ns <- ggplot(expr_df) + geom_bar(aes(x = SAMPLES, 
                                                       y = EXPRESSION), stat = "identity", colour =
                                                     "firebrick", fill = "red", alpha = 0.1 ) +
              theme(panel.background = element_blank(),
                    panel.border = element_rect(color = "grey",
                                                fill = NA))
            
            cor1 <-  try(correlation(input)[[1]][spg_promise[2],
                                                spg_promise[1]], silent = TRUE)
            if(class(cor1) == "try-error") cor1 <- NA
            plot1 <- qy_gg + geom_point(aes(x =QUERY, y = SUBJECT_1),
                                        stat = "identity", colour = "firebrick4", alpha = 0.9) + ggtitle(paste(
                                          spg_promise[1], "vs", spg_promise[2],
                                          "(correlation =", round(cor1, 3), ")" )) +
              theme(title = element_text(face = "bold",
                                         size = 15)) 
            
            if(length(spg_promise) > 2){
              cor2 <-  try(correlation(input)[[1]][spg_promise[3],
                                                  spg_promise[1]], silent = TRUE)
              if(class(cor2) == "try-error") cor2 <- NA
              plot2 <- qy_gg + geom_point(aes(x =QUERY, y = SUBJECT_2),
                                          stat = "identity", colour = "firebrick4", alpha = 0.9) + ggtitle(paste( 
                                            spg_promise[1], "vs", spg_promise[3],
                                            "(correlation =", round(cor2, 3), ")" )) +
                theme(title = element_text(face = "bold",
                                           size = 15))
            } else {
              plot2 <- qy_gg_ns + geom_bar(aes(x = SAMPLES, y = NO_SUBJECT),
                                           stat = "identity") + ggtitle(spg_promise[1]) +
                theme(title = element_text(face = "bold", size = 15))
            }
            
            
            if(length(spg_promise) > 3){
              cor3 <-  try(correlation(input)[[1]][spg_promise[4],
                                                  spg_promise[1]], silent = TRUE)
              if(class(cor3) == "try-error") cor3 <- NA
              plot3 <- qy_gg + geom_point(aes(x =QUERY, y = SUBJECT_3),
                                          stat = "identity", colour = "firebrick4", alpha = 0.9) + ggtitle(paste(
                                            spg_promise[1], "vs", spg_promise[4],
                                            "(correlation =", round(cor3, 3), ")" )) +
                theme(title = element_text(face = "bold",
                                           size = 15))
            } else {
              plot3 <- qy_gg_ns + geom_bar(aes(x = SAMPLES, y = NO_SUBJECT),
                                           stat = "identity") + ggtitle(spg_promise[1]) +
                theme(title = element_text(face = "bold", size = 15))
            }
            
            
            if(length(spg_promise) > 4){
              cor4 <-  try(correlation(input)[[1]][spg_promise[5],
                                                  spg_promise[1]], silent = TRUE)
              if(class(cor4) == "try-error") cor4 <- NA
              plot4 <- qy_gg + geom_point(aes(x =QUERY, y = SUBJECT_4),
                                          stat = "identity", colour = "firebrick4", alpha = 0.9) + ggtitle(paste(
                                            spg_promise[1], "vs", spg_promise[5],
                                            "(correlation =", round(cor4, 3), ")" )) +
                theme(title = element_text(face = "bold",
                                           size = 15))
            } else {
              plot4 <- qy_gg_ns + geom_bar(aes(x = SAMPLES, y = NO_SUBJECT),
                                           stat = "identity") + ggtitle(spg_promise[1]) +
                theme(title = element_text(face = "bold", size = 15))
            }
            
            
            if(length(spg_promise) > 5){
              cor5 <-  try(correlation(input)[[1]][spg_promise[6],
                                                  spg_promise[1]], silent = TRUE)
              if(class(cor5) == "try-error") cor5 <- NA
              plot5 <- qy_gg + geom_point(aes(x =QUERY, y = SUBJECT_5),
                                          stat = "identity", colour = "firebrick4", alpha = 0.9) + ggtitle(paste(
                                            spg_promise[1], "vs", spg_promise[6],
                                            "(correlation =", round(cor5, 3), ")" )) +
                theme(title = element_text(face = "bold",
                                           size = 15))
            } else {
              plot5 <- qy_gg_ns + geom_bar(aes(x = SAMPLES, y = NO_SUBJECT),
                                           stat = "identity") + ggtitle(spg_promise[1]) +
                theme(title = element_text(face = "bold", size = 15))
            }
            
            # arrange these plots
            # additional plot for title and explanations  
            # suppressPackageStartupMessages(require(grid))  
            #suppressPackageStartupMessages(require(png))
            
            barplot_img <- readPNG(system.file("extdata", "barplot_img.png",
                                               package ="LINC"))
            #clust_img <- readPNG("protcl_img.png")
            bar_plot <- rasterGrob(barplot_img, interpolate = TRUE)
            
            # assembly of the final plot
            customid <- ""
            if(exists("customID", envir = history(input))){
              customid <- history(input)$customID
            } 
            bar_plot <- arrangeGrob(bar_plot)
            final_plot <- grid.arrange(plot1, bar_plot, plot2,
                                       plot3, plot4, plot5,
                                       nrow = 3, ncol = 2)
            return(invisible(final_plot))
          })

setMethod(f   = "plotlinc",
          signature = c("LINCcluster",
                        "character"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            callNextMethod()
          })
setMethod(f   = "plotlinc",
          signature = c("LINCbio",
                        "character"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            callNextMethod()
          }) 



## getter methods
setMethod("results", "LINCsingle", function(x) x@results)
setMethod("assignment", "LINCsingle", function(x) x@assignment)
setMethod("correlation", "LINCsingle", function(x) x@correlation)
setMethod("express", "LINCsingle", function(x) x@expression)
setMethod("history", "LINCsingle", function(x) x@history)
setMethod("linCenvir", "LINCsingle", function(x) x@linCenvir)

## setter methods
setReplaceMethod("results", "LINCsingle",
                 function(x, value) {x@results <- value; validObject(x); x})

setReplaceMethod("assignment", "LINCsingle",
                 function(x, value) {x@assignment <- value; validObject(x); x})

setReplaceMethod("correlation", "LINCsingle",
                 function(x, value) {x@correlation <- value; validObject(x); x})

setReplaceMethod("express", "LINCsingle",
                 function(x, value) {x@expression <- value; validObject(x); x})

setReplaceMethod("history", "LINCsingle",
                 function(x, value) {x@history <- value; validObject(x); x})

setReplaceMethod("linCenvir", "LINCsingle",
                 function(x, value) {x@linCenvir <- value; validObject(x); x})

setMethod(f   = "plotlinc",
          signature = c("LINCsingle",
                        "character"),  
          def = function(
            input,
            showCor,
            returnDat = FALSE){ 
            callNextMethod()
          })  

## create a LINCsingle instance
setGeneric(name = "singlelinc",
           def  = function(input,
                           query = NULL,
                           onlycor = FALSE,
                           testFun = cor.test,
                           alternative  = "greater",
                           threshold = 0.05,
                           underth = TRUE,
                           coExprCut = NULL,
                           handleGeneIds = TRUE,
                           annotateFrom  = 'enrichGO',
                           verbose = TRUE, ...){
             standardGeneric("singlelinc")
           })

setMethod(f = "singlelinc",
          signature = c("LINCmatrix"),   
          definition  = function(input,
                                 query = NULL,
                                 onlycor = FALSE,
                                 testFun = cor.test,
                                 alternative  = "greater",
                                 threshold = 0.05,
                                 underth = TRUE,
                                 coExprCut = NULL,
                                 handleGeneIds = TRUE,
                                 annotateFrom  = 'enrichGO',
                                 verbose = TRUE, ...){
            ## method for a LINCmatrix as input
            ## PREDEFINITIONs
            # errors, warnings and messages  
            #errorm00 <- "Input object is empty"  
            errorm01 <- "'testFun' is not a valid function"
            errorm02 <- paste("'testFun' does not work; its",
                              "output must be available by '$pvalue'")
            errorm03 <- "'coExprCut' has to be a single integer"
            errorm04 <- "'threshold' has to be a single numeric value"
            errorm05 <- "this function was called without a 'query'"
            errorm06 <- "input 'query' not found; possible queries:"
            errorm07 <- "no co-expressed genes for this 'threshold'"
            
            warnim01 <- paste("'threshold' usually between 0 and 1; ",
                              "for a user-defined 'testFun' it could be different")
            warnim02 <- paste("'testFun' was supplied and 'onlycor'",
                              "equals 'TRUE', here 'onlycor' has the higher priority")
            warnim03 <- paste("'testFun' is not 'stats::cor.test'",
                              "and does not have formal arguments",
                              "'x', 'y', 'use' and 'method'")
            warnim04 <- "input not as Entrezgenes and not translated"
            inform01 <- paste("singlelinc: no test conducted, genes",
                              "selected based on correlation values")
            inform02 <- paste("singlelinc: the number of neighbour",
                              "genes was reduced by 'coExprCut'")
            inform03 <- quote(paste("singlelinc: co-expression",
                                    "analysis yielded", ql_promise, "genes"))                              
            inform04 <- paste("singlelinc: Package 'clusterProfiler' was",
                              "called for gene annotation")
            inform05 <- paste("singlelinc: compatibility of gene ids",
                              "will be checked")
            
            # get information from 'linc' 
            validObject(input)  
            
            
            
            cm_promise  <- correlation(linCenvir(input)$linc)[[1]]
            out_history <- history(linCenvir(input)$linc)
            aq_promise  <- colnames(cm_promise)
            corMethod  <- out_history$corMethod
            cor_use    <- out_history$cor_use
            
            store       <- new.env(parent = emptyenv())
            # on.exit(print("Possible queries are:"),
            #        print(paste( aq_promise)))
            
            ## SECTION0: INPUT CHECK
            # interpretation of 'onlycor', 'underth, 'handleGeneIds'
            # and 'verbose'
            onlycor <- inlogical(onlycor, direct = FALSE)
            underth <- inlogical(underth, direct = TRUE)
            handleGeneIds <- inlogical(handleGeneIds, direct = TRUE)
            verbose <- inlogical(verbose, direct = TRUE)
            if(!verbose) message <- function(x) x
            
            # interpretation of 'testFun'
            if(!is.null(testFun)){
              tF_promise <- testFun
              if(!is.function(match.fun(
                tF_promise))) stop(errorm01)
              fF_promise <- names(formals(tF_promise))
              # formals x, y, method and use
              x_promise <- any(is.element(fF_promise, "x"))
              y_promise <- any(is.element(fF_promise, "y"))
              m_promise <- any(is.element(fF_promise, "method"))
              u_promise <- any(is.element(fF_promise, "use"))
              if(!all(x_promise, y_promise, u_promise,
                      m_promise) && !identical(tF_promise, stats::cor.test)) warning(warnim03)
              ff_promise <- try(testFun(c(1:10), c(1:10), method =
                                          corMethod, use = cor_use)$pvalue)
              if(class(ff_promise) == "try-error") stop(errorm02)
              test_call <- TRUE 
            } else {
              test_call <- FALSE; onlycor <- TRUE
            }
            
            #interpretation of 'coExprCut'
            if(!is.null(coExprCut)){
              ct_promise <- try(as.integer(coExprCut))
              if(class(ct_promise) == "try-error") stop(errorm03)
              if(length(ct_promise) != 1) stop(errorm03)
            } else {
              ct_promise <- 500
            }
            
            # interpretation of the 'threshold' argument
            th_promise <- threshold
            if(length(th_promise) != 1) stop(errorm04)
            if(!is.numeric(th_promise)) stop(errorm04)
            out_history$threshold <- th_promise
            if(th_promise >= 1 | th_promise < 0) warning(warnim01)
            
            # interpretation of query
            if(is.null(query)) stop(errorm05)
            qy_promise <- query
            if(length(qy_promise) != 1) stop(errorm06)
            found <- try(any(is.element(aq_promise, qy_promise)))
            if(class(found) == "try-error")  stop(errorm06)
            if(!found) stop(errorm06)
            
            # interpretation of annotateFrom
            cP_promise <- try(match.arg(annotateFrom,
                                        c("enrichGO",
                                          "enrichKEGG",
                                          "enrichPathway",
                                          "enrichDO")))
            if(class(cP_promise) == "try-error") stop(errorm02)
            if(any(is.element(cP_promise, c("enrichGO",
                                            "enrichKEGG",
                                            "enrichPathway",
                                            "enrichDO")))) {
              #suppressPackageStartupMessages(require(clusterProfiler))  
              #suppressPackageStartupMessages(require(org.Hs.eg.db))
              message(inform04)
              annotateFun <- match.fun(cP_promise)
            }
            
            ## SECTION1: CO-EXPRESSION AND GENE SELECTION
            # argument contradiction
            if(test_call && onlycor){
              warning(warnim02)
            }
            
            # in case only correlation defined
            if(onlycor){
              message(inform01); test_call <- FALSE
              pre_selected <- cm_promise[, qy_promise]
              ps_sort <- sort(pre_selected, decreasing = !underth) #!underth)
              if(!underth)selected <- ps_sort[ps_sort > th_promise] 
              if(underth) selected <- ps_sort[ps_sort < th_promise]  #
            }
            
            # in case a correlation test should be performed
            if(!onlycor){
              # do the test for the query gene
              pc_matrix  <- out_history$pc_matrix
              nc_query   <- out_history$nc_matrix[qy_promise, ]
              query_tested  <- apply(pc_matrix, 1, function(
                x = pc_matrix, nc_query, corMethod, cor_use){
                out <- tF_promise(x, nc_query, method =
                                    corMethod, use = cor_use, alternative, ...)
                return(out$p.value)
              }, nc_query, corMethod, cor_use)
              out_history$pval_min <- min(query_tested , na.rm = TRUE) 
              
              # select from the p-values
              ps_sort <- sort(query_tested, decreasing = !underth)
              if(!underth)selected <- ps_sort[ps_sort > th_promise] 
              if(underth) selected <- ps_sort[ps_sort < th_promise]
            }
            
            # do the final selection
            ql_promise <- length(selected)
            if(ql_promise == 0) stop(errorm07)
            if(ql_promise > ct_promise){
              selected <- selected[seq_len(ct_promise)]
              ql_promise <- ct_promise
            }
            message(eval(inform03))
            qg_promise <- names(selected)
            
            # define output for correlation and the test
            if(test_call){
              out_pval <- selected  
              out_corl <- cm_promise[names(selected), qy_promise]
            }
            if(!test_call){
              out_pval <- rep(NA, ql_promise)
              out_corl <- selected
            }
            
            ## SECTION2: GENE TRANSLATION
            gn_promise <- identifyGenes(qg_promise)
            # message(inform05)
            # issue: translate to other gene ids
            if(handleGeneIds == TRUE){
              # suppressPackageStartupMessages(require(mygene))  
              message(inform05)
              transThis <- function(x){
                query <- mygene::getGenes(x, fields = "entrezgene")
                entrez_query <- query@listData$entrezgene
                return(entrez_query)
              }
            }
            if(handleGeneIds == FALSE){
              if(!any(is.element(gn_promise, "entrezgene"))){
                warning(warnim04)
              }
              transThis <- function(x) x
            }
            
            OrgDb <- 'org.Hs.eg.db'
            if(cP_promise == "enrichPathway"){ OrgDb <- "human"}
            ## SECTION3: CALL TO GENE ANNOTATION
            entrez_query <- transThis(qg_promise)  
            cP_result    <- annotateFun(gene = entrez_query, OrgDb, ...)
            if(class(cP_result) == "enrichResult"){
              qvalues     <- slot(cP_result, "result")$qvalue
              bio_terms   <- slot(cP_result, "result")$Description
              out_query   <- list(qvalues, bio_terms)
            } else {
              out_query <- list(NA, NA)
            }
            
            ## SECTION3: PREPARATION OF OUTPUT
            # checkpoint reached, do not print queries
            print <- function(x) x
            
            out_linc             <- new("LINCsingle")
            results(out_linc)     <- list(query  = qy_promise,
                                         bio   = out_query,
                                         cor   = out_corl,
                                         pval  = out_pval)  
            assignment(out_linc)  <- assignment(input)
            correlation(out_linc) <- list(single = cm_promise[, qy_promise])
            express(out_linc)  <- express(input)
            history(out_linc)     <- out_history
            out_linCenvir        <- NULL
            out_linCenvir        <- linCenvir(input)
            out_linCenvir$single <- out_linc
            #lockEnvironment(out_linCenvir, bindings = TRUE)
            linCenvir(out_linc)   <- out_linCenvir
            return(out_linc)
          })


## getter methods
setGeneric("fcustomID", function(x) standardGeneric("fcustomID"))
setMethod("fcustomID", "LINCfeature", function(x) x@customID)

setGeneric("fcustomCol", function(x) standardGeneric("fcustomCol"))
setMethod("fcustomCol", "LINCfeature", function(x) x@customCol)

setGeneric("fsetLevel", function(x) standardGeneric("fsetLevel"))
setMethod("fsetLevel", "LINCfeature", function(x) x@setLevel)

setGeneric("fshowLevels", function(x) standardGeneric("fshowLevels"))
setMethod("fshowLevels", "LINCfeature", function(x) x@showLevels)

## setter methods

setGeneric("fcustomID<-", function(x, value) standardGeneric("fcustomID<-"))
setReplaceMethod("fcustomID", "LINCfeature",
                 function(x, value) {x@customID <- value; x})

setGeneric("fcustomCol<-", function(x, value) standardGeneric("fcustomCol<-"))
setReplaceMethod("fcustomCol", "LINCfeature",
                 function(x, value) {x@customCol <- value; x})

setGeneric("fsetLevel<-", function(x, value) standardGeneric("fsetLevel<-"))
setReplaceMethod("fsetLevel", "LINCfeature",
                 function(x, value) {x@setLevel <- value; x})

setGeneric("fshowLevels<-", function(x, value) standardGeneric("fshowLevels<-"))
setReplaceMethod("fshowLevels", "LINCfeature",
                 function(x, value) {x@showLevels <- value; x})

feature <- function(setLevel   = NULL,
                    customID   = NULL,
                    customCol  = "black",
                    showLevels = FALSE){
  out_feature  <- new("LINCfeature")
  linc_classes <- c("LINCmatrix", "LINCcluster",
                    "LINCbio", "LINCsingle")
  
  # argument 'setLevel'
  if(!is.null(setLevel)){
    if(isTRUE(tryCatch(expr = (is.element(setLevel,
                                          linc_classes))))){
      fcustomCol(out_feature) <- customCol
    } else {
      stop(paste("'setLevel' not one of:",
                 paste(linc_classes, collapse = ', ')))
    }  
    fsetLevel(out_feature) <- setLevel  
  }
  
  # argument 'customID'
  if(!is.null(customID)){
    fcustomID(out_feature)    <- customID
  }
  
  # argument 'customCol'
  if(isTRUE(tryCatch(expr = (is.element(customCol,
                                        colors()))))){
    fcustomCol(out_feature) <- customCol
  } else {
    stop("invalid color name for 'customCol'")
  }
  
  # argument 'showLevels'
  fshowLevels(out_feature) <- inlogical(showLevels,
                                      direct = FALSE)
  
  return(out_feature)
}

setMethod(f = '+',
          signature = c("LINCmatrix", "LINCfeature"),   
          definition  = function(e1, e2){
            validObject(e1)  
            # set the level of e1
            levels <- mget(c("cluster", "single", "bio", "linc"),
                           envir = linCenvir(e1), ifnotfound = FALSE)
            if(length(fsetLevel(e2)) > 0){
              if(!is.logical(levels$linc) && fsetLevel(e2) ==
                 "LINCmatrix") e1  <- linCenvir(e1)$linc
              if(!is.logical(levels$cluster) && fsetLevel(e2) ==
                 "LINCcluster") e1  <- linCenvir(e1)$cluster
              if(!is.logical(levels$bio) && fsetLevel(e2) ==
                 "LINCbio") e1  <- linCenvir(e1)$bio
              if(!is.logical(levels$single) && fsetLevel(e2) ==
                 "LINCsingle") e1  <- linCenvir(e1)$single
            }
            
            # add costum  IDs and colors to the history slot  
            if(length(fcustomID(e2)) > 0){
              history(e1)$customID <- fcustomID(e2)
            }
            if(length(fcustomCol(e2)) > 0){
              history(e1)$customCol <- fcustomCol(e2)
            }
            
            if(isTRUE(fshowLevels(e2))){  
              message("levels for this object:")
              lapply(levels, function(x){
                if(!is.logical(x)) message(class(x))
              })
            }
            
            return(e1)
          })

setMethod(f = '+',
          signature = c("LINCbio", "LINCfeature"),   
          definition  = function(e1, e2){
            callNextMethod()
          })

setMethod(f = '+',
          signature = c("LINCcluster", "LINCfeature"),   
          definition  = function(e1, e2){
            callNextMethod()
          })

## intersect methods
setMethod(f = '+',
          signature = c("LINCbio", "LINCbio"),   
          definition  = function(e1, e2){
            # add the costum id to the history slot  
            
            errorm01 <- "No match found for these objects"          
            
            b1_promise <- results(linCenvir(e1)$bio)[[1]]
            b2_promise <- results(linCenvir(e2)$bio)[[1]]
            query_both <- intersect(names(b1_promise),
                                    names(b2_promise))
            if(length(query_both) == 0) stop(errorm01)
            b1_index <- match(query_both, names(b1_promise))
            b2_index <- match(query_both, names(b2_promise))
            
            bio_intersect <- mapply(function(x, y){
              bio_terms <- unlist(intersect(x[[2]], y[[2]]))
              qvalues <- unlist(x[[1]][ match(bio_terms, x[[2]]) ])
              out <- list(qvalues, bio_terms)
              query_entry <- list(out)
              return(query_entry)
            } , x = b1_promise[b1_index],
            y = b2_promise[b2_index])
            
            results(e1) <- bio_intersect
            return(e1)
          })

setGeneric(name = "overlaylinc",
           def = function(
             input1,
             input2){
             standardGeneric("overlaylinc") # do it from environment
           })
setMethod(f   = "overlaylinc",
          signature = c("LINCbio",
                        "LINCbio"),  
          def = function(
            input1,
            input2){
            
            validObject(input1); validObject(input2)  
            e1e2_intersect <- (input1 + input2)          
            to_overlay <- results(e1e2_intersect)     


            cluster  <- results(linCenvir(input1$cluster))[[1]] 
            bio_list <- results(linCenvir(input1$bio))[[1]] 
            
            ## SECTION0: INPUT CONTROL  
            # check for a cluster
            
            #+ to be added          
            
            # plot the cluster
            tree <- ggtree(cluster, colour = "firebrick") +
              coord_flip() + scale_x_reverse()
            tree <- tree + geom_tiplab(size = 3.5, angle = 90,
                                       colour = "firebrick", hjust = 1)
            
            # prepare biological terms for plotting
            # preparation of a matrix
            
            # MAKE THIS ROBUST
            term_ext <- 1
            while(term_ext < 20){
              term_ext <- (term_ext + 1)
              raw_names <- names(bio_list)
              term_crude <- lapply(bio_list, function(x, term_ext){
                x[[2]][seq_len(term_ext)] }, term_ext)
              term_unique <- unique(unlist(term_crude))
              if(length(term_unique) > 20) break
            } 
            term_unique[is.na(term_unique)] <- "NA"
            m <- length(raw_names); n <- length(term_unique)
            first_matrix <- matrix(rep(0, (m*n)), ncol = n, nrow = m )
            colnames(first_matrix) <- term_unique
            rownames(first_matrix) <- raw_names
            
            # now fill matrix with biological terms
            for(query in seq_len(m)){
              if(length(bio_list[[query]][[2]]) > 0 &&
                 is.character(bio_list[[query]][[2]]) &&
                 is.numeric(bio_list[[query]][[1]])){
                
                bio_query <- bio_list[[query]][[2]]
                pvalues <- (-log10(bio_list[[query]][[1]]))
                row_entry <- vapply(colnames(first_matrix),
                                    function(x, pvalues){
                                      if(any(x == bio_query)){
                                        pvalues[x == bio_query][1]
                                      } else {
                                        return(0)   
                                      }
                                    }, 0, pvalues)
                first_matrix[query,] <- row_entry
              }
            }
            
            # additional plot for term assignments
            term_assign <- c(seq_len(length(term_unique)))
            bio_assignment <- mapply(function(x,y){ paste(x,y) },
                                     x = term_assign, y = term_unique)
            df_assign <- data.frame(bio_assignment, y = -term_assign,
                                    x = rep(0, length(term_unique)))
            
            pty_pl <- (ggplot(df_assign, aes(x,y)) +
                         geom_point(color = "white") + xlim(0, 1) +
                         theme(axis.line = element_line(colour =
                                                          "white"), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.border = element_blank(),   
                               panel.background = element_blank()) +
                         theme(axis.title.y = element_text(color =
                                                             "white")) + theme(axis.title.x = element_text(
                                                               color = "white")) + theme(axis.text.x = 
                                                                                           element_text(color = "white")) + theme(
                                                                                             axis.text.y = element_text(color = "white")))
            
            bio_assign_plot <- pty_pl + geom_text(aes(label = 
                                                        bio_assignment),hjust=0, vjust=0,
                                                  size = 4, colour = "grey18") 
            
            over_matrix <- first_matrix
            colnames(over_matrix) <- term_unique
            over_matrix[,] <- 0
            
            for(query in seq_len(length(to_overlay))){
              if(length(to_overlay[[query]][[2]]) > 0 &&
                 is.character(to_overlay[[query]][[2]]) &&
                 is.numeric(to_overlay[[query]][[1]])){
                
                bio_query <- to_overlay[[query]][[2]]
                pvalues <- (-log10(to_overlay[[query]][[1]]))
                
                row_entry <- vapply(colnames(over_matrix),
                                    function(x, pvalues){
                                      if(any(x == bio_query)){
                                        pvalues[x == bio_query][1]
                                      } else {
                                        return(0)   
                                      }
                                    }, 0, pvalues)
                over_matrix[query,] <- row_entry
              }
            }
            
            # convert to discrete values and join with tree
            # significant in the first and second?
            plot_matrix <- first_matrix
            plot_matrix[plot_matrix > 0 ] <- 1
            over_matrix[over_matrix > 0] <- 1
            plot_matrix <- ( plot_matrix + over_matrix) 
            colnames(plot_matrix) <- term_assign
            plot_df <- as.data.frame(plot_matrix)
            clust_heat <- gheatmap(tree, plot_df, offset = 0.1,
                                   width = 0.5, colnames = TRUE,
                                   colnames_position = "top")
            interbio_heat <- clust_heat + scale_fill_gradient2(aes(values),
                                                               low = "white", mid = "cornflowerblue", high =
                                                                 "darkviolet", midpoint = 1) + guides(fill=FALSE)
            
            interbio_img <- readPNG(system.file("extdata", "interbio_img.png",
                                                package ="LINC"))
            interbio_plot <- rasterGrob(interbio_img, interpolate = TRUE)
            
            # assembly of the final plot
            customid <- ""
            if(exists("customID", envir = history(input1))){
              customid <- history(input1)$customID
            } 
            
            # suppressPackageStartupMessages(require(gridExtra))     
            right_side <- arrangeGrob(interbio_plot, bio_assign_plot,
                                      ncol = 1, bottom = customid)
            final_plot <- grid.arrange(interbio_heat, right_side,
                                       nrow = 1)
            return(invisible(final_plot))  
            
          })

# set all validation methods
setValidity("LINCmatrix", method =
              function(object){
                if(any(is.element(c("linc", "cluster",
                                    "bio", "single"), ls(linCenvir(object)))
                )){
                  return(TRUE) 
                } else {
                  stop("Not a valid instance of the 'LINC' class ")
                }
              }
)  

# function getlinc
setGeneric(name = "getlinc",
           def = function(
             input,
             subject = "queries"){
             standardGeneric("getlinc")
           })
setMethod(f   = "getlinc",
          signature = c("ANY", "character"),  
          def = function(
            input,
            subject = "queries"){
            
            # check the class          
            if(!any(is.element(class(input), 
                               c("LINCmatrix", "LINCcluster",
                                 "LINCsingle", "LINCbio")))){
              stop("input is not of a supported 'LINC' class")
            }        
            
            # one of the subject arguments         
            sj_try  <- try(any(is.element(subject,
                                          c("queries", "geneSystem",
                                            "results", "history",
                                            "customID"))), silent = TRUE)  
            if(class(sj_try) == "try-error") stop("subject invalid")
            if(length(subject) != 1 | sj_try == FALSE)
              stop("subject invalid")
            
            # go truth all arguments
            if(subject == "history"){
              print(str(mget(ls(history(input)),
                             envir = history(input))))
              return(invisible((mget(ls(history(input)),
                                     envir = history(input)))))
            }
            
            if(subject == "queries"){
              return(listQuery(input))
            }    
            
            if(subject == "geneSystem"){
              message("diagnose for the gene system:")
              assignment_ids <- try(identifyGenes(assignment(input)),
                                    silent = TRUE)
              if(class(assignment_ids) != "try-error"){
                message("from slot 'assignment':", assignment_ids)  
              }
              expression_ids <- try(identifyGenes(
                rownames(express(input))),
                silent = TRUE)
              if(class(expression_ids) != "try-error"){
                message("from slot 'expression':", expression_ids)  
              }
            }
            
            if(subject == "customID"){
              if(exists("customID", envir = history(input))){
                return(history(input)$customID)
              } else {
                message("no custom id for this object")
              }
            }
            
            if(subject == "results"){
              print(str(results(input)))
              return(invisible(results(input)))
            }
          })

# function linctable
setGeneric(name = "linctable",
           def = function(
             file_name = "linc_table",
             input){
             standardGeneric("linctable")
           })
setMethod(f   = "linctable",
          signature = c("character", "LINCbio"),  
          def = function(
            file_name = "linc_table",
            input){
            pre <- results(input)
            tab <- lapply(pre[[1]], function(y){y[[2]][1:500] })
            m_tab <- matrix(unlist(tab), ncol = 500, nrow = length(tab), byrow = TRUE )
            rownames(m_tab) <- names(tab)
            write.table(m_tab, file = file_name, col.names = FALSE, sep = "\t"  )
            message("table of enriched terms written")
          })
setMethod(f   = "linctable",
          signature = c("character", "LINCcluster"),  
          def = function(
            file_name = "linc_table",
            input){
            pre <- results(input)[[1]]$neighbours
            tab <- lapply(pre, function(y){y[1:500] })
            m_tab <- matrix(unlist(tab), ncol = 500, nrow = length(tab), byrow = TRUE )
            rownames(m_tab) <- names(tab)
            write.table(m_tab, file = file_name, col.names = FALSE, sep = "\t"  )
            message("table of co-expressed genes written")
          })

## END OF SCRIPT