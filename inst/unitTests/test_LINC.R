# Unit Tests for the R package

#require(RUnit)

# 1 computation of correlation
test_corruptedMatrix <- function() {
  
 cor_test_mat <- matrix(c(c(1:110), NA, NA, 50, Inf, 3,
 NaN, 0.5, NA, 80, 5, 4, c(1:50), rep(NA, 14), rep(9, 15)),
 ncol = 10, byrow = TRUE)
 rownames(cor_test_mat) <- as.character(c(1:20))
 colnames(cor_test_mat) <- LETTERS[1:10]
  
 ref <- cor(cor_test_mat["13", ], cor_test_mat["2", ],
 use =  "pairwise.complete.obs", method = "spearman")
 res <- linc(cor_test_mat, codingGenes = c(TRUE, FALSE,
 TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
 TRUE, TRUE, rep(FALSE, 7)))
 val <- res@correlation[[1]]["13", "2"] 
 
 checkEqualsNumeric(val, ref, tolerance = 1.0e-2)
}

# 2 finding queries in the 'LINCcluster' object
test_QueriesofCluster <- function() {
  
  cor_test_mat <- matrix(c(c(1:110), NA, NA, 50, Inf, 3,
                           NaN, 0.5, NA, 80, 5, 4, c(1:50), rep(NA, 14), rep(9, 15)),
                         ncol = 10, byrow = TRUE)
  rownames(cor_test_mat) <- as.character(c(1:20))
  colnames(cor_test_mat) <- LETTERS[1:10]
  
  linc_matrix <- linc(cor_test_mat, codingGenes = c(TRUE, FALSE,
                                            TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                            TRUE, TRUE, rep(FALSE, 7)))
 res <-  clusterlinc(linc_matrix)
 checkEquals(names(res@results$cluster$neighbours), colnames(linc_matrix@correlation[[1]]))
}

# 3 test for minimal p-values
test_pvalueCluster <- function() {
  
  cor_test_mat <- matrix(c(c(1:110), NA, NA, 50, Inf, 3,
                           NaN, 0.5, NA, 80, 5, 4, c(1:50), rep(NA, 14), rep(9, 15)),
                         ncol = 10, byrow = TRUE)
  rownames(cor_test_mat) <- as.character(c(1:20))
  colnames(cor_test_mat) <- LETTERS[1:10]
  
  linc_matrix <- linc(cor_test_mat, codingGenes = c(TRUE, FALSE,
                                                    TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                                                    TRUE, TRUE, rep(FALSE, 7)))
  res <-  clusterlinc(linc_matrix)
  checkTrue(all( res@correlation$cortest >= 1e-26))
}

# 4 queries in singlelinc
#test_QueryofSingle <- function() {
  
#  cor_test_mat <- matrix(c(c(1:110), NA, NA, 50, Inf, 3,
#                           NaN, 0.5, NA, 80, 5, 4, c(1:50), rep(NA, 14), rep(9, 15)),
#                         ncol = 10, byrow = TRUE)
#  rownames(cor_test_mat) <- as.character(c(1:20))
#  colnames(cor_test_mat) <- LETTERS[1:10]
  
#  linc_matrix <- linc(cor_test_mat, codingGenes = c(TRUE, FALSE,
#                                                    TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
#                                                    TRUE, TRUE, rep(FALSE, 7)))
#  res <-  singlelinc(linc_matrix, query = "17", onlycor = T, underth = F, threshold = 0.99)
#  checkEquals(res@results$query, "17")
#}

