\name{justlinc-methods}
\docType{methods}
\alias{justlinc-methods}
\alias{justlinc,matrix-method}
\alias{justlinc}
\title{
 Co-Expression Analysis Of ncRNA Genes In One Step
}
\description{
\code{justlinc} is a wrapper that uses different functions of the \pkg{LINC} package. It applies fixed thresholds for gene selection and requires only an (unprocessed) expression matrix as input. This enables a fast co-expression analysis with or without a list candidates.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "matrix")}}{
(see details)
}
}}

\usage{
justlinc(object,
    targetGenes = "lincRNA",
           rmPC = TRUE)
}

\arguments{
  \item{object}{
 a gene expresssion matrix with Ensembl gene ids (rows correspond to genes)
}
  \item{targetGenes}{
the gene biotype like "lincRNA" or a vector of querie(s) - one or multiple genes ids 
}
  \item{rmPC}{
remove 30 percent of the variance in the data by PCA \code{TRUE} or not \code{FALSE}.
}
}

\details{
This function was built for expression matrices which uses the Ensembl gene system with gene expression values for over 50.000 genes. The input will be matched with a static gene annotation. Genes will be selected for median expression and variance. The final correlation matrix of protein-coding genes versus the target genes (lincRNAs) considers 5000 protein-coding genes and 500 target genes. Enriched pathways will be computed for the co-expressed genes showing the highest Spearman's rank correlation.
In case \code{targetGenes = "lincRNA"}, then the function will search for the 10 best lincRNAs. Supplying a vector of gene ids the method will determine the best co-expressed genes for the given queries.
Importanly, this function is only a wrapper. Thresholds can be changed using \code{linc} and related functions.
}

\value{
depending on the input in \code{targetGenes} one or two plots and the result of the co-expression analysis
  }

\author{
Manuel Goepferich
}

\seealso{
\code{linc} 
}
\examples{
data(LIVER_EXPR)

# a gene expression matrix with > 50,000 genes
str(GTEX_LIVER_CRUDE)

# 'justlinc' will search for the 10 best candidates
justlinc(GTEX_LIVER_CRUDE)

# 'justlinc' called with queries
my_lincRNAs <- c("ENSG00000224153", "ENSG00000197813",
                 "ENSG00000179136", "ENSG00000259439",
                 "ENSG00000267462")

res <- justlinc(GTEX_LIVER_CRUDE, targetGenes = my_lincRNAs)
}  
\keyword{ ~justlinc }
\keyword{methods}