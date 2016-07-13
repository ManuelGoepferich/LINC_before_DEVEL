\name{linc-methods}
\docType{methods}
\alias{linc-methods}
\alias{linc,data.frame,ANY-method}
\alias{linc,ExpressionSet,ANY-method}
\alias{linc,LINCmatrix,missing-method}
\alias{linc,matrix,ANY-method}
\alias{linc}

\title{
 Compute A Correlation Matrix of Co-expressed Coding And Non-Coding Genes
}
\description{
The function \code{linc} can be considered as the main function of this package. It converts a given input object into a \code{LINCmatrix}. This process includes (I) statistical analysis and (II) correction of the input, (III) separation of coding and non-coding genes and (IV) computation of a correlation matrix.The input could be for instance a gene expression matrix. Rows correspond to genes; columns represent samples.Besides a suitable \code{object} a vector identifying the protein-coding genes is required. 

}
\section{Methods}{
\describe{

\item{\code{signature(object = "data.frame", codingGenes = "ANY")}}{
(see details)
}

\item{\code{signature(object = "ExpressionSet", codingGenes = "ANY")}}{
(see details)
}

\item{\code{signature(object = "LINCmatrix", codingGenes = "missing")}}{
(see details)
}

\item{\code{signature(object = "matrix", codingGenes = "ANY")}}{
(see details)
}
}}

\usage{
linc(object,
     codingGenes,
     corMethod    = "spearman",
     batchGroups,
     nsv          = 1,
     rmPC,
     outlier,
     userFun,
     verbose      = TRUE)
}

\arguments{
  \item{object}{
 a \code{matrix}, \code{data.frame} or \code{ExpressionSet} with genes corresponding to rows, preferentially the high-variance genes in a given set
}
  \item{codingGenes}{
a \code{logical} vector with the same length of the supplied genes in \code{object}. \code{TRUE} indicates that the gene is a protein-coding one. Alternatively, \code{codingGenes} can be a vector naming the protein-coding genes.
}
  \item{corMethod}{
a method for the correlation function; has to be one of \code{c("pearson", "kendall", "spearman", "user")}
}
  \item{batchGroups}{
a vector naming the batch conditions. The length of this vector has to match the number of samples supplied in \code{object}. There has to be at least two different batch conditions for the method to work.
}
  \item{nsv}{
a single \code{integer} indicating the number of hidden surrogate variables. This argument is only relevant in case \code{batchGroups} is used. 
}
  \item{rmPC}{
a vector of principle components (PCs) which should be removed. PCs are counted staring from 1 up to the maximal count of samples.  
}
  \item{outlier}{
a method for the genewise removal of single outliers; has to be one of \code{c("esd", "zscore")} 
}
  \item{userFun}{
a function or its name that should be used to calculate the correlation between coding and non-coding genes. This argument has to be used in combination with \code{corMethod = "user"} 
}
  \item{verbose}{
whether to give messages about the progression of the function \code{TRUE} or not \code{FALSE}
}

}
\details{
\code{object} can be a \code{matrix}, a \code{data.frame} or an \code{ExpressionSet} with rows 
corresponding to genes and columns to samples, the assumed co-expression conditions. Genes with duplicated names, genes having 0 variance plus genes with to many missing or infinite values will be removed from the input. For inputs showing a high inter-sample variance (ANOVA) in combination with many single outliers a warning message will appear.
By default Spearman's rank correlation will be computed between protein-coding to non-coding genes. For this method a time-efficient C++ implementation will be called. Longer computation times occur for genes > 5000 and samples > 100. Missing values are handled in a manner that only pairwise complete observations will be compared. A customized correlation function can be applied supplying the function in \code{userFun} and requires the formal arguments \code{x} and \code{x}. This has priority over \code{corMethod}.

A number of statistical methods are available in order to remove effects from a given input expression matrix which depend on the used platform or technology and may hide relevant biology.
The argument \code{batchGroups} works as a rapper of the SVA package calling \code{sva::svaseq}. The number of hidden surrogate variables is set to \code{nsv = 1} by default; it can be estimated utilizing the function \code{sva::num.sv}. For this model to work the description of at least two different batches are required in \code{batchGroups}.
Principle Component Analysis (PCA) can be performed by \code{rmPC = c(...)} where \code{...} represents a vector of principle components. The command \code{rmPC = c(2:ncol(object))} will remove the first PC from the input. This method can be used to determine whether observations are due to the main variance in the dataset i.e. main groups or subtypes.
Outliers are handled genewise. The extreme Studentized deviate (ESD) test by  Rosner, Bernard (1983) will detect one up to four outliers in a gene and replace them by \code{NA} values. The alternative \code{zscore} will perform a robust zscore test suggested by Boris Iglewicz and David Hoaglin (1993) and detect a single outlier in a gene if |Z| > 3.5.

A \code{LINCmatrix} can be recalculated with the command \code{linc(LINCmatrix, ...))} in order to change further arguments. \code{plotlinc(LINCmatrix, ...))} will plot a figure depicting the statistical analysis and correlation values. As for most objects of the \code{LINC} class manipulation of the last slot \code{linCenvir} will likely result in unexpected errors.
}

\value{
an object of the class 'LINCmatrix' (S4) with 6 \code{Slots}
  \item{results}{a \code{list} containing the original input expression matrix or a transfomed matrix if \code{rmPC}, \code{batchGroups} or  \code{outlier} was applied}
    \item{assignment}{a \code{character} vector of protein-coding genes}
  \item{correlation}{ a \code{list} of \code{$cormatrix}, the correlation of non-coding to protein-coding genes and 
  \code{$lnctolnc}, the correlation of non-coding to non-coding genes}
  \item{expression}{the original expression matrix}
  \item{history}{a storage environment of important methods, objects and parameters used to create the object}
  \item{linCenvir}{a storage environment ensuring the compatibility to other objects of the \code{LINC} class}
}

\section{Compatibility}{
\code{plotlinc(LINCmatrix, ...)}, \code{clusterlinc(LINCmatrix, ...)}, \code{singlelinc(LINCmatrix, ...)}, ...  
}

\references{
[1] https://www.bioconductor.org/packages/release/bioc/html/sva.html

[2] Rosner, Bernard (May 1983), Percentage Points for a Generalized ESD Many-Outlier Procedure,Technometrics, 25(2), pp. 165-172. 

[3] Boris Iglewicz and David Hoaglin (1993), Volume 16:How to Detect and Handle Outliers", The ASQC Basic References in Quality Control: Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
}
\author{
Manuel Goepferich
}

\seealso{
\code{justlinc} ; 
\code{clusterlinc} ; 
\code{singlelinc}
}
\examples{
data(BRAIN_EXPR)

# call 'linc' with no further arguments
crbl_matrix <- linc(cerebellum, codingGenes = pcgenes_crbl)

# remove first seven principle components
crbl_matrix_pc <- linc(cerebellum, codingGenes = pcgenes_crbl, rmPC = c(1:7))

# negative correlation by using 'userFun'
crbl_matrix_ncor <- linc(cerebellum, codingGenes = pcgenes_crbl,
                         userFun = function(x,y){ -cor(x,y) })

# remove outliers using the ESD method
crbl_matrix_esd <- linc(cerebellum, codingGenes = pcgenes_crbl, outlier = "esd")

# plot this object
plotlinc(crbl_matrix_esd)
}  

\keyword{methods}
\keyword{linc}