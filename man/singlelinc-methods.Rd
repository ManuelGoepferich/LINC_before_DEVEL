\name{singlelinc-methods}
\docType{methods}
\alias{singlelinc}
\alias{singlelinc-methods}
\alias{singlelinc,LINCmatrix-method}
\title{
 Co-Expression Analysis Of A Single ncRNA Gene
}
\description{
The function \code{singlelinc} perfroms co-expression analysis for a single query. An input \code{LINCmatrix} will be converted to a \code{LINCsingle} object. As a first step (I) a set of co-expressed protein-coding genes of a query is determined. Secondly, (II) biological terms related to the these genes are derived. The result will show the co-expression for the query.
}
\section{Methods}{
\describe{

\item{\code{signature(input = "LINCmatrix")}}{
(see details)
}
}}

\usage{
singlelinc(input,
           query,
        onlycor = FALSE,
        testFun = cor.test,
    alternative = "greater",
      threshold = 0.05,
        underth = TRUE,
      coExprCut,
  handleGeneIds = TRUE,
  annotateFrom  = 'enrichGO',
        verbose = TRUE,
        ...)
}

\arguments{
  \item{input}{
 an object of the class \code{LINCmatrix}
}
  \item{query}{
the name of the (ncRNA) gene to be evaluated; has to be present in \code{input}
}
  \item{onlycor}{
if \code{TRUE} co-expression will be decided based on absolute correlation values from the input \code{LINCmatrix} object. If \code{FALSE} co-expressed genes will be selected based on p-values from the correlation test.
}
  \item{testFun}{
a function to test the robustness of correlations. User-defined functions are allowed. The expected output is a p-value.
}
  \item{alternative}{
one of \code{c("two.sided", "less", "greater")}. This argument indicates the alternative in the correlation test. \code{"less"} instead of \code{"greater"} can be used for negative correlations.
}
  \item{threshold}{
a single number representing the threshold for selecting co-expressed genes
}
  \item{underth}{
if \code{TRUE} values lower than the threshold will be considered (intended for p-values). If \code{FALSE} values higher tahn the threshold will be considered (intended for absolute correlations)     
}
  \item{coExprCut}{
a single \code{integer} indicating the maximal number of co-expressed genes to select. In case too many genes fulfill the supplied \code{threshold} criterion their number can be reduced by this argument.
}
  \item{handleGeneIds}{
if \code{TRUE} attempt to recognize and handle issues with incompatible gene systems. If \code{FALSE} skip this section.     
}
  \item{annotateFrom}{
a function which will derive significant biological terms based on the set of co-expressed genes from a gene annotation resource.   
}
  \item{verbose}{
whether to give messages about the progression of the function \code{TRUE} or not \code{FALSE}
}
  \item{...}{
further arguments, mainly for \code{cor.test} and functions from {clusterProfiler}
}
}

\details{
In comparison to the function \code{clusterlinc} this function will provide more flexibility in terms of the selection of co-expressed genes. The option \code{onlycor = TRUE} in combination with a suitable \code{threshold} can be used to choose co-expressed protein-coding genes based on the correlation values inherited from the input \code{LINCmatrix} object. For this to work it is required to set \code{underth = FALSE} because then, values higher than the \code{threshold} will be picked. By default, co-expression depnds on the p-values from the correlation test (\code{stats::cor.test}) which demonstrate the robustness of a given correlation between two genes. A user-defined test function supplied in \code{testFun} requires the formal arguments \code{x}, \code{y}, \code{method} and \code{use}. Moreover, the p-values of the output should be accessible by \code{$pvalue}. The number of co-expressed genes can be restricted not only by \code{threshold}, but also by \code{coExprCut}. The value \code{n} for \code{coExprCut = n} will be ignored in case the number of genes which fulfill the \code{threshold} criterion is smaller than \code{n}.   
\code{handleGeneIds = FALSE} should only be used for an input with Entrez gene ids. Options for \code{annotateFrom} are for example: \code{ReactomePA::enrichPathway()}, \code{clusterProfiler::enrichGO} or \code{clusterProfiler::enrichKEGG}. Further arguments (\code{...}) are inteded to be passed to the called \code{annotateFrom} function. \code{annotateFrom = 'enrichGO', ont = "CC"} will call the subontology "Cellular Component" from GO. For more details see the documentation of\code{clusterProfiler}.  
}

\value{
an object of the class 'LINCmatrix' (S4) with 6 \code{Slots}
  \item{results}{a \code{list} of four entries: \code{$query}, the queried gene, \code{$bio}, a list of biological terms and their p-values, \code{$cor}, absolute correlations of co-expressed genes, \code{$pval}, p-values from the correlation test of co-expressed genes 
  }
    \item{assignment}{a \code{character} vector of protein-coding genes}
  \item{correlation}{ a \code{list} of \code{$single}, the correlation of the query to protein-coding genes}
  \item{expression}{the original expression matrix}
  \item{history}{a storage environment of important methods, objects and parameters used to create the object}
  \item{linCenvir}{a storage environment ensuring the compatibility to other objects of the \code{LINC} class}
}

\section{Compatibility}{
\code{plotlinc(LINCsingle, ...))}, ...  
}

\author{
Manuel Goepferich
}

\seealso{
\code{linc} ; 
\code{singlelinc}
}
\examples{
data(BRAIN_EXPR)

# selection based on absolute correlation
meg3 <- singlelinc(crbl_matrix, query = "55384", onlycor = TRUE, underth = FALSE, threshold = 0.5)
plotlinc(meg3)

# get biological processes from 'Gene Ontology'
meg3 <- singlelinc(crbl_matrix, query = "55384", threshold = 0.00005, ont = 'BP')
plotlinc(meg3)

# get enriched pathways from 'ReatomePA'
meg3 <- singlelinc(crbl_matrix, query = "55384", onlycor = TRUE, underth = FALSE, threshold = 0.5,
                   annotateFrom = 'enrichPathway')
plotlinc(meg3)
}  
\keyword{ ~singlelinc }
\keyword{methods}
