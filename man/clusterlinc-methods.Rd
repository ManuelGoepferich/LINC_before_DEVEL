\name{clusterlinc-methods}
\docType{methods}
\alias{clusterlinc}
\alias{clusterlinc-methods}
\alias{clusterlinc,LINCcluster-method}
\alias{clusterlinc,LINCmatrix-method}
\title{
 Cluster Queried ncRNAs Based On Their Interaction Partners
}
\description{
The function \code{clusterlinc} will give an overview of ncRNAs in a dataset. An input \code{LINCmatrix} will be converted to a \code{LINCcluster}. The following steps are conducted (I) computation of a correlation test, (II) setup of a distance matrix, (III) calculation of a dendrogram and (IV) selection of co-expressed genes for each query. The result is a cluster of ncRNAs and their associated protein-coding genes.
}
\section{Methods}{
\describe{

\item{\code{signature(linc = "LINCcluster")}}{
(see details)
}

\item{\code{signature(linc = "LINCmatrix")}}{
(see details)
}

}}

\usage{
clusterlinc(linc,
            distMethod    = "dicedist",
            clustMethod   = "average",
            pvalCutOff    = 0.05,
            padjust       = "none",
            coExprCut     = NULL,
            cddCutOff     = 0.05,
            verbose       = TRUE)
}

\arguments{
  \item{linc}{
 an object of the class \code{LINCmatrix}
}
  \item{distMethod}{
a method to compute the distance between ncRNAs; has to be one of \code{c("correlation", "pvalue", "dicedist")}
}
  \item{clustMethod}{
an algorithm to compute the dendrogram, has to be one of \code{c("ward.D",
                "ward.D2", "single", "complete", "average",
                "mcquitty", "median", "centroid")}
}
  \item{pvalCutOff}{
a threshold for the selection of co-expressed genes. Only protein-coding genes showing a significance in the correlation test lower than \code{pvalCutOff} will assigned to queried ncRNAs as interaction partner.
}
  \item{padjust}{
one of \code{c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}
}
  \item{coExprCut}{
a single \code{integer} indicating the number of co-expressed genes to select. If this argument is used for each ncRNA the \code{coExprCut = n} protein-coding genes with the lowest p-value in the correlation test will be assigned to queries, respectively. 
}
  \item{cddCutOff}{
a threshold that is only relevant for \code{distMethod = "dicedist"}. In this method \code{cddCutOff} defines whether a ncRNA and a protein-coding gene can be considered as interaction partners. This influences the distance matrix and the clustering process.   
}
  \item{verbose}{
whether to give messages about the progression of the function \code{TRUE} or not \code{FALSE}
}
}

\details{
As a first step \code{clusterlinc} conducts the correlation test (\code{stats::cor.test}) using the correlation method and handeling of missing values inherited from the input \code{LINCmatrix}. Resulting p-values indicate the statistical robustness of correlations instead of absolute correlation values. Co-expression of ncRNAs to protein-coding genes is assumed if the p-value from the \code{cor.test} is lower than the given \code{pvalCutOff}. An alternative way to select co-expressed genes is provided by \code{coExprCut}. This argument has priority over \code{pvalCutOff} and can be used to pick the \code{n} genes with the lowest p-value for each ncRNA. In contrast to \code{pvalCutOff}, this will result in an equal number of assigned co-expressed genes. The argument \code{padjust} can be used for multiple testing correction. In most cases this is not compatible with \code{distMethod = "dicedist"}.  

For the computation of the distance matrix of ncRNA genes three methods can be applied. The first method \code{"correlation"} uses \code{1 - correlation} as distance measure. In contrast, \code{"pvalue"} considers not the absolute correlation values, but p-values from the correlation test. A third method is termed \code{"dicedist"} and takes the Czekanovski dice distance [1] as distance measure. Here, the number of shared interaction partners between ncRNAs determines their relation to each other. The argument \code{cddCutOff} is an option to decide which p-values in the correlation matrix can be considered as interaction. A low threshold, for instance, will consider only interactions of ncRNAs and protein-coding genes supported by a p-value lower than the supplied threshold and therfore a robust correlation of these two genes. Based on the distance matrix a cluster of the ncRNAs will be computed by \code{stats::hclust}. Argument \code{clustMethod} defines which clustering method should be applied. 

A \code{LINCcluster} can be recalculated with the command \code{clusterlinc(LINCcluster, ...))} in order to change further arguments. \code{plotlinc(LINCcluster, ...))} will plot a figure that shows the cluster of ncRNAs (dendrogram) and the number of co-expressed genes with respect to different thresholds. \code{getbio(LINCcluster, ...))} will derive the biological terms associated with the co-expressed genes. Due to the correlation test longer calculation times can occur. A faster alternative to this function is \code{singlelinc()}. User-defined correlation test functions are supported for \code{singlelinc()} but not for \code{clusterlinc()}. 
}

\value{
an object of the class 'LINCmatrix' (S4) with 6 \code{Slots}
  \item{results}{a \code{list} containing an object of the class \code{"phylo"} with the additional entry \code{neighbours}, a list of queries and co-expressed genes}
    \item{assignment}{a \code{character} vector of protein-coding genes}
  \item{correlation}{ a \code{list} of \code{cormatrix}, the correlation of non-coding to protein-coding genes,
  \code{lnctolnc}, the correlation of non-coding to non-coding genes and
  \code{cortest}, p-values of the correlation test of non-coding to protein-coding genes}
  \item{expression}{the original expression matrix}
  \item{history}{a storage environment of important methods, objects and parameters used to create the object}
  \item{linCenvir}{a storage environment ensuring the compatibility to other objects of the \code{LINC} class}
}

\section{Compatibility}{
\code{plotlinc(LINCcluster, ...)}, \code{getbio(LINCcluster, ...)}, ...  
}

\references{
[1] Christine Brun, Francois Chevenet, David Martin, Jerome Wojcik, Alain Guenoche and Bernard Jacq" Functional classification of proteins for the prediction of cellular function from a protein-protein interaction network" (2003) Genome Biology, 5:R6.
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
class(crbl_matrix)

# call 'clusterlinc' with no further arguments
crbl_cluster <- clusterlinc(crbl_matrix)

# apply the distance method "correlation instead of "dicedist"
crbl_cluster_cor <- clusterlinc(crbl_matrix, distMethod = "correlation" )
# do the same as recursive call using the 'LINCcluster' object
# crbl_cluster_cor <- clusterlinc(crbl_cluster, distMethod = "correlation")

# select 25 genes with lowest p-values for each query
crbl_cluster_25 <- clusterlinc(crbl_matrix, coExprCut = 25)

# select onyl those with a p-value < 5e-5
crbl_cluster_5e5 <- clusterlinc(crbl_matrix, pvalCutOff = 5e-5)

# adjust for multiple testing
crbl_cluster_hochberg <- clusterlinc(crbl_matrix, distMethod = "correlation",
                                   padjust = "hochberg", pvalCutOff = 0.05)

# comparing two distance methods
plotlinc(crbl_cluster)
plotlinc(crbl_cluster_cor)
}  
\keyword{ ~clusterlinc }
\keyword{methods}
