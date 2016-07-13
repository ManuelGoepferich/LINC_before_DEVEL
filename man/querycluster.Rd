\name{querycluster}
\alias{querycluster}

\title{
 Cluster One ncRNA Gene Based On Its Co-Expression in Multiple Datasets
}
\description{
\code{querycluster} takes a set of \code{'LINCcluster'} objects, extracts the respective co-expressed protein-coding genes and plots a dendrogram with the distance matrix attached. This function is intended to be applied in a case where a particular ncRNA occurrs in datasets which represent different tissues, batches, statistical corrections, reduced gene sets, controls and so on. The output will show the clustering of the groups and therefore the information under which condition is the co-expression to the query most similar.  

}
\usage{
querycluster(query = NULL,
             queryTitle = NULL,
             traits = NULL, 
             method = "spearman",
             returnDat = FALSE,
             mo_promise,
             ...)
}

\arguments{
\item{query}{
the query name, i.e. the gene id of a ncRNA present in the supplied input 
}
\item{queryTitle}{
a \code{character} string used as the title of the plot
}
  \item{traits}{
\code{NULL} or a single {integer}. For \code{NULL} all co-expressed genes will be used. A number instead will be considered as maximal number of \code{traits}.    
}
  \item{method}{
a \code{character} string, has to be one of \code{c("spearman", "dicedist")}
}
  \item{returnDat}{
whether to return the data used to create the plot \code{TRUE} or not \code{FALSE}
}
  \item{mo_promise}{
\code{mo_promise = 'list'}, a list of \code{'LINCcluster'} objects (see example)
}
  \item{...}{
the \code{'LINCcluster'} objects itself, but not a combination of both, \code{mo_promise = 'list'} and supplying the objects itself (see example)
}
}


\details{
This function will search for co-expressed protein-coding genes which belong to a the defined \code{query}. Based on the co-expression in the input \code{'LINCcluster'} objects a distance matrix is computed. The method \code{"spearman"} finds the union of all interaction partners for the query und calculates the correlation between the \code{'LINCcluster'} objects. For this method the distance measure is \code{(1 - correlation)}. Alternatively, \code{method = "dicedist"} takes the Czekanovski dice distance [1] as distance mesaure of the \code{traits = n} genes. This method, however, will not work with \code{traits = NULL}. Choosing a low number for \code{n} will limit the number of different values in the distance matrix.

Apart from \code{queryTitle} the command \code{'LINCcluster' + feature(customID = ..., customCol = ...)} enables a customized plot as output. For this to work the supplied \code{'LINCcluster'} objects in \code{...} have to be modified by the function \code{feature(...)} in advance.

}

\value{
an object of the class 'gg' or a 'list' if \code{returnDat = TRUE}
}
\references{
[1] Christine Brun, Francois Chevenet, David Martin, Jerome Wojcik, Alain Guenoche and Bernard Jacq" Functional classification of proteins for the prediction of cellular function from a protein-protein interaction network" (2003) Genome Biology, 5:R6.
}
\author{
Manuel Goepferich
}

\seealso{
\code{feature} ; 
\code{clusterlinc}
}
\examples{
data(BRAIN_EXPR)

# add custom names and colors
gbm_cluster <- gbm_cluster + feature(customID = "CANCER_GBM", customCol = "red")
ctx_cluster <- ctx_cluster + feature(customID = "HEALTHY_CTX", customCol = "blue")
hpc_cluster <- hpc_cluster + feature(customID = "HEALTHY_HPC", customCol = "blue")
crbl_cluster <- crbl_cluster + feature(customID = "HEALTHY_CRBL", customCol = "blue")

# plot the dendrogram
querycluster('647979', queryTitle = 'NORAD',
             gbm_cluster,  # Glioblastoma
             ctx_cluster,  # Cortex
             hpc_cluster,  # Hippocampus
             crbl_cluster) # Cerebellum

# objects can also be supplied as a list
query_list <- list(gbm_cluster,
                   ctx_cluster,
                   hpc_cluster,
                   crbl_cluster)

# mo_promise is the (informal) argument for multiple objects
querycluster(query = '647979', queryTitle = 'NORAD', mo_promise = query_list)

# used the Czekanovski dice distance based on the 25 best 
# interaction partners in each tissue
querycluster(query = '647979', method = "dicedist", traits = 25, mo_promise = query_list)

# NOT RUN:
# querycluster(query = '647979', method = "dicedist", mo_promise = query_list)
}  
\keyword{ ~querycluster }
