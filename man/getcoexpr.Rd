\name{getcoexpr}
\alias{getcoexpr}

\title{
 Get IDs For Co-Expressed Genes from The 'LINC' Class
}
\description{
\code{getcoexpr} provides access to co-expressed genes of a query in \code{'LINC'} objects.
}
\usage{
getcoexpr(input,
          query   = NULL)
}

\arguments{
  \item{input}{
a \code{'LINCcluster'} or \code{'LINCsingle'} object
}
  \item{query}{
for a \code{'LINCcluster'} the name of the query gene
}
}

\value{
a vector containing the co-expressed genes for a query
}

\author{
Manuel Goepferich
}


\examples{
data(BRAIN_EXPR)

# Get the co-expressed genes for the gene with the entrez id "441951" (ZFAS1)
getcoexpr(crbl_cluster, "441951" )
}  
\keyword{ ~getcoexpr }
