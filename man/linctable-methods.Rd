\name{linctable-methods}
\docType{methods}
\alias{linctable-methods}
\alias{linctable,character,LINCbio-method}
\alias{linctable,character,LINCcluster-method}
\alias{linctable}
\title{Write To Table \pkg{LINC}}
\description{
 Write and save a table in \pkg{LINC}
}
\section{Methods}{
\describe{

\item{\code{signature(file_name = "character", input = "LINCbio")}}{
Table of biological terms.
}

\item{\code{signature(file_name = "character", input = "LINCcluster")}}{
Table of co-expressed genes
}
}}

\examples{
# write to a table
#linctable(file_name = "crbl_co_expr", input = crbl_cluster)
}  

\keyword{methods}
\keyword{ ~linctable }