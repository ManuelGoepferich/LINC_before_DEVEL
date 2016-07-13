\name{Arith-methods}
\docType{methods}
\alias{+}
\alias{Arith}
\alias{+-methods}
\alias{+,LINCbio,LINCbio-method} 
\alias{+,LINCbio,LINCfeature-method}
\alias{+,LINCcluster,LINCfeature-method}
\alias{+,LINCmatrix,LINCfeature-method}


\title{Plus Operator '+' In \pkg{LINC}}
\description{
\code{plotlinc} is the generic function that creates plots based on objects of the \code{LINC} class.
}
\section{Methods}{
\describe{

\item{\code{signature(e1 = "LINCbio", e2 = "LINCbio")}}{
(biological term intersection)
}

\item{\code{signature(e1 = "LINCbio", e2 = "LINCfeature")}}{
(see \code{feature()})
}

\item{\code{signature(e1 = "LINCcluster", e2 = "LINCfeature")}}{
(see \code{feature()})
}

\item{\code{signature(input = "LINCmatrix", showCor = "LINCfeature")}}{
(see \code{feature()})
}

}}

\examples{
data(BRAIN_EXPR)

# add a custom name
crbl_cluster_feat <- crbl_cluster + feature(customID = "CEREBELLUM")
plotlinc(crbl_cluster_feat)

}  

\keyword{methods}

