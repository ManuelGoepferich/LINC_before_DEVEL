\name{LINCcluster-class}
\Rdversion{1.1}
\docType{class}
\alias{LINCcluster-class}

\title{Class \code{"LINCcluster"}}
\description{
code{"LINCcluster"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LINCcluster", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{results}:}{Object of class \code{"list"} ~~ }
    \item{\code{assignment}:}{Object of class \code{"vector"} ~~ }
    \item{\code{correlation}:}{Object of class \code{"list"} ~~ }
    \item{\code{expression}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{history}:}{Object of class \code{"environment"} ~~ }
    \item{\code{linCenvir}:}{Object of class \code{"environment"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{LINCmatrix}"}, directly.
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "LINCcluster", e2 = "LINCfeature")}: ... }
    \item{clusterlinc}{\code{signature(linc = "LINCcluster")}: ... }
    \item{getbio}{\code{signature(cluster = "LINCcluster")}: ... }
    \item{linctable}{\code{signature(file_name = "character", input = "LINCcluster")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCcluster", showCor = "character")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCcluster", showCor = "missing")}: ... }
	 }
}

\examples{
showClass("LINCcluster")
}
\keyword{classes}
