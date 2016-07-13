\name{LINCmatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{LINCmatrix-class}

\title{Class \code{"LINCmatrix"}}
\description{
an object of the class 'LINCmatrix' (S4) with 6 \code{Slots}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LINCmatrix", ...)}.
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
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "LINCmatrix", e2 = "LINCfeature")}: ... }
    \item{clusterlinc}{\code{signature(linc = "LINCmatrix")}: ... }
    \item{linc}{\code{signature(object = "LINCmatrix", codingGenes = "missing")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCmatrix", showCor = "character")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCmatrix", showCor = "missing")}: ... }
    \item{singlelinc}{\code{signature(input = "LINCmatrix")}: ... }
	 }
}

\examples{
showClass("LINCmatrix")
}
\keyword{classes}
