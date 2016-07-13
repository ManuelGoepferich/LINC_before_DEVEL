\name{LINCbio-class}
\Rdversion{1.1}
\docType{class}
\alias{LINCbio-class}

\title{Class \code{"LINCbio"}}
\description{
\code{"LINCbio"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LINCbio", ...)}.
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
    \item{+}{\code{signature(e1 = "LINCbio", e2 = "LINCbio")}: ... }
    \item{+}{\code{signature(e1 = "LINCbio", e2 = "LINCfeature")}: ... }
    \item{linctable}{\code{signature(file_name = "character", input = "LINCbio")}: ... }
    \item{overlaylinc}{\code{signature(input1 = "LINCbio", input2 = "LINCbio")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCbio", showCor = "character")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCbio", showCor = "missing")}: ... }
	 }
}

\examples{
showClass("LINCbio")
}
\keyword{classes}
