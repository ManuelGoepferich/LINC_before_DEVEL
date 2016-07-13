\name{LINCsingle-class}
\Rdversion{1.1}
\docType{class}
\alias{LINCsingle-class}

\title{Class \code{"LINCsingle"}}
\description{
\code{"LINCsingle"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LINCsingle", ...)}.
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
    \item{plotlinc}{\code{signature(input = "LINCsingle", showCor = "character")}: ... }
    \item{plotlinc}{\code{signature(input = "LINCsingle", showCor = "missing")}: ... }
	 }
}

\examples{
showClass("LINCsingle")
}
\keyword{classes}
