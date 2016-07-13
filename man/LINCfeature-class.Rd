\name{LINCfeature-class}
\Rdversion{1.1}
\docType{class}
\alias{LINCfeature-class}

\title{Class \code{"LINCfeature"}}
\description{
\code{"LINCfeature"}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("LINCfeature", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{customID}:}{Object of class \code{"character"} ~~ }
    \item{\code{customCol}:}{Object of class \code{"character"} ~~ }
    \item{\code{setLevel}:}{Object of class \code{"character"} ~~ }
    \item{\code{showLevels}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "LINCbio", e2 = "LINCfeature")}: ... }
    \item{+}{\code{signature(e1 = "LINCcluster", e2 = "LINCfeature")}: ... }
    \item{+}{\code{signature(e1 = "LINCmatrix", e2 = "LINCfeature")}: ... }
	 }
}

\examples{
showClass("LINCfeature")
}
\keyword{classes}
