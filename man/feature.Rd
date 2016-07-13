\name{feature}
\alias{feature}

\title{
 Manipulate Objects Of The 'LINC' class
}
\description{
\code{feature} provides useful options intended to be used with \code{'LINC'} objects.
}
\usage{
feature(setLevel   = NULL,
        customID   = NULL,
        customCol  = "black",
        showLevels = FALSE)
}

\arguments{
  \item{setLevel}{
a \code{character} string of the class the object should be converted into
}
  \item{customID}{
a \code{character} string of a name for the object
}
  \item{customCol}{
a \code{character} string of a valid colour for the object 
}
  \item{showLevels}{
whether to show the inherited classes of the object \code{TRUE} or not \code{FALSE}
}
}

\details{
Custom ids and colours enable the identification of a particular object in plots created by \code{plotlinc} and \code{querycluster}. With \code{setLevel} it is possible to change the class of an object.

\code{feature} works in combination with the plus operator: \code{'object' + feature(customID = ..., customCol = ...)}

}

\value{
an object of the class 'LINCfeature' (S4) with 5 \code{Slots} (not shown)
}

\author{
Manuel Goepferich
}

\seealso{
\code{querycluster} ; 
\code{strlinc}
}
\examples{
data(BRAIN_EXPR)

# add a custom name
crbl_cluster_feat <- crbl_cluster + feature(customID = "CEREBELLUM")
plotlinc(crbl_cluster_feat)
}  
\keyword{ ~feature }
