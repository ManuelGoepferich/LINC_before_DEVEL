\name{overlaylinc-methods}
\docType{methods}
\alias{overlaylinc-methods}
\alias{overlaylinc}
\alias{overlaylinc,LINCbio,LINCbio-method}
\title{
 Compare Two Lists Of Biological Terms In A Plot
}
\description{
\code{overlaylinc} takes two \code{LINCbio} objects as input and computes the intersection of the first with the second. This can be useful in case the two \code{LINCbio} objects represent different statistical processing, tissues, controls, thresholds or batches.
}
\section{Methods}{
\describe{

\item{\code{signature(input1 = "LINCbio", input2 = "LINCbio")}}{
(see details)
}
}}

\usage{
overlaylinc(input1,
            input2)
}

\arguments{
  \item{input1}{
 the first \code{LINCbio} object 
}
  \item{input2}{
the second \code{LINCbio} object
}
}

\details{
The first object will be taken as reference. This function is an alternative to \code{plotlinc} which can be applied on each input object separately or on the intersection of both objects (see example).
}

\value{
an object of the class 'gtable' containing multiple 'grobs'
}

\author{
Manuel Goepferich
}

\seealso{
\code{feature} ; 
\code{plotlinc}
}
\examples{
data(BRAIN_EXPR)

# compare two brain tissues using 'overlaylinc'
ctx_cc <- ctx_cc + feature(customID = "Cortex vs. Cerebellum")
overlaylinc(ctx_cc, crbl_cc)

# alternatively, compute the intersection of both objects
brain_cc <- ctx_cc + crbl_cc
plotlinc(brain_cc)
}  
\keyword{ ~overlaylinc }
\keyword{methods}
