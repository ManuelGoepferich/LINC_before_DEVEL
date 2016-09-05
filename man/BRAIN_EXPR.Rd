\name{BRAIN_EXPR}
\alias{BRAIN_EXPR}
\alias{cerebellum}
\alias{cortex}
\alias{crbl_cc}
\alias{crbl_cluster}
\alias{crbl_matrix}
\alias{ctx_cc}
\alias{ctx_cluster}
\alias{gbm_cluster}
\alias{glioblastoma}
\alias{hippocampus}
\alias{hpc_cluster}
\alias{pcgenes_crbl}
\alias{pcgenes_ctx}
\alias{pcgenes_gbm}
\alias{pcgenes_hpc}

\docType{data}
\title{
mRNA Expression Of Normal BRAIN From GTEX And TCGA
}
\description{
These datasets represent sample mRNA expression matrices of brain tissue with normal tissue from the Genotype-Tissue Expression (GTEx) platform (ctx = Cortex, crbl = Cerebellum, phc = Hippocampus) and cancer tissue from The Cancer Genome Atlas (TCGA) platform (gbm = Glioblastoma).
}
\usage{data(BRAIN_EXPR)}
\format{
Gene expression matrices and objects of the \code{LINC} class.
}

\value{
Rows represent genes and columns samples. Gene names are given as Entrez identifiers. For the normal tissue from GTEx, expression levels are in units of [FPKM], for cancer tissue from TCGA, the unit [RESM] is used to represent expression levels.
}

\source{
http://www.gtexportal.org
Genotype-Tissue Expression (GTEx)

https://tcga-data.nci.nih.gov/docs
The Cancer Genome Atlas (TCGA)

}
\references{
Carithers et al.
Biopreservation and Biobanking. October 2015, 13(5): 311-319. doi:10.1089/bio.2015.0032.
PMID: 26484571.
}
\examples{
data(BRAIN_EXPR)
}
\keyword{datasets}


