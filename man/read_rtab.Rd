% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_rtab.R
\name{read_rtab}
\alias{read_rtab}
\title{read_rtab}
\usage{
read_rtab(file)
}
\arguments{
\item{file}{a binary tab delimited gene presence/absence matrix (.Rtab) with gene names in the first column}
}
\value{
a binary matrix
}
\description{
Loads an Rtab pangenome gene presence/absence file as produced by Panaroo, Roary and other prokaryotic pangenome programs
}
\examples{

file <- system.file("extdata", "gene_presence_absence.Rtab", package = "panstripe")
pa <- read_rtab(file)

}
