% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panfit.R
\name{new_panfit}
\alias{new_panfit}
\title{new_panfit}
\usage{
new_panfit(summary, model, data, boot_samples, tree, pa)
}
\arguments{
\item{summary}{data.frame or tibble}

\item{model}{glm}

\item{data}{data.frame or tibble}

\item{boot_samples}{data.frame or tibble}

\item{tree}{phylo}

\item{pa}{matrix, data.frame or tibble}
}
\value{
A panfit object
}
\description{
class for storing output of main panstripe fitting function
}
