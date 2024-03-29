% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_pangenomes.R
\name{compare_pangenomes}
\alias{compare_pangenomes}
\title{compare_pangenomes}
\usage{
compare_pangenomes(
  fitA,
  fitB,
  family = "Tweedie",
  intercept = FALSE,
  modeldisp = TRUE,
  ci_type = "perc",
  conf = 0.95,
  nboot = 100,
  boot_pvalue = FALSE
)
}
\arguments{
\item{fitA}{a `panfit` object generated by the `panstripe` function}

\item{fitB}{a `panfit` object generated by the `panstripe` function}

\item{family}{the family used by glm. One of 'Tweedie', 'Poisson', 'Gamma' or 'Gaussian'. (default='Tweedie')}

\item{intercept}{whether or not to include an intercept term in the GLM (default=FALSE). Adding an intercept can increase the robustness of the algorithm to errors at higher branches of the phylogeny at the expense of less sensitivity.}

\item{modeldisp}{whether or not to model the dispersion as a function of the covariates of interest if using a Tweedie family (default=TRUE)}

\item{ci_type}{the method used to calculate the bootstrap CI (default='perc'). See \link[boot]{boot.ci} for more details.}

\item{conf}{A scalar indicating the confidence level of the required intervals (default=0.95)}

\item{nboot}{the number of bootstrap replicates to perform (default=100)}

\item{boot_pvalue}{whether or not to calculate bootstrap p-values (default=FALSE)}
}
\value{
a list containing a summary of the comparison and the resulting `glm` model object
}
\description{
Compares the coefficients between two pangenomes datasets by including interaction terms in the GLM.
}
\examples{

simA <- simulate_pan(rate=1e-4, ngenomes = 200, mean_trans_size=3, fp_error_rate=1)
simB <- simulate_pan(rate=1e-3, ngenomes = 200, mean_trans_size=3, fp_error_rate=1)
fitA <- panstripe(simA$pa, simA$tree, nboot=0)
fitA$summary
fitB <- panstripe(simB$pa, simB$tree, nboot=0)
fitB$summary
comp <- compare_pangenomes(fitA, fitB, nboot=10)
comp$summary

}
