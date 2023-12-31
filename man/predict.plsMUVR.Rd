% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.plsMUVR}
\alias{predict.plsMUVR}
\title{Predict pls}
\usage{
\method{predict}{plsMUVR}(object, newdata, onlyPred = FALSE, scale = TRUE, ...)
}
\arguments{
\item{object}{a plsMUVR object}

\item{newdata}{new data}

\item{onlyPred}{Boolean for whether to report back predictions only (defaults to FALSE)}

\item{scale}{}

\item{...}{}
}
\value{
pls prediction
}
\description{
Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
}
