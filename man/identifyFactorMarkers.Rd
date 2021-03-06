% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAI_model.R
\name{identifyFactorMarkers}
\alias{identifyFactorMarkers}
\title{Identify enriched features in each factor}
\usage{
identifyFactorMarkers(object, assay, features = NULL, cutoff.W = 0.5,
  cutoff.H = 0.5, thresh.pc = 0.05, thresh.fc = 0.25,
  thresh.p = 0.05, n.top = 10)
}
\arguments{
\item{object}{scAI object}

\item{assay}{Name of assay to be analyzed}

\item{features}{a vector of features}

\item{cutoff.W}{Threshold of feature loading values}

\item{cutoff.H}{Threshold of cell loading values}

\item{thresh.pc}{Threshold of the percent of cells enriched in one factor}

\item{thresh.fc}{Threshold of Fold Change}

\item{thresh.p}{Threshold of p-values}

\item{n.top}{Number of top features to be returned}
}
\description{
Rank features in each factor by Factor loading analysis
}
