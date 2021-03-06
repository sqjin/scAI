% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scAI_model.R
\name{getAggregatedData}
\alias{getAggregatedData}
\title{Generate the aggregated epigenomic data}
\usage{
getAggregatedData(object, group = NULL)
}
\arguments{
\item{object}{an scAI object after running run_scAI}

\item{group}{cell group information if available; aggregate epigenomic data based on the available cell group information instead of the learned cell-cell similarity matrix from scAI}
}
\description{
Generate the aggregated epigenomic data
}
