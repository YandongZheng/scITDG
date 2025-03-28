% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateITDGObject.R
\name{CreateITDGObject}
\alias{CreateITDGObject}
\title{Title}
\usage{
CreateITDGObject(
  object,
  celltype.use = NULL,
  celltype = NULL,
  time.points.use = NULL,
  time.points = NULL,
  random.value = TRUE,
  seed = 176,
  sample.ncell = 200
)
}
\arguments{
\item{object}{An seurat's object.}

\item{celltype.use}{A metadata column name to store the cell type. Currently only supported for class-level (i.e. non-quantitative) attributes.}

\item{celltype}{Cell type used for analysis.}

\item{time.points.use}{A metadata column name to store the time point. Currently only supported for class-level (i.e. non-quantitative) attributes.}

\item{time.points}{Variables of time point. Must be in the same order as the time point in the studied.}

\item{random.value}{Use random sampling. Default is true.}

\item{seed}{The number of seeds used for random sampling.}

\item{sample.ncell}{The number of random samples.}
}
\description{
Title
}
