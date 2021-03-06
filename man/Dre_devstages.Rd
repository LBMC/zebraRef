% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dre_devstages-data.R
\docType{data}
\name{Dre_devstages}
\alias{Dre_devstages}
\title{Dre_devstages}
\format{
An object of class \code{list} of length 2.
}
\usage{
Dre_devstages
}
\description{
Table of reference datasets and landmark developmental stages.
}
\details{
The timing of developmental stages were compiled from the reference dataset \insertCite{domazet2010phylogenetically}{zebraRef}.
}
\section{Format }{

A list with 2 dataframes:
\itemize{
\item \verb{$devstages}
\item \verb{$datasets}
}

Both dataframes have the following fields:
\itemize{
\item \code{name}
\item \code{tstart}
\item \code{tend}
\item \code{tunit}
}
}

\references{
\insertAllCited{}
}
\keyword{datasets}
