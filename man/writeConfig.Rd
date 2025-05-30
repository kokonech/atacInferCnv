% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writeConfig.R
\name{writeConfig}
\alias{writeConfig}
\title{Write InferCNV configuration}
\usage{
writeConfig(resDir, sId, ctrlGrp = "Normal", binSize = NULL, meta = F)
}
\arguments{
\item{resDir}{Result directory path}

\item{sId}{Result name}

\item{ctrlGrp}{Name for the reference control cell type}

\item{binSize}{Size of the bin e.g. 100000 for 100Kbp}

\item{meta}{True if use meta cells, default false}
}
\description{
Write InferCNV configuration
}
