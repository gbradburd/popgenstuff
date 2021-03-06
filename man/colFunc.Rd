% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viz.R
\name{colFunc}
\alias{colFunc}
\title{Function for making standardized color range}
\usage{
colFunc(x, cols, nCols, valRange)
}
\arguments{
\item{x}{A vector of numeric values for which
you want to generate colors that indicate
value.}

\item{cols}{A palette or vector of colors that you want
to plot in}

\item{nCols}{The total number of colors you want in the absolute
(standardized) scale}

\item{valRange}{The range of values for which you want to
generate colors. If unspecified, the function will
use the range of \code{x}. Default is \code{NULL}.}
}
\value{
A vector of colors that relate the value of each
element of \code{x} within the range specified by
\code{valRange}.
}
\description{
Function for making standardized color range
}
