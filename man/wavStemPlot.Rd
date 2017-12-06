\name{wavStemPlot}
\alias{wavStemPlot}

\title{
Stem plot
}

\description{
Plot of a discrete time series with a stem extending from the abscissa and with
a circle at the head of the stem.
}

\usage{
wavStemPlot(x, y, pch=16, stem='black', head='red', ...)
}

\arguments{
  \item{x}{abscissa values}
  \item{y}{ordinate values}
  \item{pch}{Head point type}
  \item{stem}{stem color}
  \item{head}{head color}
  \item{...}{additional plot arguments}
}

\value{NULL, returned invisibly}

\examples{
wavStemPlot(1:5, stem = 'red', head = 'blue')
}