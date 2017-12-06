\name{wavEquivFilter}
\alias{wavEquivFilter}
\title{Daubechies equivalent wavelet and scaling filters}
\usage{wavEquivFilter(wavelet="s8", level=1, scaling=FALSE, normalized=TRUE)}

\description{Calculates the equivalent wavelet and scaling coefficients 
for a given filter type at the specified wavelet transform level.}

\arguments{

\item{wavelet}{a character string denoting the filter type. Supported types include:
\describe{
\item{EXTREMAL PHASE (daublet):}{"haar", "d2", "d4", "d6", "d8", "d10", "d12", "d14", "d16", "d18", "d20"}

\item{LEAST ASYMMETRIC (symmlet):}{"s2","s4", "s6", "s8", "s10", "s12", "s14", "s16", "s18", "s20"}

\item{BEST LOCALIZED:}{"l2","l4", "l6", "l14", "l18", "l20"}

\item{COIFLET:}{"c6", "c12", "c18", "c24", "c30"}}

Default: \code{"s8"}.}

\item{level}{wavelet transform decomposition level}

\item{scaling}{a logical value. If \code{TRUE}, the equivalent scaling filter
at the specified level is output. Otherwise, the equivalent wavelet filter at the
specified level is output.
}

\item{normalized}{a logical value. If \code{TRUE},
the filters are normalized by dividing each filter
coefficient by the \eqn{\sqrt{2}}{sqrt(2)}
(useful for maximum overlap wavelet transforms). Default: \code{TRUE}.}
}

\value{
a vector of numeric values defining the equivalent filter at the specified level.
}

\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.

I. Daubechies,
\emph{Orthonormal Bases of Compactly Supported Wavelets},
Communications on Pure and, Applied Mathematics, 41, 909--96.
}

\seealso{
\code{\link{wavDaubechies}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavMODWPT}}.}

\examples{
wavStemPlot(wavEquivFilter(wavelet = 'd4', level = 2))
wavStemPlot(wavEquivFilter(wavelet = 'd4', level = 5, scaling = TRUE))
}
