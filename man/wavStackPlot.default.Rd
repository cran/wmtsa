%%  WARNING: This file was automatically generated from the associated 
%%  wav_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_plot.mid file in the project MID directory. Once the 
%%  wav_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavStackPlot.default function

\name{wavStackPlot.default}
\alias{wavStackPlot.default}
\title{Plot wavelet transform crystals}
\concept{wavelettransforms}
\usage{\method{wavStackPlot}{default}(x, x.axis=TRUE, y.axis=TRUE, type="l", plot=TRUE,
    bars=FALSE, vgap=.05, grid=FALSE, times=time(x[[1]]),
    grid.lty="dashed", same.scale=NULL,
    zerocenter=FALSE, zeroline=FALSE, col=rep(1,n),
    complex.math="mod", cex.main=0.7, cex.axis=0.7, ...)}
\description{Forms a stack plot of wavelet transform crystals.}

\arguments{
\item{x}{a named list of wavelet transform crystals.}
\item{...}{additional arguments sent directly to various internal plot functions.}

\item{bars}{a logical value. If \code{TRUE}, a vertical bars are produced. Default: \code{FALSE}.}

\item{cex.axis}{axis label character expansion factor ala \code{par} function. Default: 0.7.}

\item{cex.main}{main title character expansion factor ala \code{par} function. Default: 0.7.}

\item{col}{color indices for each crystal ala \code{par} function. Default: \code{rep(1,n)}.}

\item{complex.math}{math function to perform on each crystal prior to display. Default: \code{"mod"}.}

\item{grid}{a logical value. If \code{TRUE}, grid lines are produced. Default: \code{FALSE}.}

\item{grid.lty}{grid line type ala \code{par} function. Default: \code{"dashed"}.}

\item{plot}{a logical value. If \code{TRUE}, a plot is produced. Default: \code{TRUE}.}

\item{same.scale}{a logical value. If \code{TRUE}, the crystal coefficients are scaled. Default: \code{TRUE}.}

\item{times}{a numeric vector of x-axis times to be displayed. Default: \code{time(x[[1]])}.}

\item{type}{type of plot ala \code{par} function. Default: \code{"l"}.}

\item{vgap}{vertical gap factor used to separate stacked crystals. Default: 0.05.}

\item{x.axis}{a logical value defining the visibility of the x-axis. Default: \code{TRUE}.}

\item{y.axis}{a logical value defining the visibility of the y-axis. Default: \code{TRUE}.}

\item{zerocenter}{a logical value. If \code{TRUE}, crystals are centered about zero. Default: \code{FALSE}.}

\item{zeroline}{a logical value. If \code{TRUE}, a y=0 line is drawn for reference. Default: \code{FALSE}.}
}
\details{
%
Used by various wavelet functions to plot wavelet transform crystals. In general,
the user should not call this function directly and, rather, shoudl rely on the plot methods
of \code{wavTransform} objects.
}
\seealso{
\code{\link{wavTransform}}.}

\examples{
wavStackPlot.default(wavMODWT(sunspots)$data)
}
\keyword{hplot}

