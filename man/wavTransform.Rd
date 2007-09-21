%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavTransform, [.wavTransform, [<-.wavTransform, [[.wavTransform, as.matrix.wavTransform, boxplot.wavTransform, eda.plot.wavTransform, plot.wavTransform, plot.wavTransform.crystal, print.wavTransform, print.summary.wavTransform, wavStackPlot.wavTransform, summary.wavTransform, reconstruct.wavTransform functions

\name{wavTransform}
\alias{wavTransform}
\alias{[.wavTransform}
\alias{[<-.wavTransform}
\alias{[[.wavTransform}
\alias{as.matrix.wavTransform}
\alias{boxplot.wavTransform}
\alias{eda.plot.wavTransform}
\alias{plot.wavTransform}
\alias{plot.wavTransform.crystal}
\alias{print.wavTransform}
\alias{print.summary.wavTransform}
\alias{wavStackPlot.wavTransform}
\alias{summary.wavTransform}
\alias{reconstruct.wavTransform}
\title{Constructor function for objects of class wavTransform}
\concept{class constructorwavelet}
\usage{wavTransform(data, series, n.levels, dictionary, shifted, xform)}
\description{Packs input information regarding a discrete wavelet transform
into a list.}

\arguments{
\item{data}{a list of vectors containing discrete wavelet
transform coefficients.}

\item{series}{a numeric vector or \code{signalSeries} object
representing the input series.}

\item{n.levels}{an integer denoting the number of decomposition levels.}

\item{dictionary}{an object of class \code{wavDictionary} representing
the wavelet dictionary of the transform.}

\item{shifted}{a logical value. If \code{TRUE}, it signifies that the
transform coefficients have already been shifted for approximate zero
phase alignment.}

\item{xform}{a character string denoting the type of wavelet
transform that has been performed. Typical values are
\code{"modwt"} or \code{"dwt"}.}

}

\value{
an object of class \code{wavTransform}.
}

\section{S3 METHODS}{
\describe{

\item{[}{single level data access.

Usage: x["d2"] or x[2]

Access a subset of wavelet transform crystals.}

\item{[<-}{single level data replacement method.

Usage: x["d2"] <- 1:4

Replace an entire crystal with explicitly defined coefficients.}

\item{[[}{double level data access.

Usage: x[["d2"]] or x[[2]]

Returns a vector of transform coefficients corresponding to the
specified crystal.}

\item{as.matrix}{transforms the list of wavelet transform coefficients
into a single-column matrix whose row names identify the transform coefficient,
e.g., d4(3) is the third coefficient of the \code{d4} crystal (fourth level
wavelet coefficients).

Usage: as.matrix(x)}

\item{boxplot}{plots a boxplot for each crystal in the discrete wavelet transform.

Usage: boxplot(x)}

\item{eda.plot}{extended data analysis plot. A 2x2 grid of plots containing
a stack plot, boxplot, and two energy related plots are shown.

Usage: eda.plot(x)}

\item{plot}{plot a discrete wavelet transform.
Usage: plot(x, type="h", plot.bar=TRUE, plot.pie=TRUE, add=FALSE, cex.main=ifelse1(is.R(),1,0.7), ...)

\describe{
\item{x}{A \code{wavTransform} object.}

\item{type}{A character denoting the type of line
plot to produce in a stack plot of the wavelet transform coefficients
(see primary plot() function for details). If type is the character
string \code{"energy"} an energy plot (bar or pie chart) is produced.
Default: \code{"h"}.}

\item{plot.bar}{A logical value. If \code{TRUE} and \code{type="energy"},
a bar plot of crystal energy is plotted. Default: \code{TRUE}.}

\item{plot.pie}{A logical value. If \code{TRUE} and \code{type="energy"},
a pie chart of crystal energy is plotted. Default: \code{FALSE}.}

\item{add}{A logical value. If \code{TRUE}, the plot is added to the
current plot layout without a frame ejection. Default: \code{FALSE}.}

\item{...}{Additional arguments to be sent to the stackplot routine.}}}

\item{print}{print the wavelet transform object.
Usage: print(x)}

\item{reconstruct}{reconstruct/synthesize/invert the wavelet transform.
Usage: reconstruct(x)

If the transform coefficients were not modified, the original
time series will be returned (+/- some numerical noise).}

\item{summary}{provide a statistical summary of the wavelet transform object.
Usage: z <- summary(x); print(z)}

\item{wavStackPlot}{stack plot of the wavelet transform.
Usage: wavStackPlot(x)}

}
}

\details{
%
Used internally by the \code{wavMODWT} and \code{wavDWT} functions
to package the transform contents into a list.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavDWPT}}, \code{\link{wavMODWPT}}, \code{\link{wavBoundary}}, \code{\link{wavSortCrystals}}, \code{\link{wavPacketIndices}}, \code{\link{wavShrink}}.}

\examples{
## calculate a MODWT of the sunpots series and 
## verify the class 
W <- wavMODWT(sunspots)
print(class(W))

## summarize the object 
summary(W)

## reconstruct the MODWT of the sunspots series 
## and compare to the original 
sunup <- reconstruct(W)
vecnorm(sunup - sunspots)
}
\keyword{univar}

