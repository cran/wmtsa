%%  WARNING: This file was automatically generated from the associated 
%%  wav_mrd.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_mrd.mid file in the project MID directory. Once the 
%%  wav_mrd.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavMRD, [.wavMRD, [<-.wavMRD, [[.wavMRD, as.matrix.wavMRD, boxplot.wavMRD, crystal.names.wavMRD, dotchart.wavMRD, eda.plot.wavMRD, plot.wavMRD, print.wavMRD, print.summary.wavMRD, reconstruct.wavMRD, wavStackPlot.wavMRD, summary.wavMRD functions

\name{wavMRD}
\alias{wavMRD}
\alias{[.wavMRD}
\alias{[<-.wavMRD}
\alias{[[.wavMRD}
\alias{as.matrix.wavMRD}
\alias{boxplot.wavMRD}
\alias{crystal.names.wavMRD}
\alias{dotchart.wavMRD}
\alias{eda.plot.wavMRD}
\alias{plot.wavMRD}
\alias{print.wavMRD}
\alias{print.summary.wavMRD}
\alias{reconstruct.wavMRD}
\alias{wavStackPlot.wavMRD}
\alias{summary.wavMRD}
\title{Calculate the detail sequences for wavelet transform crystals}
\concept{class constructordetail function generationmultiresolution decompositionwavelet}
\usage{wavMRD(x, level=NULL, osc=NULL)}
\description{Let \eqn{\mathbf{W}_{j,n}}{W(j,n)} be a discrete
wavelet packet crystal where \eqn{j} is the
decomposition level and \eqn{n} is the oscillation
index. The detail sequence \eqn{{\mathcal{D}}_{j,n}}{D(j,t)}
is formed (essentially) by reconstructing the transform after zeroing out all other crystals
except \eqn{\mathbf{W}_{j,n}}{W(j,n)}.
The \code{wavMRD} function
calculates the details for a DWT and MODWT in an optimized way.}

\arguments{
\item{x}{an object of class \code{wavTransform}.}
\item{level}{an integer (vector) containing the decomposition level(s) corresponding to the
crystal(s) to be decomposed. Default: If the input is of
class \code{wavTransform}, then the default is to
return the details at all levels of the transform, i.e., a full multiresolution decomposition.}

\item{osc}{an integer (vector) containing the oscillation indices corresponding to the crystal(s)
to be decomposed. Default: the default values are coordinated with that of the \code{level} argument.}
}

\value{
an object of class \code{WaveletMRD}.
}

\section{S3 METHODS}{
\describe{

\item{[}{single level data access.

Usage: x["D2"] or x["S4"]

Access a subset of wavelet transform details/smooth.}

\item{[<-}{single level data replacement method.

Usage: x["D2"] <- 1:4

Replace an entire wavelet transform details/smooth with explicitly defined coefficients.}

\item{[[}{double level data access.

Usage: x[["D2"]] or x[[2]]

Returns a vector of wavelet transform detail/smooth coefficients corresponding to the
specified crystal.}

\item{as.matrix}{transforms the list of wavelet transform details/smooth coefficients
into a single-column matrix whose row names identify the transform coefficient,
e.g., D4(3) is the third coefficient of the \code{D4} detail.

Usage: as.matrix(x)}

\item{boxplot}{plots a boxplot for each wavelet transform detail/smooth.

Usage: boxplot(x)}

\item{crystal.names}{return the crystal names for each wavelet transform detail/smooth.

Usage: crystal.names(x)}

\item{plot}{plot a stack plot of the discrete wavelet transform details/smooth.
Usage: plot(x, n.top=15, vgap=.05, col=1, show.sum=TRUE, add=FALSE, ...)

\describe{
\item{x}{A \code{wavMRD} object.}

\item{n.top}{An integer defining the (maximum) number of
top-most energetic crystals to plot. Default: \code{15}.}

\item{sort.energy}{A logical value. If \code{TRUE}, the crystals
are sorted in the display from the most energetic (top) to the
least energetic (bottom) of the specified \code{n.top} crystals.
Default: \code{FALSE}.}

\item{vgap}{A numeric scalar defining the vertical gap between plots
expressed as a fraction of the maximum value of the details/smooth
that are plotted. Default: \code{0.05}.}

\item{col}{An integer or vecto rof integers deining the color index to apply
to each detail/smooth line plot. Default: \code{1}.}

\item{show.sum}{A logical value. If \code{TRUE}, a plot of the sum
of all details/smooth is also plotted. Default: \code{TRUE}.}

\item{add}{A logical value. If \code{TRUE}, the plot is added to the
current plot layout without a frame ejection. Default: \code{FALSE}.}

\item{...}{Additional arguments to be sent to the plot routine.}}}

\item{print}{print the wavelet transform details/smooth object.
Usage: print(x)}

\item{reconstruct}{reconstruct/synthesize/invert the wavelet transform details/smooth.
Usage: reconstruct(x)

If the transform coefficients were not modified, the original
time series will be returned (+/- some numerical noise).}

\item{summary}{provide a statistical summary of the wavelet transform details/smooth object.
Usage: z <- summary(x); print(z)}

\item{wavStackPlot}{stack plot of the wavelet transform details/smooth.
Usage: wavStackPlot(x)}

}
}

\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavMRDSum}}, \code{\link{reconstruct}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
## calculate various wavelet transforms of the 
## first difference of a linear chirp sequence 
x <- make.signal("linchirp", n=1024)
x.dwt <- wavDWT(x, n.levels = 3)
x.modwt <- wavMODWT(x, n.levels = 3)

## calculate the wavelet details for all crystals 
## of the DWT and MODWT 
wavMRD(x.dwt)
wavMRD(x.modwt)

## plot the wavelet details for levels 1 and 3 of 
## the MODWT 
plot(wavMRD(x.modwt, level = c(1,3)))

## plot the wavelet details for all levels of the 
## MODWT of a linear chirp. 
plot(wavMRD(x.modwt))
}
\keyword{univar}

