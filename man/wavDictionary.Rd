%%  WARNING: This file was automatically generated from the associated 
%%  wav_dict.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_dict.mid file in the project MID directory. Once the 
%%  wav_dict.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavDictionary, as.list.wavDictionary, print.wavDictionary functions

\name{wavDictionary}
\alias{wavDictionary}
\alias{as.list.wavDictionary}
\alias{print.wavDictionary}
\title{Constructor function for objects of class wavDictionary}
\concept{class constructorwavelet}
\usage{wavDictionary(wavelet, dual, decimate, n.sample,
    attr.x, n.levels, boundary, conv, filters,
    fast, is.complex)}
\description{Packs input information regarding a discrete wavelet transform
into a dictionary list.}

\arguments{
\item{wavelet}{a character string denoting the type of wavelet used in the transform.}

\item{dual}{a logical value. If \code{TRUE}, it signifies that a dual transform
was performed.}

\item{decimate}{a logical value. If \code{TRUE}, it signifies that a decimated transform
was performed.}

\item{n.sample}{an integer representing the number of samples in
the original time series.}

\item{attr.x}{a list of additional (arbitrary) attributes to append onto
the output object.}

\item{n.levels}{an integer denoting the number of decomposition levels.}

\item{boundary}{a character string denoting the boundary extension type
used in transform. Supported values are \code{"zero"}, \code{"periodic"}, \code{"reflection"},
and \code{"continue"}.}

\item{conv}{a logical value. If \code{TRUE}, it signifies that a convolution style transform
was performed (as opposed to correlation style).}

\item{filters}{a list of vectors named \code{"scaling"} and \code{"wavelet"}
containing the scaling and wavelet filter coefficients, respectively.}

\item{fast}{a logical value. If \code{TRUE}, it signifies that a fast pyramidal scheme
was  used to develop the decimated transform as opposed to
calculating the transform coefficients via an
explicit matrix multiplication of the wavelet transform matrix and
the original time series.}

\item{is.complex}{a logical value. If \code{TRUE}, it signifies the transform
was complex-valued.}

}

\value{
an object of class \code{wavDictionary}.
}

\section{S3 METHODS}{
\describe{

\item{print}{print the dictionary.

Usage: print(x)}

}
}

\details{
%
Used internally by the \code{wavMODWT} and \code{wavDWT} functions
to package the transform contents into a dictionary list.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
## create a faux wavelet dictionary 
wavelet <- "s8"
wavDictionary(wavelet=wavelet, dual=FALSE,
    decimate=FALSE, n.sample=1024,
    attr.x=NULL, n.levels=3,
    boundary="periodic", conv=TRUE,
    filters=wavDaubechies(wavelet),
    fast=TRUE, is.complex=FALSE)
}
\keyword{univar}

