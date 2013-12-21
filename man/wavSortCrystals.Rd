%%  WARNING: This file was automatically generated from the associated 
%%  wav_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_util.mid file in the project MID directory. Once the 
%%  wav_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavSortCrystals function

\name{wavSortCrystals}
\alias{wavSortCrystals}
\title{Sort wavelet transform crystal names}
\concept{wavelet}
\usage{wavSortCrystals(x, reverse=FALSE)}
\description{Sorts the crystal names for a discrete wavelet transform by level
then by corresponding frequency content, from low to high.}

\arguments{
\item{x}{a vector of character strings containing the names of the wavelet transform crystals.}
\item{reverse}{a logical value. If \code{TRUE}, the order of the sorted names is reversed. Default: \code{FALSE}.}
}

\value{
a vector of character strings containing the sorted crystal names.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
W <- wavMODWT(1:100)
wavSortCrystals(names(W$data), reverse=TRUE)
}
\keyword{univar}

