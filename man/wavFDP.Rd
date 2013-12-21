%%  WARNING: This file was automatically generated from the associated 
%%  wav_fdp.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_fdp.mid file in the project MID directory. Once the 
%%  wav_fdp.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavFDP, eda.plot.wavFDP, plot.wavFDP, print.wavFDP, print.summary.wavFDP, summary.wavFDP functions

\name{wavFDP}
\alias{wavFDP}
\alias{eda.plot.wavFDP}
\alias{plot.wavFDP}
\alias{print.wavFDP}
\alias{print.summary.wavFDP}
\alias{summary.wavFDP}
\title{Class constructor for block- and time-dependent wavelet-based FD model parameter estimators}
\concept{class constructor}
\usage{wavFDP(estimator, delta, variance.delta,
    innovations.variance, delta.range, dictionary, levels,
    edof.mode, boundary, series, sdf.method, type)}
\description{Class constructor for block- and time-dependent wavelet-based FD model parameter estimators.}

\arguments{
\item{estimator}{character string briefly describing the estimator.}

\item{delta}{numeric value/vector denoting the estimated FD model parameter.}

\item{innovations.variance}{numeric value/vector denoting the estimated FD innovations variance.}

\item{variance.delta}{numeric value/vector defining the variance of delta.}

\item{delta.range}{two element numeric vector defining the range of delta.}

\item{dictionary}{wavelet transform dictionary used in the analysis.}

\item{levels}{vector of integers denoting the wavelet decomposition levels used in the analysis.}

\item{edof.mode}{an integer on [1,3] defining the equivalent degrees of freedom mode used in the analysis.}

\item{boundary}{a list containing named objects \code{mode} and \code{description}, containing
a logical value and a character string, respectively. The \code{mode} object should be
be \code{TRUE} if a boundary treatment was used, and \code{description} should contain
a description of the boundary treatment.}

\item{series}{a signSeries object containing the input series.}

\item{sdf.method}{a character string defining the SDF method used in
the analysis, e.g., \code{"Integration lookup table"}.}

\item{type}{a character string defining the type of estimator,
e.g., \code{""instantaneous""} or \code{"block"}.}

}

\section{S3 METHODS}{
\describe{

\item{eda.plot}{extended data analyis plot of the data. Available options are:
\describe{
\item{mean.delta}{mean value of delta, plotted as a horizontal reference
line for instantaneous delat estimations. Default: \code{NULL} (no reference line).}

\item{xlab}{character string defining the x-axis label. Default: \code{"Time"}.}

\item{ylab}{character string defining the y-axis label. Default: biased/unbiased and
estimator condition.}

\item{title.str}{character string defining the main title of the plot. Default: \code{NULL} (no title).}

\item{type}{character string defining type of plot ala \code{par} function. Default: \code{"l"} (solid line).}}}

\item{plot}{plots a summary of the results. Available options are:
\describe{
\item{mean.delta}{mean value of delta, plotted as a horizontal reference
line for instantaneous delat estimations. Default: \code{NULL} (no reference line).}

\item{xlab}{character string defining the x-axis label. Default: \code{"Time"}.}

\item{ylab}{character string defining the y-axis label. Default: biased/unbiased and
estimator condition.}

\item{title.str}{character string defining the main title of the plot. Default: \code{NULL} (no title).}

\item{type}{character string defining type of plot ala \code{par} function. Default: \code{"l"} (solid line).}

\item{show.key}{a logical value. If \code{TRUE}, a key of the plot is shown. Default: \code{TRUE}.}

\item{conf.color}{color index ala \code{par} for the confidence intervals. Default: \code{16}.}}}

\item{print}{prints the object. Available options are:
\describe{
\item{digits}{number of digits to use in displaying numeric values. Default: \code{5}.}}}

\item{print.summary}{print a summary of the results.}

\item{summary}{create a summary of the results.}

}
}

\seealso{
\code{\link{wavFDPBlock}}, \code{\link{wavFDPTime}}.}

\examples{
## create a faux dictionary 
dictionary <- wavDictionary(wavelet="s8",
    dual=FALSE, decimate=FALSE, n.sample=512,
    attr.x=NULL, n.levels=5,
    boundary="periodic", conv=TRUE,
    filters=wavDaubechies("s8"),
    fast=TRUE, is.complex=FALSE)

## construct a faux wavFDP object 
z <- wavFDP(estimator="wlse",
    delta=0.45,
    variance.delta=1.0,
    innovations.variance=1.0,
    delta.range=c(-10.0,10.0),
    dictionary=dictionary,
    levels=c(1,3:4),
    edof.mode=2,
    boundary=list(mode=TRUE,description="unbiased"),
    series=create.signalSeries(fdp045),
    sdf.method="Integration lookup table",
    type="block")

## print the result 
print(z)
}
\keyword{models}

