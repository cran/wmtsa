%%  WARNING: This file was automatically generated from the associated 
%%  wav_sig.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_sig.mid file in the project MID directory. Once the 
%%  wav_sig.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the create.signalSeries function

\name{create.signalSeries}
\alias{create.signalSeries}
\title{Converts various time series to an object of class}
\concept{utilities}
\usage{create.signalSeries(x=NULL, position=list(from=1,by=1,units=character()),
    from=NULL, by=NULL, to=NULL, length.out=NULL,
    units=character(), title.data=character(), documentation=character(), na.rm=TRUE)}
\description{Converts numeric data to an object of class \code{} containing
one dimensional data. The input data is assumed to be uniformly sampled.}

\arguments{

\item{by}{a numeric containing the sampling rate at
which the values in \code{data} should be extracted. This
parameter must coordinate with the \code{position} arguments and can be
used in combination with the \code{by}, \code{to}, or
\code{length.out} arguments. This argument is not the same as the
\code{position}\$\code{by} argument which denotes the sampling
rate of the original \code{data}. Default: \code{NULL}.}

\item{documentation}{a string used to describe the input
\code{data}. Default: \code{character()}.}

\item{from}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{NULL}.}

\item{length.out}{an integer containing the maximum number
of values to extract from \code{data}. Because \code{data} is a
finite length sequence, the actual number of values
returned may be less than that specified by this argument depending
upon the conditions imposed by the \code{from} and \code{by} arguments.
The \code{length.out} argument should not be specified if both the
\code{from} and \code{to} arguments are specified.
Default: \code{NULL}.}

\item{na.rm}{a logical flag used to indicate if NaN values should be removed from the input. Default: \code{TRUE}.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{title.data}{a string representing the name of the input
\code{data}. Default: \code{character()}.}

\item{to}{a numeric containing the end point
in \code{data} from which the values should be extracted. This
parameter must coordinate with the \code{position} arguments and can be
used in combination with the \code{by} and \code{to} arguments. The
\code{length.out} argument should not be specified if both the
\code{from} and \code{to} arguments are specified. Default: \code{NULL}.}

\item{units}{a string denoting the units of the time series. Default: \code{character()}.}

\item{x}{a numeric vector, matrix or an object of class \code{ts}
(uniform sampling assumed). If a matrix, the longest first row or column is
extracted. Default: \code{NULL} (no data).}
}

\value{
an object of class \code{signalSeries}.
}
\seealso{
\code{\link{make.signal}}.}

\examples{
## convert an explicitly developed numeric vector 
x <- 1:10
create.signalSeries(x)

## now impose hypothetical positions on the data 
create.signalSeries(x, pos=list(from=0.3, by=0.1))

## extract the values from position 0.5 onward 
create.signalSeries(x, pos=list(from=0.3, by=0.1), from=0.5)

## extract the values from position 0.5 onward, 
## but keep only the first 3 values of the 
## extraction 
create.signalSeries(x, pos=list(from=0.3, by=0.1), from=0.5, length=3)

## extract the values from position 0.5 onward and 
## skip every other point (sample the data at 
## 0.2 position intervals) 
create.signalSeries(x, pos=list(from=0.3, by=0.1), from=0.5, by=0.2)

## simply return the first 4 values, and supply a 
## title and some documentation comments to the 
## data 
create.signalSeries(x, length=4, title="Faux Data", doc="An example")
}
\keyword{ts}

