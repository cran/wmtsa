######################################################
## S+Fractal Holder spectrum functionality
##
##   holderSpectrum
##
## Class: wtmmTree
## Constructor function: wtmmTree
## Methods:
##
##   [.wtmmTree
##   plot.wtmmTree
##
######################################################

###
# holderSpectrum
###

"holderSpectrum" <- function(chains, n.scale.min=3, fit=lmsreg)
{
  n.chain         <- length(chains)
  holder.time     <- vector(mode="integer", length=n.chain)
  holder.exponent <- vector(mode="numeric", length=n.chain)

  count <- 0

  for (i in seq(chains)){

    times    <- chains[[i]]$time
    logscale <- logb(chains[[i]]$scale, base=2)
    logwtmm  <- logb(chains[[i]]$wtmm, base=2)

    if (length(logscale) > 1){

      # order data from small to large scale
      if (logscale[2] < logscale[1]){
        times    <- rev(times)
        logscale <- rev(logscale)
        logwtmm  <- rev(logwtmm)
      }

      # obtain the appropriate scaling range
      # choose the range of scale which corresponds
      # to the smallest scales
      breaks <- linearSegmentation(logscale, logwtmm)
      div    <- ifelse1(length(breaks), breaks[1], length(times))

      if (div >= n.scale.min){

        count <- count + 1
        holder.time[count] <- times[1]
        datalist <- list(logwtmm=logwtmm, logscale=logscale, div=div)

        model <- ifelse1(div > 5, fit(logwtmm ~ logscale, subset=seq(div), data=datalist),
          lm(logwtmm ~ logscale, subset=seq(div), data=datalist))

        holder.exponent[count] <- coef(model)["logscale"]
      }
    }
  }

  holder.exponent <- holder.exponent[1:count]
  holder.time     <- holder.time[1:count]

  itime <- order(holder.time)
  list(exponent=holder.exponent[itime], time=holder.time[itime])
}

###
# wtmmTree
###

"wtmmTree" <- function(x, bridge.gaps=FALSE, strength.min= 0.0,
  n.octave.min=2, tolerance=0.0, border=FALSE)
{
  # define local functions
  "WTMM" <- function(x, tolerance=NULL){

    if (!is(x,"wavCWT"))
      stop("Input object must be of class wavCWT")

    # obtain attributes
    x.attr   <- attributes(x)
    times    <- x.attr$time
    scales   <- x.attr$scale
    n.sample <- x.attr$n.sample
    series   <- x.attr$series

    # check tolerances
    if (is.null(tolerance)){

      # use Donoho and Johnstone universal thresholding
      # tol=sqrt(2 * sigma^2 * log(N)) where sigma^2
      # is the variance of the noise, and N is the length
      # of the original time series. Since we do not
      # know sigma^2 a priori, we use the median absolute deviation
      # the level 1 DWT coefficients using the Haar filter as an
      # approxiation. Finally, we divide by the sqrt(scale) to make
      # the threshold scale dependent

      # noise.variance <- (median(abs((diff(series)[seq(2,n.sample-1,by=2)] / sqrt(2))))/ 0.6745)^2
      # vprint(noise.variance)
      # tolerance <- sqrt(2 * noise.variance * log(n.sample)) / sqrt(scales)
      tolerance <- mad(Mod(x[,1])) / scales
    }

    if (length(tolerance) < length(scales))
      tolerance <- tolerance[1] / sqrt(scales)

    wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima",
      as.matrix(x), tolerance,
      CLASSES=c("matrix","numeric"),
      COPY=rep(FALSE,2),
      PACKAGE="ifultools")

    itime  <- as.vector(wtmmz[[1]]) + 1
    iscale <- as.vector(wtmmz[[2]]) + 1

    return(list(time=times[itime],
      scale =scales[iscale],
      itime =itime,
      iscale=iscale,
      tolerance=tolerance))
  }

  # obtain attributes
  x.attr   <- attributes(x)
  times    <- x.attr$time
  scales   <- x.attr$scale
  n.sample <- x.attr$n.sample
  sampling.interval <- x.attr$sampling.interval
  border.times <- range(times) + sampling.interval * c(1,-1)

  # locate the the WTMM in the CWT matrix
  wtmm <- WTMM(x, tolerance=tolerance)

  # form the WTMM tree
  tree <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima_tree",
    as.integer(wtmm$itime - 1), as.integer(wtmm$iscale - 1),
    as.matrix(x), as.numeric(times), as.numeric(scales),
    as.logical(bridge.gaps), as.numeric(n.octave.min), as.numeric(strength.min),
    COPY=rep(FALSE,8),
    CLASSES=c("integer", "integer", "matrix", "numeric",
      "numeric", "logical", "numeric", "numeric"),
    PACKAGE="ifultools")

  z        <- list()
  ibranch  <- 1
  finetime <- NULL

  # prune border branches
  for (i in seq(tree)){

    branch <- tree[[i]]
    n.link <- ncol(branch)
    itime  <- as.integer(branch[1,]) + 1

    if (!border){

      # if the difference between successive times is (say)
      # greater than 1/4 the range of times, then the branch
      # is circularly wrapped in time. As the data is stored
      # from long to short scales, remove the longer scale data
      # of the wrapped branch by removing the first portion of
      # the corresponding vectors
      ibad  <- which(abs(diff(itime)) > n.sample/4)
      istart <- ifelse1(length(ibad) > 0, min(ibad[1] + 1, n.link), 1)
    }
    else
      istart <- 1

    igood <- seq(istart, n.link)

    z[[ibranch]] <- list(itime=itime[igood],
      iscale=as.integer(branch[2,igood]) + 1,
      time=branch[3,igood],
      scale=branch[4,igood],
      wtmm=branch[5,igood])

    finetime[ibranch] <- branch[3, n.link]
    ibranch <- ibranch + 1
  }

  # sort branch list with respect to time
  endtime <- finetime[1:(ibranch-1)]
  isort   <- order(endtime)

  z <- z[isort]

  names(z) <- seq(z)
  attr(z, "endtime") <- endtime[isort]
  attr(z, "time")  <- times
  attr(z, "scale") <- scales
  attr(z, "wtmm")  <- wtmm

  oldClass(z) <- "wtmmTree"

  z
}

###
# plot.wtmmTree
###

"plot.wtmmTree" <- function(x, add=FALSE, pch="o",  label=TRUE, log.="y",
  xlab=NULL, ylab=NULL, wtmm=FALSE, zoom=NULL,
  fit=FALSE, models= c("lm", "lmsreg", "ltsreg"),
  cex=0.8, col.skip=max(1,round(256/length(x))),  ...)
{

  branches <- names(x)

  if (!length(x))
    stop("No branches found in input object")

  if (fit){

    n.branch <- min(length(x), 4)

    nrow <- 1
    ncol <- 1

    if (n.branch > 1)
      nrow <- 2
    if (n.branch > 2)
      ncol <- 2

    frame()
    gap <- 0.11
    old.plt <- splitplot(nrow,ncol,1,gap=gap)
    on.exit(par(old.plt))

    lwd <- c(2, 1, 2, rep(1,10))
    lty <- c(1, 1, 4, seq(5,length=10))

    # initialize arguments
    slope <- vector(mode="numeric", length=length(models))

    for (i in seq(n.branch)){

      if (i > 1)
        splitplot(nrow,ncol,i,gap=gap)

      x.branch <- rev(log(x[[i]]$scale))
      y.branch <- rev(log(x[[i]]$wtmm))
      breaks   <- linearSegmentation(x.branch, y.branch)

      if (is.null(breaks))
        breaks <- length(x.branch)

      plot(x.branch, y.branch, xlab="log(scale)", ylab="log(|WTMM|)", type="b",
        pch="o", cex=cex)

      abline(v=x.branch[breaks], lty=2, xpd=FALSE)

      cut      <- seq(breaks[1]-1)
      x.branch <- x.branch[cut]
      y.branch <- y.branch[cut]

      # the lmsreg() and ltsreg() models will return an
      # error if less than 6 points are fit, so restrict
      # the models used accodingly
      if (length(x.branch) > 5)
        admissible.models <- models
      else
        admissible.models <- "lm"

      datalist <- list(x.branch=x.branch, y.branch=y.branch)

      for (imodel in seq(admissible.models)){

        eval(parse(text=paste("fit <- coef(", admissible.models[imodel],
          "(y.branch ~ x.branch,data=datalist))", sep="")))
        abline(fit, lty=lty[imodel], lwd=lwd[imodel], xpd=FALSE)
        slope[imodel] <- fit[2]
      }

      imodel <- seq(admissible.models)

      key.cex <- 0.7

      if (is.R()){
	mdlName <- upperCase(admissible.models)
	   legend("bottomright",
           paste(format(mdlName), "\t: slope=", format(round(slope[imodel], 4))),
    	     lty=lty[imodel],
    	     lwd=lwd[imodel],
    	     cex=key.cex)

			} else {

	      key(corner=c(1,0),
	        text=list(upperCase(admissible.models),adj=1, cex=key.cex),
	        line=list(lty=lty[imodel], lwd=lwd[imodel]),
	        text=list(paste("Slope :: ", as.character(round(slope[imodel], 4))),
	        adj=0, cex=key.cex))
			}
      mtext(paste("Branch", branches[i]), adj=1, cex=0.85, line=0.5)
    }

    return(NULL)
  }

  if (!add)
    frame()

  x.attr <- attributes(x)
  x.wtmm <- x.attr$wtmm

  if (wtmm){
    data <- scaleZoom(x.wtmm$time, x.wtmm$scale, zoom=zoom, log=log., xy.linked=TRUE,
      xlab="Time", ylab="Scale")
    if (!add)
      plot(data$x, data$y, col=6, pch="o", xlab=data$xlab, ylab=data$ylab, ...)
    else
      points(data$x, data$y, col=6, pch="o", ...)

    return(NULL)
  }

  if (!add){

    times  <- range(x.attr$time)
    scales <- range(x.attr$scale)

    if(is.element(log., "y")){
      scales <- logb(scales, base=2)
      if (is.null(ylab))
        ylab <- "log2(scale)"
    }
    else if (is.null(ylab))
      ylab <- "Scale"

    if(is.element(log., "x")){
      times <- logb(times, base=2)
      if (is.null(xlab))
        xlab <- "log2(time)"
    }
    else if (is.null(xlab))
      xlab <- "Time"

    plot(times, scales, type="n", xlab=xlab, ylab=ylab)
  }

  for (i in seq(along=x)){

    times  <- x[[i]]$time
    scales <- x[[i]]$scale

    if(is.element(log., "y"))
      scales <- logb(scales, base=2)
    if(is.element(log., "x"))
      times <- logb(times, base=2)

    if (label){

      imaxscale <- order(scales)[length(scales)]
      lines(times[-imaxscale], scales[-imaxscale], col=i*col.skip, pch=pch,
        cex=0.5, type="b", ...)
      text(times[imaxscale], scales[imaxscale], as.character(branches[i]),
        cex=1, col=i*col.skip)
    }
    else
      lines(times, scales, col=i, lty=1, pch=pch, ...)
  }

  invisible(NULL)
}

###
# [.wtmmTree
###

"[.wtmmTree" <- function(x, i, ..., time=NULL, range=NULL)
{
  ax    <- attributes(x)
  times <- ax$endtime

  if (!missing(range)){

    if (length(range) == 2){

      i <- which(times >= range[1] & times <= range[2])
    }
    else{
      i <- seq(length(x))
    }
  }
  else if(!missing(time)){

    i <- sort(time)

    min.scale <- median(diff(times))

    itime <- NULL

    for (j in i){

      itime <- c(itime, which(times >= j - min.scale & times <= j + min.scale))
    }

    i <- itime
  }

  z <- oldUnclass(x)[i]

  attributes(z) <- c(attributes(z),
    ax[ setdiff(names(ax), c("names", "dim", "dimnames")) ])

  z
}

