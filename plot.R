skeweta = function(eta,limit=.1)
{
	ind = which(eta < limit)
	if (length(ind) > 0) eta[ind] = seq(min(eta),max(eta[ind]),length.out=length(ind))
	eta
}

scale10 = function(x,scientific=TRUE)
{
	xmax = max(abs(x),na.rm=TRUE)
	if (length(xmax) == 0 || xmax == 0) return(1)

	xlog = log10(xmax/5)

	if (scientific) {
		10^(xlog%/%3*3)
	} else {
		10^round(xlog-.5)
	}
}

plotb = function(x,y,range=3,lty=1,col="bisque",pch="+",ylim=NULL,...)
{
	if (is.null(ylim)) ylim = range(boxplot(y,plot=FALSE)$stats)
	plot(x,colMeans(y),type="l",ylim=ylim,...)
	w = diff(par("usr")[1:2])/(1.5*length(x)+2)
	boxplot(y,range=range,col=col,lty=lty,pch=pch,pars=list(boxwex=w),add=TRUE,at=x,
		xaxt="n",yaxt="n")
}

plotbv = function(x,y,range=3,lty=1,col="bisque",pch="+",xlim=NULL,...)
{
	if (is.null(xlim)) xlim = range(boxplot(x,range=FALSE,plot=FALSE)$stats)
	plot(colMeans(x),y,type="l",xlim=xlim,...)
	w = diff(par("usr")[3:4])/(1.5*length(y)+2)
	boxplot(x,range=range,col=col,lty=lty,pch=pch,pars=list(boxwex=w),horizontal=TRUE,
		add=TRUE,at=y,xaxt="n",yaxt="n")
}

plotv = function(x,y,z,breaks,ylim=rev(range(y,finite=TRUE)),xaxs="i",yaxs="i",
	palette="YlOrRd",long=FALSE,...)
{
	if (par("mar")[1] > 2.5 || par("mar")[4] < 2) stop("mar does not fit")

	if (missing(breaks) || length(breaks) < 2) {
		breaks = prettyBreaks(z,crop=TRUE)$breaks
		dbr = diff(range(breaks))*.Machine$double.eps
		stopifnot(all(min(breaks)-dbr <= z & z <= max(breaks)+dbr))
	}

	if (length(breaks) == 2) breaks = c(breaks,breaks[2]+diff(breaks))

	rev = regexpr("\\+$",palette) < 0
	cols = hcl.colors(length(breaks)-1,sub("\\+$","",palette),rev=rev)
	if (long) {
		i = which(diff(x) < 0)
		if  (length(i) > 0) x[-(1:i)] = x[-(1:i)]+360
	}

	image(x,y,z,ylim=ylim,col=cols,breaks=breaks,xaxs=xaxs,yaxs=yaxs,...)
	lev = sprintf("% .3g",breaks)
	maplegend(lev,col=cols)
}

elems = function(n=1,np=1,nrow=1,ncol=1)
{
	# n elements to plot in np figures each on an [nrow,ncol] layout
	stopifnot(length(n) == 1)
	stopifnot(n > 0)
	stopifnot(length(ncol) == 1 && length(nrow) == 1)
	stopifnot(nrow > 0 && ncol > 0)

	# nelem: nb of elements per page
	nelem = nrow*ncol/np
	stopifnot(nelem == (nrow*ncol)%/%np)

	ni = (n-1)%/%nelem+1
	e = list()
	for (i in seq(ni)) e[[i]] = 1:min(n-nelem*(i-1),nelem)+nelem*(i-1)

	e
}
