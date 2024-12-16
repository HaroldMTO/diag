plotb = function(x,y,lty=1,col="bisque",pch="+",ylim=NULL,...)
{
	if (is.null(ylim)) ylim = range(boxplot(y,outline=FALSE,plot=FALSE)$stats)
	plot(x,colMeans(y),type="l",ylim=ylim,...)
	w = diff(par("usr")[1:2])/(1.5*length(x)+2)
	boxplot(y,col=col,lty=lty,pch=pch,pars=list(boxwex=w),add=TRUE,at=x,xaxt="n",yaxt="n")
}

plotbv = function(x,y,lty=1,col="bisque",pch="+",xlim=NULL,...)
{
	if (is.null(xlim)) xlim = range(boxplot(x,outline=FALSE,plot=FALSE)$stats)
	plot(colMeans(x),y,type="l",xlim=xlim,...)
	w = diff(par("usr")[3:4])/(1.5*length(y)+2)
	boxplot(x,col=col,lty=lty,pch=pch,pars=list(boxwex=w),horizontal=TRUE,add=TRUE,at=y,
		xaxt="n",yaxt="n")
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

