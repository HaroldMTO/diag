Gdiag = "~/util/diag"

library(maps)
library(mapproj)
library(parallel)
library(mffield)

Gvp0 = 101325

readDom = function(filename,domd=list())
{
	df = read.table(filename,header=TRUE)

	ind = which(df$west == "-" | df$east == "-" | df$south == "-" | df$north == "-")
	if (length(ind) == 0) {
		domd = list()
	} else {
		indd = match(df$name[ind],names(domd))
		if (any(is.na(indd))) stop(paste("unknown domains in",filename))

		domd = domd[indd]
		df = df[-ind,]
		if (dim(df)[1] == 0) return(domd)

		df[c("west","east","south","north")] = lapply(df[c("west","east","south","north")],
			as.numeric)
	}

	stopifnot(all(df$east > df$west & df$north > df$south))
	stopifnot(all(-90 < df$north & df$north <= 90))
	stopifnot(all(-90 <= df$south & df$south < 90))

	doms = list()

	for (i in seq(dim(df)[1])) {
		doms[[i]] = new("Domain",xlim=c(df$west[i],df$east[i]),
			ylim=c(df$south[i],df$north[i]),proj=df$proj[i])
		if (df$param[i] != "-") doms[[i]]@param = df$param[i]
	}

	names(doms) = df$name
	c(domd,doms)
}

getGrid = function(con)
{
	nx = readBin(con,"integer",1,endian="swap")
	type = readBin(con,"integer",nx/4,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	# int is 4 bytes, num is nfp bytes (= 4 or 8)
	nfp = type[1]

	# 0: global (Gauss) grid, 1: LAM grid
	if (type[2] == 0) {
		nx = readBin(con,"integer",1,endian="swap")
		gem = readBin(con,"numeric",nx/nfp,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		ndgnh = nx/4
		nloeng = readBin(con,"integer",ndgnh,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		ndglg = 2*ndgnh
		nloeng[(ndgnh+1):ndglg] = rev(nloeng)

		grid = list(nfp=nfp,pole=gem[1:2],stretch=gem[3],nlat=ndglg,nlong=nloeng,
			npdg=sum(nloeng),gauss=TRUE)
	} else if (type[2] == 1) {
		ndglg = type[3]
		ndlon = type[4]
		grid = list(nfp=nfp,nlat=ndglg,nlong=rep(ndlon,ndglg),npdg=ndglg*ndlon,gauss=FALSE)
	} else {
		stop("unknown grid type")
	}

	grid
}

getVCoord = function(con)
{
	nx = readBin(con,"integer",1,endian="swap")
	nl1 = nx/nfp
	Ah = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))
	nx = readBin(con,"integer",1,endian="swap")
	Bh = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	eta = ((Ah[-1]+Ah[-nl1])/Gvp0+(Bh[-1]+Bh[-nl1]))/2
	#nlevel = nl1-1

	eta
}

getTime = function(con)
{
	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/4 == 2)
	base = readBin(con,"integer",2,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/nfp == 1)
	step = readBin(con,"numeric",1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	base = as.POSIXct(as.Date(as.character(base[1]),"%Y%m%d"))+base[2]
	new("FcTime",base=base,step=as.integer(step))
}

is.equal = function(grid,g4)
{
	if (! identical(grid$nlong,g4@nlong)) return(FALSE)

	if (grid$gauss != is(g4,"GaussGrid")) return(FALSE)

	if (grid$gauss) {
		if (! isTRUE(all.equal(grid$pole,g4@pole))) return(FALSE)
		if (! isTRUE(all.equal(grid$stretch,g4@stretch))) return(FALSE)
	}
}

getVars = function(con,grid,vars,quiet=FALSE)
{
	lvar = vector("list",length(vars))

	cvars = sub("^\\.?","",vars)

	while (TRUE) {
		nx = readBin(con,"integer",1,endian="swap")
		if (length(nx) == 0 || nx == 0) break

		varname = trimws(readChar(con,nx))
		if (! quiet) cat("reading variable",varname,nx,"\n")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		stopifnot(nx == 4)
		nl = readBin(con,"integer",1,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))

		nx = readBin(con,"integer",1,endian="swap")
		npdg = sum(grid$nlong)
		stopifnot(nx/grid$nfp == nl*npdg)

		i = match(sub("^\\.?","",varname),cvars)
		if (is.na(i)) {
			seek(con,nl*npdg*grid$nfp+4,"current")
			next
		}

		lvar[[i]] = readBin(con,"numeric",nl*npdg,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))
		dim(lvar[[i]]) = c(npdg,nl)
		names(lvar)[i] = varname
		if (all(! sapply(lvar,is.null))) break
	}

	lvar
}

setDiagData = function(grid,eta,ilev,data)
{
	if (dim(data)[2] == 1) {
		f = new("Field",grid=grid,eta=1)
		f = setDatPart(f,data)
	} else {
		nl = dim(data)[2]
		f = new("Field",grid=grid,eta=eta[ilev])
		if (nl == length(ilev)) {
			f = setDatPart(f,data)
		} else {
			if (nl != length(ilev)+1 && nl != length(eta)) {
				warning("inconsistent levels for bin data: ",nl," != ",length(eta))
			}

			cat("--> change 'bin' levels:",dim(data)[2],"-->",length(ilev),"\n")
			f = setDatPart(f,data[,ilev,drop=FALSE])
		}
	}

	f
}

getFrame = function(file)
{
	out = system(sprintf("epy_dump.py %s -f frame -o %s",file,Gficbin),intern=TRUE)
	out = grep("dimensions",out,ignore.case=TRUE,value=TRUE)
	cat(out[1],"\n")

	con = file(Gficbin,"rb")
	dims = readBin(con,"integer",5,size=8)

	# global (Gauss grid) or LAM (Cartesian grid)
	gauss = dims[5] == 0
	if (gauss) {
		nlong = readBin(con,"integer",dims[1],size=8)
		nwave = readBin(con,"integer",dims[2],size=8)
		ngp = dims[4]
		nl = dims[3]
	} else {
		ngp = prod(dims[1:2])
		nl = dims[5]
	}

	base = readBin(con,"numeric",1)
	base = as.POSIXct(base,origin="1970-01-01")
	step = readBin(con,"numeric",1)
	if (regexpr(".+\\+0*([0-9]+)",file) > 0) {
		ech = as.integer(gsub(".+\\+0*([0-9]+).*","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	fc = new("FcTime",base=base,step=as.integer(step))

	lats = readBin(con,"numeric",ngp)
	longs = readBin(con,"numeric",ngp)
	longs = (longs+180)%%360-180

	Ah = readBin(con,"numeric",nl+1)
	Bh = readBin(con,"numeric",nl+1)
	eta = ((Ah[-1]+Ah[-(nl+1)])/Gvp0+(Bh[-1]+Bh[-(nl+1)]))/2

	if (gauss) gem = readBin(con,"numeric",3)

	close(con)

	if (gauss) {
		nlat = length(nlong)
		theta = 90-180*(seq(nlat)-.5)/nlat
		g4 = new("GaussGrid",lat=lats,long=longs,nlong=nlong,pole=gem[2:3],stretch=gem[1],
			theta=theta)

		#frame = list(nlat=nlat,nwave=dims[2],nlevel=nl,npdg=ngp,nlong=nlong,theta=theta,
		#	lat=lats,long=longs,eta=eta,base=fc@base,step=fc@step,mucen=gem[2],locen=gem[3],
		#	rstret=gem[1],gauss=TRUE,lam=FALSE,g4=g4,fc=fc)
		frame = list(nwave=dims[2],nlevel=nl,npdg=ngp,eta=eta,g4=g4,fc=fc)
	} else {
		nlong = rep(dims[2],dims[1])
		g4 = new("LAMGrid",lat=lats,long=longs,nlong=nlong,center=rep(0,2))

		#frame = list(nlat=dims[1],nlong=nlong,nwavex=dims[3],nwavey=dims[4],nlevel=nl,
		#	npdg=ngp,lat=lats,long=longs,eta=eta,base=fc@base,step=fc@step,gauss=FALSE,
		#	lam=TRUE,g4=g4,fc=fc)
		frame = list(nwavex=dims[3],nwavey=dims[4],nlevel=nl,npdg=ngp,eta=eta,g4=g4,fc=fc)
	}

	frame
}

formatStep = function(step,hour="")
{
	if (step %% 3600 == 0) {
		sprintf("%d%s",step/3600,hour)
	} else if (step %% 60 == 0) {
		sprintf("%dm",step/60)
	} else {
		if (step != round(step)) warning("step is rounded to an int")
		sprintf("%ds",round(step))
	}
}

getField = function(fic,param,symbol,frame,frlow=frame,cache.alt=FALSE)
{
	stopifnot(! is.null(frlow$ilev))

	fc = frame$fc
	fsave = sprintf("%s/e%s/%s.RData",dirname(fic),formatStep(fc@step),symbol)
	if (cache.alt) fsave = sub("\\.RData",".2.RData",fsave)

	if (file.exists(fsave) && file.info(fsave)$mtime < file.info(fic)$mtime) {
		cat("--> cache file older than data file, removed\n")
		file.remove(fsave)
	}

	grid = frlow$g4
	npdg = sum(grid@nlong)

	lread = TRUE
	if (file.exists(fsave)) {
		ilev = 0
		fic.save = ""
		noms = load(fsave)
		stopifnot(all(c("data","ilev","step","fic.save") %in% noms))

		# read data if different
		lread = ! identical(ilev,frlow$ilev) || fic.save != fic || dim(data)[1] != npdg
		if (lread) {
			cat("--> cache file does not fit, change for 2nd cache file\n")
			fsave = sprintf("%s/e%s/%s.2.RData",dirname(fic),formatStep(fc@step),symbol)
			if (cache.alt) fsave = sub("\\.2\\.RData",".RData",fsave)
			if (file.exists(fsave) && file.info(fsave)$mtime < file.info(fic)$mtime) {
				cat("--> 2nd cache file older than data file, removed\n")
				file.remove(fsave)
			}
		}
	}

	if (file.exists(fsave) && lread) {
		load(fsave)

		# read data if different
		lread = ! identical(ilev,frlow$ilev) || fic.save != fic || dim(data)[1] != npdg

		if (fic.save != fic) {
			cat("--> different origin file",fic.save," read data again\n")
		} else if (! identical(ilev,frlow$ilev)) {
			cat("--> different levels, read data again\n")
		} else if (dim(data)[1] != npdg) {
			cat("--> different grid, read data again\n")
		} else {
			cat("--> cache file read:",fsave,"\n")
		}
	}

	if (lread) {
		grid = frame$g4
		data <- getGPFields(fic,param,sum(grid@nlong))
	}
	stopifnot(dim(data)[1] == length(grid))

	ilev = frlow$ilev
	if (dim(data)[2] == 1) {
		eta = 1
	} else if (dim(data)[2] == length(ilev)) {
		eta = frlow$eta[ilev]
	} else {
		# data should have full (or reduced) nlev (if not, we don't know its levels)
		eta = frame$eta
		if (! is.null(frame$ilev)) eta = eta[frame$ilev]
	}
	stopifnot(dim(data)[2] == length(eta))

	f = new("Field",grid=grid,eta=eta)
	f = setDataPart(f,data)
	if (! lread) return(f)

	h = hist(f,plot=FALSE)
	if (length(which(h$counts > 0)) == 2) {
		stopifnot(length(table(f)) == 2)
		meth = "ppp"
	} else {
		meth = "linear"
	}

	if (sum(grid@nlong) != npdg) {
		nlat = length(grid@nlong)
		cat("--> interpolation from",sum(grid@nlong),nlat,"gridpoints/lat - method:",
			meth,"\n")
		it = system.time(f <- interp(f,frlow$g4,meth))
		cat("interpolation time:",it,"\n")
	}

	if ("ind" %in% names(frlow)) {
		cat("--> selection of gridpoints:",length(frlow$g4),"/",npdg,"\n")
		f = setDataPart(f,f[frlow$ind,,drop=FALSE])
	}

	if (dim(f)[2] > 1 && length(f@eta) != length(ilev)) {
		cat("--> interpolation from",length(f@eta),"levels - method:",meth,"\n")
		f = interpAB(f,frlow$eta[ilev],meth)
	}

	data = getDataPart(f)
	step = fc@step
	fic.save = fic
	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	save("data","ilev","step","fic.save",file=fsave)

	f
}

getGPFields = function(file,field,npdg)
{
	cmd = sprintf("epy_dump.py %s -f '%s' -o %s",file,field,Gficbin)
	rt = system.time(out <- system(cmd,intern=TRUE))
	cat("FAreading time:",rt,"\n")
	out = grep("date|time|step|gridpoint size",out,ignore.case=TRUE,value=TRUE)
	cat(paste(out,collapse=", "),"\n")

	con = file(Gficbin,"rb")
	nl = readBin(con,"integer",1,size=8)
	stopifnot(length(nl) == 1)

	data = matrix(nrow=npdg,ncol=nl,dimnames=list(NULL,seq(nl)))
	rt = system.time(for (j in seq(nl)) data[,j] = readBin(con,"numeric",npdg))
	cat("Binreading time:",rt,"\n")
	stopifnot(nl == readBin(con,"integer",1,size=8))
	close(con)

	data
}

getSPFields = function(file,field,nwave)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",1,size=8)
	data = readBin(con,"numeric",nsp2)
	close(con)

	stopifnot(nsp2 == nwave*(nwave-1))

	matrix(data,nrow=nwave)
}

ucs = function(nord,x,y)
{
	x*nord$m-y*nord$l
}

vcs = function(nord,x,y)
{
	y*nord$m+x*nord$l
}

equalize = function(mapf,nc=5,offset=1.5)
{
	# nc equally numbered groups of latitudes (mapf: normalized grid-point density)
	mapc = quantile(mapf,prob=seq(nc-1)/nc)
	cat(nc-1,"quantiles of mapf:\n")
	print(mapc)

	# mean of mapf for each group
	im = findInterval(mapf,mapc)+1
	mapm = tapply(mapf,im,mean)

	# index of last element of each group
	indl = c(match(seq(2,nc),im)-1,length(mapf))

	n = as.integer((max(mapm)+offset)/(mapm+offset))

	ilat = integer()
	i1 = 1
	for (i in seq(nc)) {
		# equalizing density (offset to limit n)
		cat("--> lat group/from/to/by:",i,i1,indl[i],n[i],"- mean:",mapm[i],"\n")
		ilat = c(ilat,seq(i1,indl[i]+1,by=n[i]))
		i1 = indl[i]+1
	}

	if (max(ilat) > length(mapf)) ilat = ilat[-length(ilat)]
	unique(ilat)
}

dilat = function(grid)
{
	clats = c(0,cumsum(grid@nlong))
	nlat = length(grid@nlong)

	ddm = numeric(nlat)

	for (i in seq(nlat)) {
		ip = clats[i]+seq(grid@nlong[i])

		dlon = abs(diff(grid@long[ip[1:2]]))
		a = cos(grid@lat[ip[1]]*pi/180)
		if (grid@pole[1] < 1) {
			# check that 1st points are in zonal direction (ie lats vary less at start)
			stopifnot(abs(diff(grid@lat[ip[1:2]])) < abs(diff(grid@lat[ip[3:4]])))
		}

		# distances in zonal and meridian directions
		dlzon = a*pmin(dlon,360-dlon)
		dlmer = diff(grid@lat[ip[1:2]])

		# rough estimate of "true distance" between points (ie mapping factor, roughly)
		ddm[i] = sqrt(dlzon^2+dlmer^2)
	}

	# normalization to 1
	ddm/max(ddm[1:(nlat/2)])
}

prettyBreaks = function(x,breaks="Sturges",nmin=5,n=8,crop=FALSE,split=FALSE)
{
	h = hist(x,breaks,plot=FALSE)

	if (is.character(breaks)) {
		stopifnot(1 < nmin && nmin <= n)

		have = h$counts > 0
		nbin = length(h$counts)

		# case for binary values (e.g. 0/1): 2 non-empty bins that are the extreme ones
		if (length(which(have)) == 2 && all(which(have) %in% c(1,nbin))) {
			h = hist(x,2,plot=FALSE)
			cat("--> only 2 populated extreme bins, try to limit to 2:",length(h$counts),"\n")
		} else if (nbin > 1.25*n) {
			h = hist(x,n,plot=FALSE)
			if (length(h$counts) > 1.25*n) {
				cat("--> failed to limit bins from",length(which(have)),nbin,"to",n,":",
					length(h$counts),"\n")
			}

			# still too many bins: try n-1
			if (FALSE && length(h$counts) > 1.2*n && length(which(h$counts > 0)) > 2) {
				h = hist(x,n-1,plot=FALSE)
				cat("--> try to limit bins from",length(which(h$counts > 0)),
					length(h$counts),"to",n,":",length(h$counts),"\n")
			}
		}

		if (length(which(h$counts > 0)) > 2 && length(h$counts) < nmin) {
			hh = hist(x,nmin,plot=FALSE)
			cat("--> too few bins, try minimum of",nmin,"bins:",length(hh$counts),"\n")
			# yes, 0s can be fewer with more bins...
			if (length(which(hh$counts == 0)) <= length(which(h$counts == 0))) h = hh
		}
	}

	if (split && length(h$counts) > 2) {
		indx = order(h$density,decreasing=TRUE)
		dy = diff(h$density[indx])
		ix = which.min(dy)
		if (length(which(h$counts > 0)) > 1 && ix < length(which(h$density > 0)) &&
			(r=h$density[indx[ix]]/h$density[indx[ix+1]]) > 3*ix) {
			# n is 3, 4 or 5 (ie 2, 3 or 4)
			nb = 1+min(4,as.integer(sqrt(r/(3*ix))+1))
			br = seq(h$breaks[indx[1]],h$breaks[indx[1]+1],length.out=nb)
			if (nb > 3) cat("--> splitting bin",ix,"/",length(h$counts),"into",nb,"\n")
			h = hist(x,unique(sort(c(h$breaks,br))),plot=FALSE)
		}
	}

	if (crop && length(h$counts) > 2) {
		br = h$breaks
		dxn = diff(br[1:2])
		dxx = diff(br[-(1:(length(br)-2))])
		# crop the extreme bins by 1/5-th of their width or by their half width
		if (min(x) > br[2]-dxn/5 && max(x) < br[length(br)-1]+dxx/5) {
			h$breaks[1] = br[2]-dxn/5
			h$breaks[length(br)] = br[length(br)-1]+dxx/5
		} else if (min(x) > br[2]-dxn/2 && max(x) < br[length(br)-1]+dxx/2) {
			h$breaks[1] = br[2]-dxn/2
			h$breaks[length(br)] = br[length(br)-1]+dxx/2
		}
	}

	# values in last bin are stuck to left bound (if right=F, may be?)
	if (length(h$counts) > 1 && h$breaks[length(h$breaks)-1] == max(x,na.rm=TRUE)) {
		cat("--> get rid of the last bin",length(h$breaks),h$counts[length(h$counts)],"\n")
		h = hist(x,h$breaks[-length(h$breaks)])
	}

	h
}

mapxy = function(dom,...)
{
	if (dom@proj == "-") {
		l = map(xlim=dom@xlim,ylim=dom@ylim,...)
		p = .Last.projection()
		if (p$projection != "") {
			cat("--> reset projection\n")
			p = .Last.projection(list(projection="",parameters=NULL,orientation=NULL))
			stopifnot(p$projection == "")
		}

		# axes for rectangular projection only
		axis(1)
		axis(2)
	} else {
		l = map(projection=dom@proj,parameters=dom@param,xlim=dom@xlim,ylim=dom@ylim,...)
	}

	l
}

mapdom = function(dom,points,data,main=NULL,breaks="Sturges",palette="YlOrRd",pch=20,
	cex=.6,ppi=72,quiet=FALSE,...)
{
	if (par("mar")[1] > 2.5 || par("mar")[4] < 2.5) stop(paste(par("mar"),collapse=" "))
	if (ppi > 144) stop("ppi > 144")

	l = mapxy(dom,mar=par("mar"),new=TRUE)
	box()

	h = prettyBreaks(data,breaks,crop=TRUE)

	f = par("fin")
	npmax = as.integer(min(prod(f*ppi/(4*cex)),.Machine$integer.max))
	if (length(data) < npmax/100) {
		cat("--> very few points, magnify plotting symbol (x2)\n")
		cex = 2*cex
	} else if (length(data) > 1.2*npmax) {
		cex = max(.2,round(cex*sqrt(npmax/length(data)),3))
		npmax = as.integer(min(prod(f*ppi/(4*cex)),.Machine$integer.max))
	}

	if (length(data) > 1.2*npmax) {
		if (! quiet) {
			cat("--> reducing xy points from",length(data),"to",npmax,"and cex to",cex,"\n")
		}

		ind = select(points,npmax)

		b2 = cut(data[ind],h$breaks)

		ilost = which(h$counts > 0 & table(b2) == 0)
		if (length(ilost) > 0) {
			b = cut(data,h$breaks)
			ind1 = which(b %in% levels(b)[ilost])
			ind = c(ind,ind1)
			stopifnot(all(! duplicated(ind)))
			if (! quiet) {
				cat("--> selecting back",length(ind1),"lost points in",length(ilost),
					"data bins\n")
			}
		}

		points = points[ind]
		data = data[ind]
	}

	if (cex < .2) cex = .2

	br = h$breaks

	ind = findInterval(data,br,rightmost.closed=TRUE)
	rev = regexpr("\\+$",palette) < 0
	cols = hcl.colors(length(br),sub("\\+$","",palette),rev=rev)

	tind = table(ind)

	if (dom@proj == "-") {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			points(points@long[ii],points@lat[ii],col=cols[ind[ii]],pch=pch,cex=cex,...)
		}
	} else {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			mp = mapproject(points@long[ii],points@lat[ii])
			points(mp$x,mp$y,col=cols[ind],pch=pch,cex=cex,...)
		}
	}

	levels = sprintf("% .3g",br)
	if (any(duplicated(levels))) levels = sprintf("% .4g",br)
	DOPfill.legend(levels,col=cols)

	lines(l)
	title(main)
}

mapsegments = function(dom,lat,long,...)
{
	# S->N cross-section
	map.grid(c(long,long,dom@ylim[1],dom@ylim[2]),nx=1,labels=FALSE,pretty=FALSE,...)

	# W->E cross section
	map.grid(c(dom@xlim[1],dom@xlim[2],lat,lat),ny=1,labels=FALSE,pretty=FALSE,...)
}

mapdom2 = function(dom,points,zx,zy,main=NULL,breaks="Sturges",colvec=TRUE,length=.05,
	angle=15,ppi=6,quiet=FALSE,...)
{
	if (ppi > 96) stop("ppi > 96\n")

	l = mapxy(dom,mar=par("mar"),new=TRUE)
	box()

	f = par("fin")
	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (length(zx) > 1.2*npmax) {
		if (! quiet) cat("--> reducing xy2 plot from",length(zx),"to",npmax,"points\n")
		ind = select(points,npmax)

		points = points[ind]
		zx = zx[ind]
		zy = zy[ind]
	}

	ffz = sqrt(zx^2+zy^2)

	u = par("usr")
	ux = diff(u[1:2])
	uy = diff(u[3:4])

	# scale of wind vectors: a fraction of the smallest side
	asp = min(ux,uy)/max(ux,uy)
	ffu = asp*sqrt((ux^2+uy^2)/length(zx))

	if (colvec) {
		fx = zx/ffz*ffu
		fy = zy/ffz*ffu

		br = prettyBreaks(ffz,crop=TRUE)$breaks
		ic = findInterval(ffz,br,rightmost.closed=TRUE)
		palette = rainbow(length(br),start=2/3,end=0)
		cols = palette[ic]
	} else {
		ff1 = ffu*scale(ffz,center=FALSE)
		fx = zx/ffz*ff1
		fy = zy/ffz*ff1
		cols = rep(1,length(ff1))
	}

	x2 = points@long+fx
	y2 = points@lat+fy

	ffi = sqrt((fx/ux*f[1])^2+(fy/uy*f[2])^2)

	ind = ffi > 2.e-3

	if (dom@proj == "-") {
		arrows(points@long,points@lat,x2,y2,length,angle,col=cols[ind],...)
	} else {
		mp = mapproject(points@long[ind],points@lat[ind])
		mp2 = mapproject(x2[ind],y2[ind])
		arrows(mp$x,mp$y,mp2$x,mp2$y,length,angle,col=cols[ind],...)
	}

	if (any(! ind)) {
		if (dom@proj == "-") {
			points(points@long[! ind],points@lat[! ind],pch=1,col=1)
		} else {
			mp = mapproject(points@long[! ind],points@lat[! ind])
			points(mp$x,mp$y,pch=1,col=1)
		}
	}

	if (colvec) {
		levels = sprintf("% .3g",br)
		if (any(duplicated(levels))) levels = sprintf("% .4g",br)
		DOPfill.legend(levels,col=palette)
	}

	lines(l)
	title(main)
}

DOPfill.legend = function(breaks,col,...)
{
	u = par("usr")
	p = par("plt")

	width = (1-p[2])/6*diff(u[1:2])/diff(p[1:2])
	x = u[2]+width/3

	nl = length(breaks)
	height = diff(u[3:4])
	dy = height/(nl-1)
	y = u[3]
	ybas = y + dy*(seq(nl-1)-1)
	yhaut = ybas + dy

	rect(x,ybas,x+width,yhaut,col=col,border=NA,xpd=TRUE)

	op = par(las=2,yaxt="s")
	axis(4,c(ybas[1],yhaut),breaks,tick=FALSE,pos=x+width/12,mgp=c(1,.7,0),...)
	par(op)
}

hist.annot = function(data,nmin=10,n=20,split=TRUE,...)
{
	if (par("mar")[1] < 1.5 || par("mar")[4] > 3) stop(paste(mar,collapse=" "))

	h = prettyBreaks(data,nmin=nmin,n=n,split=split)
	h = hist(data,h$breaks,...)

	rug(range(h$breaks),ticksize=-.03,lwd=.8)
	if (h$equidist) {
		y = max(h$counts)
	} else {
		y = max(h$density)
	}

	znx = range(data)
	text(znx[1],-.02*y,sprintf("| %.3g",min(data,na.rm=TRUE)),adj=0,col=4)
	text(znx[2],-.02*y,sprintf("%.3g |",max(data,na.rm=TRUE)),adj=1,col=4)

	m = mean(data,na.rm=TRUE)

	s = sd(data,na.rm=TRUE)
	segments(m,.9*y,m,y,col=4)
	arrows(m-s,.95*y,m+s,.95*y,length=.03,angle=90,code=3,col=4)

	dx = diff(range(h$breaks))/15

	if (m+dx > h$breaks[1]+.9*diff(range(h$breaks))) {
		text(m-dx,.97*y,sprintf("m: %.3g",m),adj=1,col=4)
		text(m-dx,.92*y,sprintf("m/sd: %.3g",m/s),adj=1,col=4)
	} else {
		text(m+dx,.97*y,sprintf("m: %.3g",m),adj=0,...)
		text(m+dx,.92*y,sprintf("m/sd: %.3g",m/s),adj=0,col=4)
	}
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
	DOPfill.legend(lev,col=cols)
}

plotsectiongeo = function(field,dom,main,...)
{
	longm = mean(dom@xlim)
	latm = mean(dom@ylim)

	#br = prettyBreaks(field,crop=TRUE)$breaks

	yz = sectiongeo(field,longm,dom@ylim)
	stopifnot(all(is.finite(yz$data)))

	main[2] = sprintf("S->N section, longitude %g",round(longm,1))
	plotv(yz$lats,field@eta,yz$data,main=main,xlab="",ylab="",...)

	xz = sectiongeo(field,dom@xlim,latm)
	stopifnot(all(is.finite(xz$data)))

	main[2] = sprintf("W->E section, latitude %g",round(latm,1))
	plotv(xz$longs,field@eta,xz$data,main=main,xlab="",ylab="",long=TRUE,...)
}

mapexi = function(f,grmap,desc,doms,prefix=character(),main,np=4,prob=seq(0,100)/100,
	mnx=TRUE,force=FALSE,...)
{
	nl = length(f@eta)
	qv = array(dim=c(length(prob),nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))

	for (i in seq(along=doms)) {
		d4 = doms[[i]]
		ind = inDomain(f@grid,d4)
		if (length(which(ind)) < np) {
			cat("--> less than",np,"points in domain",names(doms)[i],"\n")
			next
		}

		if (is(f@grid,"GaussGrid") && diff(d4@xlim) > 350 && diff(d4@ylim) > 170) {
			# no domain, no zoom
			domin = d4
			xy = f
			if (! is.null(grmap)) {
				xy@grid = grmap
				xy = setDataPart(xy,f[grmap@ind,,drop=FALSE])
			}
		} else {
			pts = f@grid[ind]
			domin = new("Domain",long=pts@long,lat=pts@lat,proj=d4@proj,param=d4@param)
			if (area(domin,d4) < .1) {
				cat("--> area too small for domain",names(doms)[i],"\n")
				next
			}

			xy = zoom(f,d4)
		}

		stopifnot(all(is.finite(xy)))
		qv[,,i] = apply(xy,2,quantile,prob=prob)

		if (length(prefix) == 0) next

		dom = names(doms)[i]
		cat(".. domain:",dom,"- points/total:",length(xy@grid),dim(f)[1],"\n")

		ficpng = sprintf("%s/map%s%s_%s.png",pngd,prefix,dom,desc$symbol)
		if (file.exists(ficpng) && ! force) next

		png(ficpng)

		tt = main

		if (nl == 1) {
			op = par(c(Gparm,list(cex=.8)))
			mapdom(doms[[i]],xy@grid,xy[,1],main=tt,palette=desc$palette,...)
		} else {
			op = par(c(Gparm,list(mfrow=c(2,2),cex=.7)))

			# be careful of nl=2...
			indl = as.integer(seq(1,nl,length.out=min(4,nl)))
			if (length(Gindl) == 2) indl = c(1,Gindl,nl)
			for (il in indl) {
				tt[2] = sprintf("%s level %d (eta %g)",main[2],il,round(f@eta[il],2))
				mapdom(doms[[i]],xy@grid,xy[,il],main=tt,palette=desc$palette,
					quiet=il!=indl[1],...)
				mapsegments(domin,lat=mean(domin@ylim),long=mean(domin@xlim),col="darkgrey")
			}
		}

		dev.off()

		ficpng = sprintf("%s/hist%s%s_%s.png",pngd,prefix,dom,desc$symbol)

		png(ficpng)

		tt[2] = main[2]

		if (nl == 1) {
			op = par(c(Gparh,list(cex=.8)))
			hist.annot(xy,col="whitesmoke",main=tt,xlab=NULL,ylab=NULL)
		} else {
			op = par(c(Gparm,list(mfcol=c(2,2)),cex=.7))
			plotsectiongeo(f,d4,main=desc$longname,palette=desc$palette)

			par(mar=Gparh$mar)
			hist.annot(xy,col="whitesmoke",main=tt,xlab=NULL,ylab=NULL)
			plotbv(qv[,,i],xy@eta,main=tt[1],xlab="",ylab="",ylim=yeta)
			abline(v=0,col="darkgrey",lty=2)
		}

		dev.off()

		if (! mnx || nl < 4) next

		il = (nl+1)%/%2

		ficpng = sprintf("%s/mapn%s%s_%s.png",pngd,prefix,dom,desc$symbol)

		png(ficpng)
		op = par(mfrow=c(2,2),cex=.7)
		par(Gparm)
		datan = apply(xy[,1:il],1,min)
		tt[2] = sprintf("%s min. of lev. [1,%d]",main[2],il)
		mapdom(d4,xy@grid,datan,main=tt,palette="Heat+",...)
		par(mar=Gparh$mar)
		hist.annot(datan,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

		par(Gparm)
		datan = apply(xy[,-(1:il)],1,min)
		tt[2] = sprintf("%s min. of lev. [%d,%d]",main[2],il+1,nl)
		mapdom(d4,xy@grid,datan,main=tt,palette="Heat+",...)
		par(mar=Gparh$mar)
		stopifnot(all(par("mar") == Gparh$mar))
		hist.annot(datan,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
		dev.off()

		ficpng = sprintf("%s/mapx%s%s_%s.png",pngd,prefix,dom,desc$symbol)

		png(ficpng)
		op = par(mfrow=c(2,2),cex=.7)
		par(Gparm)
		datax = apply(xy[,1:il],1,max)
		tt[2] = sprintf("%s max. of lev. [1,%d]",main[2],il)
		mapdom(d4,xy@grid,datax,main=tt,palette="Heat",...)
		par(mar=Gparh$mar)
		hist.annot(datax,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

		par(Gparm)
		datax = apply(xy[,-(1:il)],1,max)
		tt[2] = sprintf("%s max. of lev. [%d,%d]",main[2],il+1,nl)
		mapdom(d4,xy@grid,datax,main=tt,palette="Heat",...)
		par(mar=Gparh$mar)
		hist.annot(datax,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
		dev.off()
	}

	qv
}

mapex2 = function(u,v,doms,prefix=NULL,main,np=4,palette,...)
{
	nl = length(u@eta)
	qv = array(dim=c(101,nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))

	for (i in seq(along=doms)) {
		d4 = doms[[i]]$d4
		ind = inDomain(u@grid,d4)
		if (length(which(ind)) < np) next

		if (is(u@grid,"GaussGrid") && diff(d4@xlim) > 350 && diff(d4@ylim) > 170) {
			# world: no domain, no zoom
			domin = d4
			x = u
			y = v
		} else {
			pts = u@grid[ind]
			domin = new("Domain",long=pts@long,lat=pts@lat,proj=d4@proj,param=d4@param)
			if (area(domin,d4) < .1) next

			x = zoom(u,d4)
			y = zoom(v,d4)
		}

		ffdom = sqrt(x^2+y^2)
		stopifnot(all(is.finite(ffdom)))
		qv[,,i] = apply(ffdom,2,quantile,prob=seq(0,100)/100)

		if (is.null(prefix)) next

		dom = names(doms)[i]
		tt = "u/v wind"
		cat(".. domain:",dom,"\n")

		ficpng = sprintf("%s/map%s%s_ff.png",pngd,prefix,dom)
		if (file.exists(ficpng)) next

		png(ficpng)

		if (nl == 1) {
			op = par(c(Gparm,list(cex=.8)))
			mapdom2(d4,x,x,y,main=main,...)
		} else {
			op = par(c(Gparm,list(mfcol=c(2,2),cex=.7)))

			# if nl=2...
			indl = as.integer(seq(1,nl,length.out=min(4,nl)))
			if (length(Gindl) == 2) indl = c(1,Gindl,nl)
			for (il in indl) {
				tt[2] = sprintf("%s level %d (eta %g)",main[2],il,round(u@eta[il],2))
				mapdom2(doms[[i]],x,x[,il],y[,il],main=tt,...)
				mapsegments(domin,lat=mean(domin@ylim),long=mean(domin@xlim),col="darkgrey")
			}
		}

		dev.off()

		ficpng = sprintf("%s/hist%s%s_ff.png",pngd,prefix,dom)

		png(ficpng)
		tt = "wind speed"
		if (nl == 1) {
			op = par(c(Gparh,list(cex=.8)))
			hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab="wind speed")
		} else {
			ff = sqrt(u^2+v^2)
			ff = setDataPart(u,ff)
			op = par(c(Gparm,list(mfcol=c(2,2),cex=.7)))
			plotsectiongeo(ff,doms[[i]],main=tt[1],palette=desc$palette)

			par(mar=Gparh$mar)
			hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
			plotbv(qv[,,i],u@eta,main=tt[1],xlab="",ylab="",ylim=yeta)
			abline(v=0,col="darkgrey",lty=2)
		}

		dev.off()
	}

	qv
}

mapgpole = function(grid,u,v,prefix)
{
	nord = compass(g4)
	#if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
	u = ucs(nord,u,v)
	v = vcs(nord,u,v)

	p = pole(grid)
	pole = new("Domain",xlim=p[2]+c(-6,6),ylim=p[1]+c(-5,5),proj="stereo")
	x = zoom(u,pole)
	y = zoom(v,pole)
	png(sprintf("%s/%s_ff.png",pngd,prefix))
	mapdom2(pole,x,x[,1],y[,1],main="u/v wind components")
	dev.off()
}

stat2array = function(lstat)
{
	# can be list() (ie no dim) because of bin files or missing dates
	dims = NULL

	for (l in lstat) {
		dims = dim(l[[1]])
		if (! is.null(dims)) break
	}

	stopifnot(! is.null(dims))

	data = array(NA_real_,dim=c(dims,length(lstat[[1]]),length(lstat)))

	for (i in seq(along=lstat)) {
		for (id in seq(along=lstat[[i]])) {
			if (is.list(lstat[[i]][[id]])) next

			data[,,,id,i] = simplify2array(lstat[[i]][[id]])
		}
	}

	data
}

rms = function(x,...)
{
	sqrt(mean(x^2,...))
}

rmx = function(x,...)
{
	c(rms=rms(x,...),ave=mean(x,...),q90=quantile(abs(x),prob=.9,...))
}

plotb = function(x,y,lty=1,col="bisque",pch="+",ylim=NULL,...)
{
	if (is.null(ylim)) ylim = range(boxplot(y,outline=FALSE,plot=FALSE)$stats)
	plot(x,apply(y,2,mean),type="l",ylim=ylim,...)
	w = diff(par("usr")[1:2])/25
	boxplot(y,col=col,lty=lty,pch=pch,pars=list(boxwex=w),add=TRUE,at=x,xaxt="n",yaxt="n")
}

plotbv = function(x,y,lty=1,col="bisque",pch="+",xlim=NULL,...)
{
	if (is.null(xlim)) xlim = range(boxplot(x,outline=FALSE,plot=FALSE)$stats)
	plot(apply(x,2,mean),y,type="l",xlim=xlim,...)
	w = diff(par("usr")[3:4])/30
	boxplot(x,col=col,lty=lty,pch=pch,pars=list(boxwex=w),horizontal=TRUE,add=TRUE,at=y,
		xaxt="n",yaxt="n")
}

matplott = function(x,y,col=1,pch="+",...)
{
	matplot(x,y,type="o",lty=1:3,col=col,pch=pch,...)
	abline(h=0,col="darkgrey")
	legend("topleft",c("RMS","bias","errx2"),col=col,lty=1:3,pch=pch,bg="transparent")
}

matplotv = function(x,y,col=1,pch="+",...)
{
	matplot(x,y,type="o",lty=1:3,col=col,pch=pch,...)
	abline(v=0,col="darkgrey")
	legend("topleft",c("RMS","bias","errx2"),col=col,lty=1:3,pch=pch,bg="transparent")
}

statarray = function(data,x,y)
{
	stopifnot(! is.null(dim(data)))
	array(dim=c(dim(data),length(x),length(y)),dimnames=c(dimnames(data),list(x,y)))
}

addrmxd = function(rmxd,data,nrmxd,id)
{
	if (! is.matrix(rmxd)) rmxd = matrix(rmxd,ncol=4)

	if (nrmxd == 0) {
		rmxd[,1] = data^2
		rmxd[,2] = data
		rmxd[,3] = abs(data)
		rmxd[,4] = 1
	} else {
		rmxd[,1] = rmxd[,1]+data^2
		rmxd[,2] = rmxd[,2]+data
		ind = which(abs(data) > rmxd[,3])
		rmxd[,3] = pmax(rmxd[,3],abs(data))
		rmxd[ind,4] = id
	}

	rmxd
}

saveZmean = function(filename,framep,paramsp,datesp,htimep,lzmeanp)
{
	if (! file.exists(filename)) {
		frlow = framep
		params = paramsp
		dates = datesp
		htime = htimep
		lzmeand = lzmeanp
		save(frlow,params,dates,htime,lzmeand,file=filename)
		return()
	}

	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lzmeand") %in% noms))

	if (class(frlow$g4) != class(framep$g4)) {
		cat("--> different classes:",class(frlow$g4),class(framep$g4),"\n")
		return()
	} else if (! identical(frlow$g4@nlong,framep$g4@nlong)) {
		cat("--> different grid\n")
		return()
	} else if (frlow$nlevel != framep$nlevel || ! identical(frlow$ilev,framep$ilev)) {
		cat("--> different levels, in file:",frlow$nlevel,framep$nlevel,"\n")
		return()
	} else if (! identical(dates,datesp)) {
		cat("--> different dates, in file:",dates,"\n")
		return()
	}

	indt = match(htimep,htime)

	if (identical(params,paramsp)) {
		inew = is.na(indt)

		for (j in seq(along=lzmeand)) {
			if (! all(inew)) {
				lzmeand[[j]][,,,na.omit(indt)] = lzmeanp[[j]][,,,which(! inew)]
			}

			if (any(inew)) {
				cat("--> add times to file:",htimep[inew],"\n")
				d = dim(lzmeand[[j]])
				dn = dimnames(lzmeand[[j]])
				zmeand = c(lzmeand[[j]],lzmeanp[[j]][,,,which(inew)])
				d[4] = d[4]+length(which(inew))
				if (! is.null(dn)) dn[[4]] = c(dn[[4]],dimnames(lzmeanp[[j]])[[4]][inew])
				lzmeand[[j]] = array(zmeand,dim=d,dimnames=dn)
			}
		}

		if (any(inew)) htime = c(htime,htimep[inew])
	} else {
		indp = match(paramsp,params)
		cat("--> different params, in file:",params,"\n")
		if (! identical(htime,htimep)) {
			cat("--> different steps, in file:",htime,"-> no saving\n")
			return()
		}

		inew = is.na(indp)
		if (! all(inew)) {
			cat("--> replace params already in file:",paramsp[! inew],"\n")
			lzmeand[na.omit(indp)] = lzmeanp[which(! inew)]
		}

		if (any(inew)) {
			cat("--> add params to file:",paramsp[inew],"\n")
			lzmeand = c(lzmeand,lzmeanp[inew])
			params = c(params,paramsp[inew])
		}
	}

	save(frlow,params,dates,htime,lzmeand,file=filename)
}

saveStat = function(filename,framep,paramsp,domsp,datesp,htimep,lstatdp)
{
	if (! file.exists(filename)) {
		frlow = framep
		params = paramsp
		doms = domsp
		dates = datesp
		htime = htimep
		lstatd = lstatdp
		save(frlow,params,doms,dates,htime,lstatd,file=filename)
		return()
	}

	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lstatd") %in% noms))
	if (! "doms" %in% noms) {
		cat("doms absent from file, use domsp\n")
		stopifnot(dim(lstatd[[1]])[3] == length(domsp))
		doms = domsp
	}

	if (class(frlow$g4) != class(framep$g4)) {
		cat("--> different classes:",class(frlow$g4),class(framep$g4),"\n")
		return()
	} else if (! identical(frlow$g4@nlong,framep$g4@nlong)) {
		cat("--> different grid\n")
		return()
	} else if (frlow$nlevel != framep$nlevel || ! identical(frlow$ilev,framep$ilev)) {
		cat("--> different levels, in file:",frlow$nlevel,framep$nlevel,"\n")
		return()
	} else if (! identical(doms,domsp)) {
		cat("--> different doms, in file:",doms,"\n")
		return()
	} else if (! identical(dates,datesp)) {
		cat("--> different dates, in file:",dates,"\n")
		return()
	}

	indt = match(htimep,htime)

	if (identical(params,paramsp)) {
		inew = is.na(indt)

		for (j in seq(along=lstatd)) {
			if (! all(inew)) {
				cat("--> replace times already in file:",htimep[! inew],"\n")
				lstatd[[j]][,,,,na.omit(indt)] = lstatdp[[j]][,,,,which(! inew)]
			}

			if (any(inew)) {
				cat("--> add times to file:",htimep[inew],"\n")
				d = dim(lstatd[[j]])
				dn = dimnames(lstatd[[j]])
				statd = c(lstatd[[j]],lstatdp[[j]][,,,,which(inew)])
				d[5] = d[5]+length(which(inew))
				if (! is.null(dn)) dn[[5]] = c(dn[[5]],dimnames(lstatdp[[j]])[[5]][inew])
				lstatd[[j]] = array(statd,dim=d,dimnames=dn)
			}
		}

		if (any(inew)) htime = c(htime,htimep[inew])
	} else {
		indp = match(paramsp,params)
		cat("--> different params, in file:",params,"\n")
		if (! identical(htime,htimep)) {
			cat("--> different steps, in file:",htime,"-> no saving\n")
			return()
		}

		inew = is.na(indp)
		if (! all(inew)) {
			cat("--> replace params already in file:",paramsp[! inew],"\n")
			lstatd[na.omit(indp)] = lstatdp[which(! inew)]
		}

		if (any(inew)) {
			cat("--> add params to file:",paramsp[inew],"\n")
			lstatd = c(lstatd,lstatdp[inew])
			params = c(params,paramsp[inew])
		}
	}

	save(frlow,params,doms,dates,htime,lstatd,file=filename)
}

loadZmean = function(filename,framep,paramsp,datesp,htimep)
{
	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lzmeand") %in% noms))

	if (class(frlow$g4) != class(framep$g4)) {
		cat("--> different classes:",class(frlow$g4),class(framep$g4),"\n")
		return(NULL)
	} else if (! identical(frlow$g4@nlong,framep$g4@nlong)) {
		cat("--> different grid\n")
		return(NULL)
	} else if (frlow$nlevel != framep$nlevel || ! identical(frlow$ilev,framep$ilev)) {
		cat("--> different levels:",frlow$nlevel,framep$nlevel,"\n")
		return(NULL)
	}

	indp = match(paramsp,params)
	if (any(is.na(indp))) {
		cat("--> some params missing:",paramsp[is.na(indp)],"\n")
		return(NULL)
	}

	indd = match(datesp,dates)
	if (any(is.na(indd))) {
		cat("--> some dates missing:",datesp[is.na(indd)],"\n")
		return(NULL)
	}

	indt = match(htimep,htime)
	if (any(is.na(indt))) {
		cat("--> some time-steps missing:",htimep[is.na(indt)],"\n")
		return(NULL)
	}

	lapply(lzmeand[indp],function(x) x[,,indd,indt,drop=FALSE])
}

loadStat = function(filename,framep,paramsp,domsp,datesp,htimep)
{
	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lstatd") %in% noms))
	if (! exists("doms")) {
		cat("doms absent, use domsp\n")
		stopifnot(dim(lstatd[[1]])[3] == length(domsp))
		doms = domsp
	} else if (is.list(doms)) {
		doms = names(doms)
	}

	if (class(frlow$g4) != class(framep$g4)) {
		cat("--> different classes:",class(frlow$g4),class(framep$g4),"\n")
		return(NULL)
	} else if (! identical(frlow$g4@nlong,framep$g4@nlong)) {
		cat("--> different grid\n")
		return(NULL)
	} else if (frlow$nlevel != framep$nlevel || ! identical(frlow$ilev,framep$ilev)) {
		cat("--> different levels:",frlow$nlevel,framep$nlevel,"\n")
		return(NULL)
	}

	indp = match(paramsp,params)
	if (any(is.na(indp))) {
		cat("--> some params missing:",paramsp[is.na(indp)],"\n")
		return(NULL)
	}

	indz = match(domsp,doms)
	if (any(is.na(indz))) {
		cat("--> some domains missing:",domsp[is.na(indz)],"\n")
		return(NULL)
	}

	indd = match(datesp,dates)
	if (any(is.na(indd))) {
		cat("--> some dates missing:",datesp[is.na(indd)],"\n")
		return(NULL)
	}

	indt = match(htimep,htime)
	if (any(is.na(indt))) {
		cat("--> some time-steps missing:",htimep[is.na(indt)],"\n")
		return(NULL)
	}

	lapply(lstatd[indp],function(x) x[,,indz,indd,indt,drop=FALSE])
}

Gparm = list(mar=c(1.8,1.8,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)
Gparh = list(mar=c(1.8,1.8,2.5,1),mgp=c(1.5,.5,0),tcl=-.3)
Gpart = list(mar=c(2.5,2.5,2.5,1),mgp=c(1.5,.5,0),tcl=-.3,cex=.8)
Gparmt = list(mar=c(2.5,2.5,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)
Gficbin = tempfile(fileext=".bin")

doms = readDom(sprintf("%s/config/domain.txt",Gdiag))
if (file.exists("config/domain.txt")) doms = readDom("config/domain.txt",doms)

descall = read.table(sprintf("%s/config/params.txt",Gdiag),header=TRUE)

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

nlatmax = 800
if (file.exists("config/setting.R")) {
	cat("--> reading constants:\n")
	source("config/setting.R",echo=TRUE,spaced=FALSE)
}

if ("fic" %in% names(cargs)) {
	dff = read.table(cargs$fic,header=TRUE)
} else {
	dff <- read.table("config/file.txt",header=TRUE)
}

dates = NA
dd = data.frame(date=dates,graph=TRUE)

if (file.exists("config/date.txt")) {
	dd = read.table("config/date.txt",skip=1,col.names=c("date","graph"))
	res = scan("config/date.txt",integer(),n=1,quiet=TRUE,comment.char="#")
	dates = dd$date
}

if (! "params" %in% names(cargs)) {
	params = scan("config/param.txt",what=character(),quiet=TRUE,comment.char="#")
} else if (length(cargs$params) == 1 && file.exists(cargs$params)) {
	params = scan(cargs$params,what=character(),quiet=TRUE,comment.char="#")
} else {
	params = cargs$params
}

ind = match(params,descall$faname)
if (any(is.na(ind))) stop("unknown parameters, see config/params.txt")

desc = descall[ind,]
npar = dim(desc)[1]

ind2 = match(c("ff","gradsp","gradl"),descall$symbol)
if (any(is.na(ind2))) stop("unknown parameters 2, see config/params.txt")
desc2 = descall[ind2,]

ilev = 0
if ("level" %in% names(cargs)) {
	if (length(cargs$level) == 1 && file.exists(cargs$level)) {
		ilev = scan(cargs$level,what=numeric(),quiet=TRUE)
	} else {
		ilev = as.numeric(cargs$level)
	}
} else if (file.exists("config/level.txt")) {
	ilev = scan("config/level.txt",what=numeric(),quiet=TRUE)
}

Gindl = which(as.integer((10*ilev)%%10) > 0)
ilev = as.integer(ilev)
cat("--> level indices/number (0: all levels):",ilev,"\n")

pngd = "."
if ("png" %in% names(cargs)) pngd = cargs$png

ref = "."
if ("ref" %in% names(cargs)) ref = cargs$ref

mapstat = TRUE
if ("mapstat" %in% names(cargs)) mapstat = cargs$mapstat == 1

if ("png" %in% names(cargs) || ! capabilities("X11")) {
	cat("--> no X11 device, sending plots to PNG files\n")
	if (! file.exists(pngd)) dir.create(pngd,recursive=TRUE)
} else {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
}

if (interactive()) browser()

lstat2 = lstat2d = lerr2 = vector("list",dim(desc2)[1])
names(lstat2) = names(lstat2d) = desc2$symbol
lstat = lstatd = lerr = lzmeand = vector("list",npar)
lstatr = NULL
if (! is.null(cargs$cmp)) {
	stopifnot(file.exists(cargs$cmp) && file.info(cargs$cmp)$isdir)
	lstatr = list()
}

htime = integer(length(dff$file))
tstep = rep(NA,length(dff$file))
frame.save = NULL
quiet = FALSE

nerr = matrix(0,nrow=npar,ncol=length(dff$file))
nerr2 = matrix(0,nrow=dim(desc2)[1],ncol=length(dff$file))

for (id in seq(along=dates)) {
	if (is.na(dates[id])) {
		fics = dff$file
	} else {
		fics = sprintf("%02d/%d/%s",res,dates[id],dff$file)
		cat("\nReading files for base:",dates[id],"\n")
	}

	if (any(! file.exists(fics))) {
		cat("-->",length(which(! file.exists(fics))),"files missing:\n",
			head(fics[! file.exists(fics)]),"...\n")
		stop("Error")
	}

	isbin = regexpr("\\.bin$",fics) > 0
	if (isbin[1]) stop("1st file is bin, not possible")

	for (i in seq(along=fics)) {
		cat("\nFile",fics[i],":\n")

		if (isbin[i]) {
			con = file(fics[i],"rb")
			grid = getGrid(con)
			eta = getVCoord(con)
			fc = getTime(con)
			ldata = getVars(con,grid,desc$faname,quiet)
			close(con)
			quiet = TRUE
			stopifnot(exists("frame",mode="list"))
			stopifnot(is.equal(grid,frame$g4) && length(eta) != frame$nlevel)
			frame$fc = fc
			frame$eta = eta
		} else {
			frame = getFrame(fics[i])
		}

		fc = frame$fc
		base = as.character(fc@base,"%Y%m%d R%H")
		if (! is.na(dates[id]) && as.character(fc@base,"%Y%m%d") != dates[id]) {
			message("dates mismatch: ",fc@base," != ",dates[id])
			stop("date mismatch")
		}

		g4 = frame$g4
		cat("--> grid type:",class(g4),"\n")
		hh = formatStep(fc@step,"h")
		s = sprintf("base/step: %s +%s - file: %s - graph: %s",base,hh,fics[i],dff$graph[i])

		if (is.na(tstep[i])) {
			cat("-->",s,"\n")
			tstep[i] = s
		} else {
			re = "base.+? \\+(.+?) -.+"
			stopifnot(sub(re,"\\1",s) == sub(re,"\\1",tstep[i]))
			# weird, but tstep is extended by n files after the 1st day
			tstep = c(tstep,s)
		}

		gg = dd$graph[id] && dff$graph[i]
		grmap = NULL

		if (gg && is(g4,"GaussGrid")) {
			mapf = dilat(g4)
			ilateq = equalize(mapf,offset=.5)
			#if (length(ilateq) < length(g4@nlong)) grmap = degrade(g4,ilateq)

			png(sprintf("%s/stretch%d.png",pngd,i))
			C = sqrt(mapf[length(mapf)]/mapf[1])
			tt = c("Stretching coefficient",sprintf("C = %g",round(C,2)))
			plot(mapf,type="l",main=tt,xlab="Latitude index",ylab="Stretching coef.",log="y")
			#points(ilateq,mapf[ilateq],pch=20,cex=.5)
			grid(col="darkgrey")
			dev.off()
		}

		if (is.null(frame.save)) {
			frame.save = frame
		} else {
			stopifnot(identical(frame$eta,frame.save$eta))
			stopifnot(frame$nlevel == frame.save$nlevel)
			if (! g4 == frame.save$g4) cat("--> test for different frame!\n")
		}

		if (length(g4@nlong) > nlatmax) {
			cat("--> halving initial resolution:",length(g4@nlong),max(g4@nlong),"\n")
			frlow = frame
			frlow$g4 = degrade(g4,nlat=2,nlon=2)
			frlow$npdg = sum(frlow$g4@nlong)
		} else {
			frlow = frame
		}

		stopifnot(all(abs(ilev) <= length(frlow$eta)))

		nlev = length(frlow$eta)
		if (length(ilev) == 1 && ilev < 0) {
			stopifnot(-ilev < nlev)

			cat("--> uniform eta selection of levels:",-ilev,"among",nlev,"\n")
			e = seq(min(frlow$eta),max(frlow$eta),length.out=-ilev)
			frlow$ilev = sapply(e,function(x) which.min(abs(x-frlow$eta)))
		} else if (length(ilev) < nlev && ilev > 0) {
			stopifnot(all(ilev %in% seq(nlev)))

			cat("--> selection of levels:",ilev,"\n")
			frlow$ilev = ilev
		} else {
			frlow$ilev = seq(nlev)
		}

		etai = frlow$eta[frlow$ilev]
		yeta = rev(range(frlow$eta))

		date = fc@base+fc@step
		cat("--> validity:",strftime(date,"%Y%m%d %H"),"\n")
		htime[i] = fc@step

		ficref = sprintf("%s/%s/%s",ref,strftime(date,"%H/%Y%m%d"),dff$ref[i])
		if (! file.exists(ficref) && fc@step > 0) {
			cat("--> no ref file",ficref,"- try base time\n")
			ficref = sprintf("%s/%s/%s",ref,strftime(fc@base,"%H/%Y%m%d"),dff$ref[i])
		}

		if (! file.exists(ficref)) {
			cat("--> no ref file",ficref,"\n")
		} else {
			cat("Ref file:",ficref,"\n")
			framo = getFrame(ficref)
		}

		datax = datay = dataz = datal = datam = datagpl = NULL
		dataox = dataoy = dataoz = NULL

		for (j in seq(npar)) {
			if (desc$ltype[j] == "-") {
				patt = desc$faname[j]
			} else if (length(frlow$ilev) > 1) {
				# selection requires same nlevel (interpolation impossible)
				if (frlow$nlevel == frame$nlevel && length(frlow$ilev) < frlow$nlevel) {
					patt = paste(sprintf("%03d",frlow$ilev),collapse="|")
					patt = sprintf("%s(%s)%s",desc$ltype[j],patt,desc$faname[j])
				} else {
					patt = sprintf("%s\\d+%s",desc$ltype[j],desc$faname[j])
				}
			} else if (desc$ltype[j] == "S") {
				patt = sprintf("S%03d%s",frlow$ilev,desc$faname[j])
			} else {
				cat("--> level type not 'S' for reading specified level\n")
				next
			}

			ss = desc$symbol[j]
			nom = desc$longname[j]

			if (isbin[i] && regexpr("^\\.",desc$faname[j]) > 0) {
				cat(". get field as diag",desc$faname[j],"\n")
				if (is.null(ldata[[j]])) {
					cat("--> field not found in .bin file\n")
					next
				}

				f = setDiagData(frame$g4,frlow$eta,frlow$ilev,ldata[[j]])
			} else if (! isbin[i] && regexpr("^\\.",desc$faname[j]) < 0) {
				cat(". get fields matching pattern",patt,"\n")
				f = getField(fics[i],patt,ss,frame,frlow)
				if (is.null(f)) {
					cat("--> field not found in FA/GRIB file\n")
					next
				}
			} else {
				cat("--> no field",desc$faname[j],"in",fics[i],"(OK)\n")
				next
			}

			if (! is.null(cargs$test)) {
				browser()
				ft = f
				#ft[,] = c(-1:1,0,-1)
				ft[,] = ft@grid@lat
				glow = degrade(frlow$g4)
				fl = interp(ft,glow,"linear")
				fp = interp(ft,glow,"ppp")
				nlat = length(glow@nlong)
				# suppress 1 lat line and 1 point on each remaining lat
				clats = c(0,cumsum(glow@nlong))
				glow@nlong = glow@nlong[-nlat]
				glow@nlong = as.integer(glow@nlong-1)
				glow@theta = 90-180*(seq(nlat-1)-.5)/(nlat-1)
				ind = lapply(seq(nlat-1),function(ilat) clats[ilat]+seq(glow@nlong[ilat]))
				ind = unlist(ind)
				glow@lat = glow@lat[ind]
				glow@long = glow@long[ind]
				fl = interp(ft,glow,"linear")
				fp = interp(ft,glow,"ppp")
				if (exists("framo")) {
					fi = interp(ft,framo$g4,"linear")
					gi = fi@grid
					stopifnot(isTRUE(identical(gi@pole,framo$g4@pole)))
					stopifnot(isTRUE(identical(gi@nlong,framo$g4@nlong)))
				}
			}

			nl = length(f@eta)
			if (desc$conv[j] != "") {
				cat("Convert data by:",desc$conv[j],"\n")
				fconv = eval(parse(text=desc$conv[j]))
				f[] = fconv(f)
			}

			cat("Primary checks\n")
			stopifnot(f@grid == frlow$g4)
			if (length(f@eta) > 1) stopifnot(all(f@eta == frlow$eta[frlow$ilev]))
			stopifnot(all(apply(f,2,function(x) all(! is.infinite(x)) && any(! is.na(x)))))

			#if (! all(apply(f,2,function(x) length(unique(na.omit(x))) > 1))) {
			if (length(unique(na.omit(as.vector(f)))) == 1) {
				cat("--> constant value of ",na.omit(as.vector(f))[1],"\n")
				next
			}

			cat("Compute stats and graphics for domains\n")
			prefix = character()
			if (gg) prefix = sprintf("%s+%s",as.character(fc@base,"%Y%m%d%H"),hh)
			tt = desc$longname[j]
			tt[2] = sprintf("on %s +%s",base,hh)
			qv = mapexi(f,grmap,desc[j,],doms,prefix,main=tt)
			if (is.null(lstat[[j]])) lstat[[j]] = statarray(qv,dates,dff$file)
			stopifnot(all(dim(qv) == dim(lstat[[j]])[1:3]))
			lstat[[j]][,,,id,i] = qv

			if (desc$save[j] == "X") {
				datax = f
			} else if (desc$save[j] == "Y") {
				datay = f
			} else if (desc$save[j] == "Z") {
				dataz = f
			} else if (desc$save[j] == "L") {
				datal = f
			} else if (desc$save[j] == "M") {
				datam = f
			} else if (desc$save[j] == "GPL" && is(g4,"GaussGrid")) {
				datagpl = gradl(f,f@grid,mc.cores=8)
			}

			if (! file.exists(ficref)) {
				rm(f,qv)
				gc()
				next
			}

			cat(". get ref field",desc$faname[j],"\n")
			if (regexpr("^\\.",desc$faname[j]) > 0) {
				if (! isbin[i]) next
				if (is.null(ldata[[j]])) next

				fo = setDiagData(framo$g4,frlow$eta,frlow$ilev,ldata[[j]])
			} else {
				if (isbin[i]) next

				if (desc$ltype[j] != "-") {
					# selection requires same nlevel (interpolation impossible)
					if (framo$nlevel == frlow$nlevel) {
						framo$ilev = frlow$ilev
					} else {
						ind = findInterval(frlow$eta[frlow$ilev],framo$eta)
						stopifnot(all(0 < ind & ind < framo$nlevel))
						# for linear interpolation, each level and the following one...
						framo$ilev = unique(sort(c(ind,ind+1)))
					}

					if (length(framo$ilev) < framo$nlevel) {
						patt = paste(sprintf("%03d",framo$ilev),collapse="|")
						patt = sprintf("%s(%s)%s",desc$ltype[j],patt,desc$faname[j])
					} else {
						patt = sprintf("%s\\d+%s",desc$ltype[j],desc$faname[j])
					}

					cat("--> new pattern for ref:",patt,"\n")
				}

				fo = getField(ficref,patt,ss,framo,frlow,cache.alt=TRUE)
			}

			if (desc$conv[j] != "") {
				cat("Convert data with:",desc$fconv[j],"\n")
				fo[] = fconv(fo)
			}

			cat("Primary checks\n")
			stopifnot(fo@grid == frlow$g4)
			if (length(fo@eta) > 1) stopifnot(all(fo@eta == frlow$eta[frlow$ilev]))
			stopifnot(all(apply(fo,2,function(x) all(! is.infinite(x)) && any(! is.na(x)))))

			if (desc$save[j] == "X") {
				dataox = fo
			} else if (desc$save[j] == "Y") {
				dataoy = fo
			} else if (desc$save[j] == "Z") {
				dataoz = fo
			}

			stopifnot(all(is.na(f) == is.na(fo)))
			fd = setDataPart(f,f-fo)

			cat("Compute stats and graphics for domains\n")
			tto = c(desc$longname[j],strftime(date,"valid: %Y%m%d %Hh"))
			mapexi(fo,grmap,desc[j,],doms,sprintf("diff%s",prefix),main=tto,prob=0:1)
			qv = mapexi(fd,grmap,desc[j,],doms)

			cat("Compute errors and zonal mean\n")
			zm = zonalmean(f)-zonalmean(fo)
			if (is.null(lstatd[[j]])) {
				lstatd[[j]] = statarray(qv,dates,dff$file)
				lzmeand[[j]] = statarray(zm,dates,dff$file)
				lerr[[j]] = statarray(fd,c("rmse","bias","errx","dayx"),dff$file)
			}

			stopifnot(all(dim(qv) == dim(lstatd[[j]])[1:3]))
			lstatd[[j]][,,,id,i] = qv
			lzmeand[[j]][,,id,i] = zm
			stopifnot(all(dim(fd) == dim(lerr[[j]])[1:2]))
			lerr[[j]][,,,i] = addrmxd(lerr[[j]][,,,i],fd,nerr[j,i],id)
			stopifnot(! is.null(dimnames(lstatd[[j]])))
			stopifnot(! is.null(dimnames(lerr[[j]])))

			nerr[j,i] = nerr[j,i]+1

			rm(fo,fd,qv,zm)
			gc()
		}

		if (! is.null(datax) && ! is.null(datay)) {
			cat(". compound field wind speed\n")
			ff = setDataPart(datax,sqrt(datax^2+datay^2))
			print(summary(as.vector(ff)))

			if (gg && is(g4,"GaussGrid") && pole(g4)[1] < 90) {
				cat("rotating wind to geo North\n")
				mapgpole(g4,datax,datay,sprintf("map%dpole",i))
			} else {
				u = datax
				v = datay
			}

			j2 = match("ff",desc2$symbol)
			tt = desc2$longname[j2]
			tt[2] = sprintf("on %s +%s",as.character(fc@base,"%Y%m%d R%H"),hh)
			qv = mapex2(datax,datay,doms,prefix,main=tt,palette="PuRd")
			if (is.null(lstat2[[j2]])) lstat2[[j2]] = statarray(qv,dates,dff$file)
			lstat2[[j2]][,,,id,i] = qv
			rm(datax,datay)
		}

		if (! is.null(datal) && ! is.null(datam)) {
			cat(". compound field modulus\n")
			gsp = setDataPart(datal,sqrt(datal^2+datam^2))

			j2 = match("gradsp",desc2$symbol)
			tt = desc2$longname[j2]
			tt[2] = sprintf("on %s +%s",as.character(fc@base,"%Y%m%d R%H"),hh)
			qv = mapexi(gsp,grmap,desc2[j2,],doms,prefix,main=tt)
			if (is.null(lstat2[[j2]])) lstat2[[j2]] = statarray(qv,dates,dff$file)
			lstat2[[j2]][,,,id,i] = qv
			rm(datal,datam)
		}

		if (! is.null(datagpl)) {
			cat(". computed field gradl\n")
			print(summary(as.vector(datagpl)))

			j2 = match("gradl",desc2$symbol)
			tt = desc2$longname[j2]
			tt[2] = sprintf("on %s +%s",as.character(fc@base,"%Y%m%d R%H"),hh)
			qv = mapexi(datagpl,grmap,desc2[j2,],doms,prefix,main=tt)
			if (is.null(lstat2[[j2]])) lstat2[[j2]] = statarray(qv,dates,dff$file)
			lstat2[[j2]][,,,id,i] = qv
			rm(datagpl)
		}

		if (! is.null(dataox) && ! is.null(dataoy)) {
			cat(". compound field wind speed\n")
			ffo = setDataPart(dataox,sqrt(dataox^2+dataoy^2))
			print(summary(as.vector(ffo)))

			if (is(g4,"GaussGrid") && pole(g4)[1] < 90) {
				cat("rotating wind to geo North\n")
				mapgpole(g4,dataox,dataoy,sprintf("mapdiff%dpole",i))
			}

			ffd = setDataPart(ff,ff-ffo)

			j2 = match("ff",desc2$symbol)
			tt = desc2$longname[j2]
			tt[2] = sprintf("on %s +%s",as.character(fc@base,"%Y%m%d R%H"),hh)
			qv = mapexi(ffd,grmap,desc2[j2,],doms,sprintf("diff%s",prefix),main=tt)
			if (is.null(lstat2d[[j2]])) {
				lstat2d[[j2]] = statarray(qv,dates,dff$file)
				lerr2[[j2]] = statarray(ffd,c("rmse","bias","errx","dayx"),dff$file)
			}

			lstat2d[[j2]][,,,id,i] = qv
			lerr2[[j2]][,,,i] = addrmxd(lerr2[[j2]][,,,i],ffd,nerr2[j2,i],id)
			nerr2[j2,i] = nerr2[j2,i]+1
			rm(dataox,dataoy)
		}

		gc()
	}
}

for (j in seq(along=lstat)) {
	stopifnot(! is.null(dimnames(lstat[[j]])))
	stopifnot(! is.null(dimnames(lstatd[[j]])))
	stopifnot(! is.null(dimnames(lerr[[j]])))
	stopifnot(! is.null(dimnames(lzmeand[[j]])))
}

if (! is.null(cargs$png)) {
	tstepgg = grep(" TRUE",tstep,value=TRUE)
	cat("Write the steps file:",length(tstep),"steps\n")
	write(tstepgg,file=sprintf("%s/steps.txt",cargs$png))
}

# those i2 (for desc2) may not be the same but the 1st one involves the others
# so desc should be OK
i2 = sapply(lstat2,is.array)
if (any(i2)) {
	cat("Add lstat2 to the statistics of forecasts\n")
	lstat = c(lstat,lstat2[i2])
	desc = rbind(desc,desc2[i2,])
}

i2 = sapply(lstat2d,is.array)
if (any(i2)) {
	cat("Add lstatd2 to the statistics of errors\n")
	lstatd = c(lstatd,lstat2d[i2])
	lerr = c(lerr,lerr2[i2])
	nerr = rbind(nerr,nerr2[i2,])
}

if (mapstat && any(sapply(lerr,is.array))) {
	cat("Maps and statistics of error (mean bias, RMSE, max and dayx error)\n")
	descb = desc
	descb$palette = "Blue-Red+"
	tt = desc$longname

	for (j in seq(along=lerr)) {
		cat(". param",desc$symbol[j],"\n")
		if (dim(lerr[[j]])[2] == 1) {
			fe = new("Field",grid=frlow$g4,eta=1)
		} else {
			fe = new("Field",grid=frlow$g4,eta=frlow$eta[frlow$ilev])
		}

		for (i in seq(along=dff$file)) {
			if (! dff$graph[i]) next

			frlow$fc = frame.save$fc
			frlow$fc@step = htime[i]
			fc = frlow$fc
			tt = descb$longname[j]
			tt[2] = sprintf("step +%s",formatStep(fc@step))

			prefix = character()
			if (dff$graph[i]) prefix = formatStep(fc@step,"h")

			cat(". bias, RMSE and errmax for",dff$file[i],"\n")
			bias = array(lerr[[j]][,,2,i]/nerr[j,i],dim(lerr[[j]])[1:2])
			fe = setDataPart(fe,bias)
			invisible(mapexi(fe,grmap,descb[j,],doms,sprintf("bias%s",prefix),main=tt,mnx=FALSE))

			rmse = array(sqrt(lerr[[j]][,,1,i]/nerr[j,i]),dim(lerr[[j]])[1:2])
			fe = setDataPart(fe,rmse)
			invisible(mapexi(fe,grmap,desc[j,],doms,sprintf("rmse%s",prefix),main=tt,mnx=FALSE))

			errx = array(lerr[[j]][,,3,i],dim(lerr[[j]])[1:2])
			fe = setDataPart(fe,errx)
			invisible(mapexi(fe,grmap,desc[j,],doms,sprintf("errx%s",prefix),main=tt,mnx=FALSE))

			cat(". index of maximum error\n")
			dayx = array(lerr[[j]][,,4,i],dim(lerr[[j]])[1:2])
			fe = setDataPart(fe,dayx)
			invisible(mapexi(fe,grmap,desc[j,],doms,sprintf("dayx%s",prefix),main=tt,mnx=FALSE))
		}
	}
}

cat("Save statistics\n")
saveStat("diag.RData",frlow,params,names(doms),dates,htime,lstatd)
saveZmean("diagz.RData",frlow,params,dates,htime,lzmeand)
for (j in seq(along=lstatd)) {
	ss = desc$symbol[j]
	fsave = sprintf("diag_%s.RData",ss)
	saveStat(fsave,frlow,params[j],names(doms),dates,htime,lstatd[j])
}

lzmeanr = NULL
if (! is.null(lstatr)) {
	ficsave = sprintf("%s/diag.RData",cargs$cmp)
	cat("Load cmp file",ficsave,"\n")
	lstatr = loadStat(ficsave,frlow,params,names(doms),dates,htime)
	if (is.null(lstatr)) cat("--> no stat in common from cmp\n")
	ficsave = sprintf("%s/diagz.RData",cargs$cmp)
	lzmeanr = loadZmean(ficsave,frlow,params,dates,htime)
}

if (length(htime) > 2 && diff(htime)[1] == 0) {
	htime[1] = -diff(htime[2:3])/3
	cat("--> changing time for graphics:",htime[1:2],"\n")
}

# units and labels
tlab = sprintf("Lead time (%s)",c("s","mn","h","days"))
unit = c(1,60,3600,86400)
mindiff = 5*c(0,60,3600,86400)
tunit = c("s","mn","h","d")
iu = findInterval(diff(range(htime)),mindiff)
tunit = tunit[iu]
tlab = tlab[iu]
ht = htime/unit[iu]
tfreq = 1
if (unit[iu] == 3600 && diff(range(ht)) > 12) tfreq = 3

ht6 = pretty(ht/tfreq,7)*tfreq
if (all((ht6-ht6[1])%%5 == 0)) {
	tfreq = 6
	ht6 = pretty(ht/tfreq,7)*tfreq
}

lat20 = pretty(sort(frlow$g4@theta),7)
etai = frlow$eta[frlow$ilev]
yeta = rev(range(frlow$eta))

ndom = length(doms)
nt = length(htime)
if (nt < length(ht)) length(ht) = nt

cat("Statistics of forecasts over",ndom,"domains\n")
for (j in seq(along=lstat)) {
	ss = desc$symbol[j]
	nom = desc$longname[j]
	cat(". param",nom,"\n")

	if (is.null(lstat[[j]])) {
		cat("--> no statistics of values to plot\n")
		next
	}

	nl = dim(lstat[[j]])[2]
	stopifnot(all(dim(lstat[[j]]) == c(101,nl,ndom,length(dates),nt)))

	# dims are [val5,lev,dom,time]
	qfc = apply(lstat[[j]],c(2,3,5),quantile,prob=c(0,2,5,8,10)/10,na.rm=TRUE)

	tt = sprintf("Statistics of %s",nom)
	nr = min(3,nl)
	nc = min(2,ndom)
	indl = round(seq(0,nl,length.out=nr+1)[-1])

	for (i in seq((ndom-1)%/%nc+1)-1) {
		png(sprintf("%s/stat%d_%s.png",pngd,i,ss))
		op = par(c(Gpart,list(mfcol=c(nr,nc))))

		for (id in 1:min(ndom-nc*i,nc)+nc*i) {
			tt[2] = sprintf("domain %s",names(doms)[id])
			il0 = 1

			for (il in indl) {
				if (nl > 1) {
					tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il0],
						frlow$ilev[il])
				}
				qfct = matrix(qfc[,il0:il,id,],ncol=length(ht))
				plotb(ht,qfct,main=tt,xlab=tlab,ylab=ss,xaxt="n")
				axis(1,ht6)
				il0 = il
			}
		}

		dev.off()
	}

	if (nl == 1) next

	indt = which(apply(qfc,4,function(x) any(! is.na(x))))
	if (length(indt) > 3) indt = indt[c(1,length(indt)%/%2,length(indt))]
	nr = min(2,ndom)
	nc = length(indt)

	for (i in seq((ndom-1)%/%nr+1)-1) {
		png(sprintf("%s/statv%d_%s.png",pngd,i,ss))
		op = par(c(Gpart,list(mfrow=c(nr,nc))))

		for (id in 1:min(ndom-nr*i,nr)+nr*i) {
			for (it in indt) {
				tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],tunit)
				plotbv(qfc[,,id,it],etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
			}
		}

		dev.off()
	}
}

if (any(sapply(lstatd,is.array))) {
	cat("Scores of forecasts over",ndom,"domains\n")
	for (j in seq(along=lstatd)) {
		ss = desc$symbol[j]
		nom = desc$longname[j]
		cat(". param",ss,"\n")

		if (is.null(lstatd[[j]])) {
			cat("--> no statistics of errors to plot\n")
			next
		}

		nl = dim(lstatd[[j]])[2]
		stopifnot(all(dim(lstatd[[j]]) == c(101,nl,ndom,length(dates),nt)))

		# dims are [qvalues,lev,dom,time]
		qe = apply(lstatd[[j]],c(2,3,5),quantile,prob=seq(0,50)/50,na.rm=TRUE)

		tt = sprintf("Error of %s",nom)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = round(seq(0,nl,length.out=nr+1)[-1])

		# qe err: ll box 3lx2
		for (i in seq((ndom-1)%/%nc+1)-1) {
			png(sprintf("%s/err%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				il0 = 1
				tt[2] = sprintf("domain %s",names(doms)[id])
				for (il in indl) {
					if (nl > 1) {
						tt[2] = sprintf("domain %s, levels [%d,%d]",names(doms)[id],
							frlow$ilev[il0],frlow$ilev[il])
					}
					qet = matrix(qe[,il0:il,id,],ncol=length(ht))
					plotb(ht,qet,main=tt,xlab=tlab,ylab=ss,xaxt="n")
					axis(1,ht6)
					il0 = il
				}
			}

			dev.off()
		}

		if (nl > 1) {
			# qe errv: box v 2x3t back to err, after scores
			tt = sprintf("Error of %s",nom)
			indt = which(apply(lstatd[[j]],5,function(x) any(! is.na(x))))
			if (length(indt) > 3) indt = indt[seq(1,length(indt),len=3)]
			nr = min(2,ndom)
			nc = length(indt)

			for (i in seq((ndom-1)%/%nr+1)-1) {
				png(sprintf("%s/errv%d_%s.png",pngd,i,ss))
				op = par(c(Gpart,list(mfrow=c(nr,nc))))

				for (id in 1:min(ndom-nr*i,nr)+nr*i) {
					it0 = 1
					for (it in indt) {
						tt[2] = sprintf("domain %s, lead-time [%g,%g]%s",names(doms)[id],
							ht[it0],ht[it],tunit)
						qet = qe[,,id,it0:it]
						if (length(dim(qet)) == 3) qet = matrix(aperm(qet,c(3,1:2)),ncol=nl)
						plotbv(qet,etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
						it0 = it
					}
				}

				dev.off()
			}
		}

		rmxv = apply(lstatd[[j]],c(2,3,5),rmx,na.rm=TRUE)
		rmxv = aperm(rmxv,c(2:4,1))
		if (! is.null(lstatr)) {
			rmxvr = apply(lstatr[[j]],c(2,3,5),rmx,na.rm=TRUE)
			rmxvr = aperm(rmxvr,c(2:4,1))
			rmxv = array(c(rmxv,rmxvr),c(dim(rmxv)[1:3],2*dim(rmxv)[4]))
		}

		tt = sprintf("Score of %s",nom)
		cols = rep(c("black","royalblue1"),each=3)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = as.integer(seq(1,nl,length.out=nr))

		# rmxv score: ll plot 3lx2
		for (i in seq((ndom-1)%/%nc+1)-1) {
			png(sprintf("%s/score%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				tt[2] = sprintf("domain %s",names(doms)[id])
				for (il in indl) {
					if (nl > 1) {
						tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il])
					}
					rmxt = matrix(rmxv[il,id,,],ncol=dim(rmxv)[4])
					matplott(ht,rmxt,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
					axis(1,ht6)
				}
			}

			dev.off()
		}

		if (nl > 1) {
			nr = min(2,ndom)
			nc = length(indt)

			# rmxv scorev: v plot 2x3t
			for (i in seq((ndom-1)%/%nr+1)-1) {
				png(sprintf("%s/scorev%d_%s.png",pngd,i,ss))
				op = par(c(Gpart,list(mfrow=c(nr,nc))))

				for (id in 1:min(ndom-nr*i,nr)+nr*i) {
					for (it in indt) {
						tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],
							tunit)
						matplotv(rmxv[,id,it,],etai,xlab=ss,ylab="eta",main=tt,
							xlim=range(rmxv[,id,,]),ylim=yeta,col=cols)
					}
				}

				dev.off()
			}
		}

		if (length(lzmeand) >= j) {
			# dims are [lat,lev,date,time]
			# rmxv score: ll plot 6l
			zrmxv = apply(lzmeand[[j]],c(2,4),rmx,na.rm=TRUE)
			zrmxv = aperm(zrmxv,c(2:3,1))
			if (! is.null(lzmeanr)) {
				zrmxvr = apply(lzmeanr[[j]],c(2,4),rmx,na.rm=TRUE)
				zrmxvr = aperm(zrmxvr,c(2:3,1))
				zrmxv = array(c(zrmxv,zrmxvr),c(dim(zrmxv)[1:2],2*dim(zrmxv)[3]))
			}

			tt = sprintf("Score zon. mean, %s",nom)
			nr = min(3,nl%/%2)
			nc = min(2,nl%/%2)
			if (nl == 1) nr = nc = 1
			indl = as.integer(seq(1,nl,length.out=nr*nc))

			png(sprintf("%s/scorez1_%s.png",pngd,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (il in indl) {
				if (nl > 1) tt[2] = sprintf("level %d",frlow$ilev[il])
				zt = matrix(zrmxv[il,,],nrow=length(ht))
				matplott(ht,zt,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
				axis(1,ht6)
			}

			dev.off()

			if (nl > 1) {
				indt = which(apply(lzmeand[[j]],4,function(x) any(! is.na(x))))
				nr = min(2,length(indt)%/%2)
				nc = min(3,length(indt)%/%2)
				if (length(indt) == 1) nr = nc = 1
				if (length(indt) > nr*nc) indt = indt[seq(1,length(indt),len=nr*nc)]

				png(sprintf("%s/scorevz1_%s.png",pngd,ss))
				op = par(c(Gpart,list(mfrow=c(nr,nc))))

				for (it in indt) {
					tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
					matplotv(zrmxv[,it,],etai,main=tt,ylim=yeta,xlab=ss,ylab="eta",col=cols)
				}

				dev.off()

				nr = min(3,length(indt))
				nc = 1
				indt = as.integer(seq(1,length(htime),len=nr))
				lat = sort(frlow$g4@theta)
				zrmsev = apply(lzmeand[[j]],c(1:2,4),rms,na.rm=TRUE)

				png(sprintf("%s/scorevz2_%s.png",pngd,ss))
				op = par(c(Gparmt,list(mfrow=c(nr,nc))))

				for (it in indt) {
					tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
					plotv(lat,etai,zrmsev[,,it],main=tt,ylim=yeta,xlab="Latitude",ylab="eta",
						xaxt="n",las=0)
					axis(1,lat20)
				}

				dev.off()
			}
		}

		if (length(dates) < 2) next

		cat("Scores for days\n")
		indt = which(apply(qe,4,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[seq(1,length(indt),len=3)]

		dd = as.Date(as.character(dates),format="%Y%m%d")
		indd = order(dd)
		nr = length(indt)
		nc = min(2,ndom)
		rmxt = apply(lstatd[[j]],3:5,rmx,na.rm=TRUE)
		rmxt = aperm(rmxt,c(2:4,1))
		if (! is.null(lstatr)) {
			rmxr = apply(lstatr[[j]],3:5,rmx,na.rm=TRUE)
			rmxr = aperm(rmxr,c(2:4,1))
			rmxt = array(c(rmxt,rmxr),c(dim(rmxt)[1:3],2*dim(rmxt)[4]))
		}

		# rmxt scoret: day t plot 3tx2
		for (i in seq((ndom-1)%/%nc+1)-1) {
			png(sprintf("%s/scoret%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				ylim = range(rmxt[id,,,],na.rm=TRUE)

				for (it in indt) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],
						tunit)
					matplott(dd[indd],rmxt[id,indd,it,],ylim=ylim,main=tt,xlab="Date",
						ylab=ss,col=cols,xaxt="n")
					axis.Date(1,dd,format="%Y%m%d")
				}
			}

			dev.off()
		}

		if (nl > 1) {
			# rmxvd rmsev: day v image 3tx2
			ttd = sprintf("RMSE of %s",nom)
			rmxvd = apply(lstatd[[j]],2:5,rms,na.rm=TRUE)

			for (i in seq((ndom-1)%/%nc+1)-1) {
				png(sprintf("%s/rmsevt%d_%s.png",pngd,i,ss))
				op = par(c(Gparmt,list(mfcol=c(nr,nc))))

				for (id in 1:min(ndom-nc*i,nc)+nc*i) {
					br = prettyBreaks(rmxvd[,id,,indt[1:nr]],crop=TRUE)$breaks

					for (it in indt) {
						ttd[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],
							tunit)
						plotv(dd[indd],etai,t(rmxvd[,id,indd,it]),br,main=ttd,ylim=yeta,
							xlab="Date",ylab="eta",xaxt="n",las=0)
						axis.Date(1,dd,format="%Y%m%d")
					}
				}

				dev.off()
			}
		}
	}
}
