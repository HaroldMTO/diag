Gdiag = "~/util/diag"

library(maps)
library(mapproj)
library(parallel)
library(mffield)

Gvp0 = 101325

source(sprintf("%s/plot.R",Gdiag))

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
	ficbin = tempfile(fileext=".bin")
	out = system(sprintf("epy_dump.py %s -f frame -o %s",file,ficbin),intern=TRUE)
	out = grep("dimensions",out,ignore.case=TRUE,value=TRUE)
	cat(out[1],"\n")

	con = file(ficbin,"rb")
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
	file.remove(ficbin)

	if (gauss) {
		theta = equiLat(length(nlong))
		g4 = new("GaussGrid",lat=lats,long=longs,nlong=nlong,pole=gem[2:3],stretch=gem[1],
			theta=theta)
		g4@theta = csLat(g4)

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

FApattern = function(ltype,faname,ilev,selev=FALSE)
{
	if (ltype == "-") {
		patt = faname
	} else if (length(ilev) > 1) {
		# selection requires same nlevel (interpolation impossible)
		if (selev) {
			patt = paste(sprintf("%03d",ilev),collapse="|")
			patt = sprintf("%s(%s)%s",ltype,patt,faname)
		} else {
			patt = sprintf("%s\\d+%s",ltype,faname)
		}
	} else if (ltype == "S") {
		patt = sprintf("S%03d%s",ilev,faname)
	} else {
		cat("--> level type not 'S' for reading specified level\n")
		return(NULL)
	}

	patt
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

getField = function(fic,param,symbol,frame,frlow=frame,cache.alt=FALSE,mc.cores=1)
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
	# no, don't test alternative file
	if (FALSE && file.exists(fsave)) {
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
		ilev = 0
		fic.save = ""
		noms = load(fsave)
		stopifnot(all(c("data","ilev","step","fic.save") %in% noms))

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
		data = getGPFields(fic,param,sum(grid@nlong))
	}
	stopifnot(dim(data)[1] == length(grid))

	ilev = frlow$ilev
	selev = FALSE
	if (dim(data)[2] == 1) {
		eta = 1
	} else if (dim(data)[2] == length(ilev)) {
		eta = frlow$eta[ilev]
	} else {
		# data contain a selection of levels (ilev, built by way of param, ie pattern)
		selev = TRUE
		# data should have full (or reduced) nlev (if not, we don't know its levels)
		eta = frame$eta
		if (! is.null(frame$ilev)) eta = eta[frame$ilev]
		stopifnot(dim(data)[2] == length(eta))
	}

	f = new("Field",grid=grid,eta=eta)
	f = setDataPart(f,data)
	if (! lread) return(f)

	# ppp is just for fun...
	meth = "linear"
	if (length(f) < 2e6) {
		h = hist(f,plot=FALSE)
		if (length(which(h$counts > 0)) == 2) {
			stopifnot(length(table(f)) == 2)
			meth = "ppp"
		}
	}

	if (sum(grid@nlong) != npdg) {
		nlat = length(grid@nlong)
		cat("--> interpolation from",sum(grid@nlong),nlat,"grid-points/lats - method:",
			meth,"\n")
		it = system.time(f <- interp(f,frlow$g4,meth,mc.cores))
		if (prof) cat("interpolation time:",it,"\n")
	}

	if ("ind" %in% names(frlow)) {
		cat("--> selection of grid-points:",length(frlow$g4),"/",npdg,"\n")
		f = setDataPart(f,f[frlow$ind,,drop=FALSE])
	}

	# if selev, dim(f)[2] > 1 && length(f@eta) != length(ilev)
	if (selev) {
		cat("--> interpolation from",length(f@eta),"levels - method:",meth,"\n")
		f = interpAB(f,frlow$eta[ilev],method="quad")
	}

	data = getDataPart(f)
	step = fc@step
	fic.save = fic
	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	st <- system.time(save("data","ilev","step","fic.save",file=fsave))
	if (prof) cat("caching time:",st,"\n")

	f
}

getGPFields = function(file,field,npdg)
{
	ficbin = tempfile(fileext=".bin")
	cmd = sprintf("epy_dump.py %s -f '%s' -o %s",file,field,ficbin)
	rt = system.time(out <- system(cmd,intern=TRUE))
	if (prof) cat("FAreading time:",rt,"\n")
	out = grep("date|time|step|gridpoint size",out,ignore.case=TRUE,value=TRUE)
	cat(paste(out,collapse=", "),"\n")

	con = file(ficbin,"rb")
	nl = readBin(con,"integer",1,size=8)
	stopifnot(length(nl) == 1)

	data = matrix(nrow=npdg,ncol=nl,dimnames=list(NULL,seq(nl)))
	rt = system.time(for (j in seq(nl)) data[,j] = readBin(con,"numeric",npdg))
	if (prof) cat("Binreading time:",rt,"\n")
	stopifnot(nl == readBin(con,"integer",1,size=8))
	close(con)
	file.remove(ficbin)

	data
}

getSPFields = function(file,field,nwave)
{
	ficbin = tempfile(fileext=".bin")
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,ficbin),ignore.stdout=TRUE)
	con = file(ficbin,"rb")
	nsp2 = readBin(con,"integer",1,size=8)
	data = readBin(con,"numeric",nsp2)
	close(con)
	file.remove(ficbin)

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

mapsegments = function(dom,lat,long,...)
{
	# S->N cross-section
	map.grid(c(long,long,dom@ylim[1],dom@ylim[2]),nx=1,labels=FALSE,pretty=FALSE,...)

	# W->E cross section
	map.grid(c(dom@xlim[1],dom@xlim[2],lat,lat),ny=1,labels=FALSE,pretty=FALSE,...)
}

mapdom = function(dom,grid,ind,data,main=NULL,scale=TRUE,mar=c(2,2,3,5),mgp=c(2,1,0),...)
{
	# mar must be set before calling map AND passed to map
	# (because map resets mar internally and resets it on exit: this is then fake mar!)
	par(mar=mar,mgp=mgp)
	l = mapxy(dom,mar=mar,new=TRUE)
	box()

	if (scale) {
		ymax = max(abs(data),na.rm=TRUE)
		scal = 10^-round(log10(ymax/1.5))
		if (.001 <= scal && scal < 1 || is.infinite(scal)) scal = 1
		if (scal != 1) {
			data = scal*data
			if (length(main) > 0) {
				main[1] = sprintf("%s - scaling: %g",main[1],scal)
			} else {
				message(sprintf("data scaling: %g",scal))
			}
		}
	}

	mappoints(grid,ind,data,...)

	lines(l)
	title(main)
}

mapdom2 = function(dom,points,zx,zy,main=NULL,colvec=TRUE,length=.05,angle=15,ppi=6,
	mar=c(2,2,3,5),mgp=c(2,1,0),quiet=FALSE,...)
{
	# mar must be set before calling map AND passed to map
	# (because map resets mar internally and resets it on exit: this is then fake mar!)
	par(mar=mar,mgp=mgp)
	l = mapxy(dom,mar=mar,new=TRUE)
	box()

	if (ppi > 96) stop("ppi > 96")

	nppi = prod(par("fin"))*ppi
	npmax = min(as.integer(nppi),.Machine$integer.max)

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
		maplegend(levels,col=palette)
	}

	lines(l)
	title(main)
}

hist.annot = function(data,nmin=10,n=20,split=TRUE,info=TRUE,...)
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

	# print min/max under the plot
	znx = range(data)
	text(znx[1],-.02*y,sprintf("| %.3g",min(data,na.rm=TRUE)),adj=0,col=4)
	text(znx[2],-.02*y,sprintf("%.3g |",max(data,na.rm=TRUE)),adj=1,col=4)

	if (! info) return(NULL)

	# print mean and sd near the top of the plot
	m = mean(data,na.rm=TRUE)

	s = sd(data,na.rm=TRUE)
	segments(m,.9*y,m,y,col=4)
	arrows(m-s,.95*y,m+s,.95*y,length=.03,angle=90,code=3,col=4)

	dx = diff(range(h$breaks))/15

	if (m+dx > h$breaks[1]+.9*diff(range(h$breaks))) {
		text(m-dx,.97*y,sprintf("m: %.3g",m),adj=1,col=4)
		text(m-dx,.92*y,sprintf("m/sd: %.3g",m/s),adj=1,col=4)
	} else {
		text(m+dx,.97*y,sprintf("m: %.3g",m),adj=0,col=4)
		text(m+dx,.92*y,sprintf("m/sd: %.3g",m/s),adj=0,col=4)
	}
}

plotsectiongeo = function(field,dom,main,...)
{
	longm = round(mean(dom@xlim),1)
	latm = round(mean(dom@ylim),1)

	yz = sectiongeo(field,longm,dom@ylim)
	stopifnot(all(is.finite(yz$data)))

	main[2] = sprintf("S->N section, longitude %g",longm)
	plotv(yz$lats,field@eta,yz$data,main=main,xlab="",ylab="",...)
	stopifnot(all(abs(range(yz$lats)) <= 90))

	xz = sectiongeo(field,dom@xlim,latm)
	stopifnot(all(is.finite(xz$data)))

	main[2] = sprintf("W->E section, latitude %g",latm)
	plotv(xz$longs,field@eta,xz$data,main=main,xlab="",ylab="",long=TRUE,...)
	if (any(abs(range(xz$longs)) > 180)) cat("--> longs:",rangeLong(xz$longs),"\n")
}

mapext = function(d4,xy,prefix,dom,desc,ext,main,mc.cores=1,...)
{
	if (ext == "min") {
		pal = "Heat+"
	} else if (ext == "max") {
		pal = "Heat"
	} else {
		stop(sprintf("function '%s' must be min or max",ext))
	}

	fun = get(ext)
	ficpng = sprintf("%s/map%s%s%s_%s.png",pngd,substr(ext,3,3),prefix,dom,desc$symbol)
	if (file.exists(ficpng)) return()

	nl = dim(xy)[2]
	il = (nl+1)%/%2

	png(ficpng)
	par(mfrow=c(2,2),cex=.7)

	par(Gparm)
	if (mc.cores == 1) {
		data = apply(xy[,1:il],1,fun)
	} else {
		lxy = mclapply(seq(dim(xy)[1]),function(i) fun(xy[i,1:il]),mc.cores=mc.cores)
		data = simplify2array(lxy)
	}

	tt[2] = sprintf("%s %s. of lev. [1,%d]",main,ext,il)
	mapdom(d4,xy@grid,xy@ind,data,main=tt,palette=pal,...)

	par(mar=Gparh$mar)
	hist.annot(data,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

	par(Gparm)
	if (mc.cores == 1) {
		data = apply(xy[,-(1:il)],1,fun)
	} else {
		lxy = mclapply(seq(dim(xy)[1]),function(i) fun(xy[i,-(1:il)]),mc.cores=mc.cores)
		data = simplify2array(lxy)
	}

	tt[2] = sprintf("%s %s. of lev. [%d,%d]",main,ext,il+1,nl)
	mapdom(d4,xy@grid,xy@ind,data,main=tt,palette=pal,...)

	par(mar=Gparh$mar)
	stopifnot(all(par("mar") == Gparh$mar))
	hist.annot(data,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

	dev.off()
}

mapexi = function(f,desc,doms,prob=seq(0,100)/100,prefix=character(),main,np=4,
	mnx=TRUE,force=FALSE,mc.cores=1,...)
{
	nl = length(f@eta)
	qv = array(dim=c(length(prob),nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))

	for (i in seq(along=doms)) {
		d4 = doms[[i]]
		dom = names(doms)[i]
		ind = inDomain(f@grid,d4)
		if (length(which(ind)) < np) {
			cat("--> less than",np,"points in domain\n")
			next
		}

		if (is(f@grid,"GaussGrid") && diff(d4@xlim) > 350 && diff(d4@ylim) > 170) {
			# no domain, no zoom
			domin = d4
			xy = f
		} else {
			pts = f@grid[ind]
			domin = new("Domain",long=pts@long,lat=pts@lat,proj=d4@proj,param=d4@param)
			if (area(domin,d4) < .1) {
				cat("--> area too small for domain\n")
				next
			}

			xy = zoom(f,d4)
		}

		stopifnot(all(is.finite(xy)))
		#qv[,,i] = apply(xy,2,quantile,prob=prob)
		lqv = mclapply(seq(dim(xy)[2]),function(i) quantile(xy[,i],prob),mc.cores=mc.cores)
		qv[,,i] = simplify2array(lqv)

		if (length(prefix) == 0) next

		cat(".. domain:",dom,"- points/total:",dim(xy)[1],dim(f)[1],"\n")

		ficpng = sprintf("%s/map%s%s_%s.png",pngd,prefix,dom,desc$symbol)
		if (! file.exists(ficpng) || force) {
			tt = main
			png(ficpng)

			if (nl == 1) {
				op = par(c(Gparm,list(cex=.8)))
				mapdom(d4,xy@grid,xy@ind,xy[,1],main=tt,palette=desc$palette,ppi=96,...)
			} else {
				op = par(c(Gparm,list(mfrow=c(2,2),cex=.7)))

				# be careful of nl=2...
				indl = as.integer(seq(1,nl,length.out=min(4,nl)))
				if (length(Gindl) == 2) indl = c(1,Gindl,nl)
				pt = system.time(for (il in indl) {
					tt[2] = sprintf("%s lev %d (eta %g)",main[2],il,round(f@eta[il],2))
					mapdom(d4,xy@grid,xy@ind,xy[,il],main=tt,palette=desc$palette,quiet=TRUE,...)
					mapsegments(domin,lat=mean(domin@ylim),long=mean(domin@xlim),col="darkgrey")
				})
				if (prof) cat("map time:",pt,"\n")
			}

			dev.off()
		}

		ficpng = sprintf("%s/hist%s%s_%s.png",pngd,prefix,dom,desc$symbol)
		if (! file.exists(ficpng) || force) {
			tt[2] = main[2]
			png(ficpng)

			if (nl == 1) {
				op = par(c(Gparh,list(cex=.8)))
				hist.annot(xy,col="whitesmoke",main=tt,xlab=NULL,ylab=NULL)
			} else {
				op = par(c(Gparm,list(mfcol=c(2,2)),cex=.7))
				st = system.time(plotsectiongeo(f,d4,main=desc$longname,palette=desc$palette))
				if (prof) cat("section time:",st,"\n")

				par(mar=Gparh$mar)
				hist.annot(xy,col="whitesmoke",main=tt,xlab=NULL,ylab=NULL)
				plotbv(qv[,,i],xy@eta,main=tt[1],xlab="",ylab="",ylim=rev(range(xy@eta)))
				abline(v=0,col="darkgrey",lty=2)
			}

			dev.off()
		}

		if (! mnx || nl < 4) next

		et = system.time(mapext(d4,xy,prefix,dom,desc,ext="min",main[2],mc.cores=mc.cores,...))
		if (prof) cat("extreme time:",et,"\n")
		mapext(d4,xy,prefix,dom,desc,ext="max",main[2],mc.cores=mc.cores,...)
	}

	qv
}

mapex2 = function(u,v,doms,prefix=NULL,main,np=4,palette,force=FALSE,...)
{
	nl = length(u@eta)
	qv = array(dim=c(101,nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))

	for (i in seq(along=doms)) {
		d4 = doms[[i]]$d4
		dom = names(doms)[i]
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

		tt = "u/v wind"
		cat(".. domain:",dom,"\n")

		ficpng = sprintf("%s/map%s%s_ff.png",pngd,prefix,dom)
		if (! file.exists(ficpng) || force) {
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
					tt[2] = sprintf("%s lev %d (eta %g)",main[2],il,round(u@eta[il],2))
					mapdom2(doms[[i]],x,x[,il],y[,il],main=tt,...)
					mapsegments(domin,lat=mean(domin@ylim),long=mean(domin@xlim),col="darkgrey")
				}
			}

			dev.off()
		}

		ficpng = sprintf("%s/hist%s%s_ff.png",pngd,prefix,dom)
		if (! file.exists(ficpng) || force) {
			png(ficpng)

			tt = "wind speed"
			if (nl == 1) {
				par(c(Gparh,list(cex=.8)))
				hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab="wind speed")
			} else {
				ff = sqrt(u^2+v^2)
				ff = setDataPart(u,ff)
				par(c(Gparm,list(mfcol=c(2,2),cex=.7)))
				plotsectiongeo(ff,doms[[i]],main=tt[1],palette=desc$palette)

				par(mar=Gparh$mar)
				hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
				plotbv(qv[,,i],u@eta,main=tt[1],xlab="",ylab="",ylim=rev(range(u@eta)))
				abline(v=0,col="darkgrey",lty=2)
			}

			dev.off()
		}
	}

	qv
}

statarray = function(data,x,y)
{
	stopifnot(! is.null(dim(data)))
	array(dim=c(dim(data),length(x),length(y)),dimnames=c(dimnames(data),list(x,y)))
}

addrmxd = function(rmxd,data,nrmxd,id,parallel=FALSE)
{
	if (! is.matrix(rmxd)) rmxd = matrix(rmxd,ncol=4)

	if (nrmxd == 0) {
		rmxd[,1] = data^2
		rmxd[,2] = data
		rmxd[,3] = abs(data)
		rmxd[,4] <- 1
	} else if (parallel) {
		# dependence on rmxd[,3]: rmxd[,4] must be done before
		ind = which(abs(data) > rmxd[,3])
		rmxd[ind,4] = id

		t1 = mcparallel(rmxd[,1]+data^2)
		t2 = mcparallel(rmxd[,2]+data)
		t3 = mcparallel(pmax(rmxd[,3],abs(data)))

		lc = mccollect(list(t1,t2,t3))
		try(parallel:::mckill(list(t1,t2,t3),15),silent=TRUE)
		for (i in 1:3) rmxd[,i] = lc[[i]]
		rm(lc)
	} else {
		# dependence on rmxd[,3]: rmxd[,4] must be done before
		ind = which(abs(data) > rmxd[,3])
		rmxd[ind,4] = id
		rmxd[,1] = rmxd[,1]+data^2
		rmxd[,2] = rmxd[,2]+data
		rmxd[,3] = pmax(rmxd[,3],abs(data))
	}

	rmxd
}

saveStat = function(filename,framep,paramsp,domsp,datesp,htimep,lstatdp,lzmeanp)
{
	if (! file.exists(filename)) {
		frlow = framep
		params = paramsp
		doms = domsp
		dates = datesp
		htime = htimep
		lstatd = lstatdp
		lzmeand = lzmeanp
		save(frlow,params,doms,dates,htime,lstatd,lzmeand,file=filename)
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
				lzmeand[[j]][,,,na.omit(indt)] = lzmeanp[[j]][,,,which(! inew)]
			}

			if (any(inew)) {
				cat("--> add times to file:",htimep[inew],"\n")
				d = dim(lstatd[[j]])
				dn = dimnames(lstatd[[j]])
				statd = c(lstatd[[j]],lstatdp[[j]][,,,,which(inew)])
				d[5] = d[5]+length(which(inew))
				if (! is.null(dn)) dn[[5]] = c(dn[[5]],dimnames(lstatdp[[j]])[[5]][inew])
				lstatd[[j]] = array(statd,dim=d,dimnames=dn)
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
			lstatd[na.omit(indp)] = lstatdp[which(! inew)]
			lzmeand[na.omit(indp)] = lzmeanp[which(! inew)]
		}

		if (any(inew)) {
			cat("--> add params to file:",paramsp[inew],"\n")
			lstatd = c(lstatd,lstatdp[inew])
			lzmeand = c(lzmeand,lzmeanp[inew])
			params = c(params,paramsp[inew])
		}
	}

	save(frlow,params,doms,dates,htime,lstatd,lzmeand,file=filename)
}

interptest = function(f)
{
	f[,] = f@grid@lat
	glow = degrade(f@grid)
	fl = interp(f,glow,"linear")
	fp = interp(f,glow,"ppp")
	nlat = length(glow@nlong)

	# suppress 1 lat line and 1 point on each remaining lat
	clats = c(0,cumsum(glow@nlong))
	glow@nlong = glow@nlong[-nlat]
	glow@nlong = as.integer(glow@nlong-1)
	glow@theta = equiLat(nlat-1)
	ind = lapply(seq(nlat-1),function(ilat) clats[ilat]+seq(glow@nlong[ilat]))
	ind = unlist(ind)
	glow@lat = glow@lat[ind]
	glow@long = glow@long[ind]

	# compare f to fl and fp (plot, hist or anything else)
	fl = interp(f,glow,"linear")
	fp = interp(f,glow,"ppp")

	return(NULL)
}

interptesto = function(f,grid)
{
	fi = interp(f,grid,"linear")
	gi = fi@grid
	stopifnot(isTRUE(identical(gi@pole,grid@pole)))
	stopifnot(isTRUE(identical(gi@nlong,grid@nlong)))
}

Gparm = list(mar=c(1.8,1.8,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)
Gparh = list(mar=c(1.8,1.8,2.5,1),mgp=c(1.5,.5,0),tcl=-.3)
Gpart = list(mar=c(2.5,2.5,2.5,1),mgp=c(1.5,.5,0),tcl=-.3,cex=.8)
Gparmt = list(mar=c(2.5,2.5,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)

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
cat("--> parameters:",params,"\n")

ind = match(params,descall$faname)
if (any(is.na(ind))) stop("unknown parameters, see config/params.txt")

desc = descall[ind,]
npar = dim(desc)[1]
descb = desc
descb$palette = "Blue-Red+"

ind2 = match(c("ff","gradsp","gradl"),descall$symbol)
if (any(is.na(ind2))) stop("unknown parameters 2, see config/params.txt")
desc2 = descall[ind2,]
desc2b = desc2
desc2b$palette = "Blue-Red+"

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

ref = "."
if ("ref" %in% names(cargs)) ref = cargs$ref

mapstat = TRUE
if ("mapstat" %in% names(cargs)) mapstat = cargs$mapstat == 1

prof = FALSE
if ("prof" %in% names(cargs)) prof = as.logical(cargs$prof)

wide = FALSE
if ("wide" %in% names(cargs)) wide = as.logical(cargs$wide)

etahigh = 0
if ("etahigh" %in% names(cargs)) etahigh = as.numeric(cargs$etahigh)

pngd = ""
if ("png" %in% names(cargs)) pngd = as.character(cargs$png)

if (! is.na(pngd) && nzchar(pngd) || ! capabilities("X11")) {
	cat("--> sending plots to PNG files\n")
	if (! file.exists(pngd)) dir.create(pngd,recursive=TRUE)
} else {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
}

ff = list(vars=c("X","Y"),formula=function(x,y) sqrt(x^2+y^2),symbol="ff")
gsp = list(vars=c("L","M"),formula=function(x,y) sqrt(x^2+y^2),symbol="gradsp")
#gl = list(vars="GPL",formula=gradl,symbol="gradl")
compound = list("wind speed"=ff,"gradient module"=gsp)

if (interactive()) browser()

lstat2 = lstat2d = lerr2 = vector("list",dim(desc2)[1])
names(lstat2) = names(lstat2d) = desc2$symbol
lstat = lstatd = lerr = lzmeand = vector("list",npar)

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
			stopifnot(exists("frame.save",mode="list"))
			stopifnot(is.equal(grid,frame.save$g4) && isTRUE(all.equal(eta),frame.save$eta))
			frame$fc = fc
			frame$eta = eta
		} else {
			frame = getFrame(fics[i])
			ldata = vector("list",length(desc$faname))
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
		s = sprintf("base/step: %s +%s - file: %s",base,hh,fics[i])

		date = fc@base+fc@step
		cat("--> validity:",strftime(date,"%Y%m%d %H"),"\n")
		htime[i] = fc@step

		ficref = sprintf("%s/%s/%s",ref,strftime(date,"%H/%Y%m%d"),dff$ref[i])
		if (! file.exists(ficref) && fc@step > 0) {
			cat("--> no ref file",ficref,"- try base time\n")
			ficref = sprintf("%s/%s/%s",ref,strftime(fc@base,"%H/%Y%m%d"),dff$ref[i])
		}

		if (file.exists(ficref)) s = sprintf("%s - ref: %s",s,dff$ref[i])
		s = sprintf("%s - graph: %s",s,dff$graph[i])

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
		prefixp = character()
		if (gg) prefixp = sprintf("%s+%s",as.character(fc@base,"%Y%m%d%H"),hh)

		if (gg && is(g4,"GaussGrid") && ! is.null(cargs$test)) {
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

		frlow = frame

		if (length(g4@nlong) > nlatmax) {
			cat("--> halving initial resolution:",length(g4@nlong),max(g4@nlong),"\n")
			frlow$g4 = degrade(g4,nlat=2,nlon=2)
			frlow$npdg = sum(frlow$g4@nlong)
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

		if (! file.exists(ficref)) {
			cat("--> no ref file",ficref,"\n")
			framo = NULL
		} else {
			cat("Ref file:",ficref,"\n")
			framo = getFrame(ficref)

			# selection requires same nlevel (interpolation impossible)
			if (isTRUE(all.equal(framo$eta,frlow$eta))) {
				framo$ilev = frlow$ilev
			} else {
				ind = findInterval(frlow$eta[frlow$ilev],framo$eta)
				stopifnot(all(0 < ind & ind <= framo$nlevel))
				# for linear interpolation, each level and the following one
				# for quad interpolation, each level and a few surrounding ones
				ind = unique(sort(c(ind-1,ind,ind+1,ind+2)))
				framo$ilev = ind[0 < ind & ind <= framo$nlevel]
			}

			fco = framo$fc
			dato = fco@base+fco@step
			if (dato != date) {
				cat("--> ref file has different time:",strftime(dato,"%H/%Y%m%d"),"\n")
				framo = NULL
			}

			prefixo = sprintf("ref%s",prefixp)
			prefixd = character()
			if (! mapstat) prefixd = sprintf("diff%s",prefixp)
		}

		lc = list()
		selev = length(frlow$ilev) < length(frlow$eta)
		for (j in seq(npar)) {
			if (isbin[i] || regexpr("^\\.",desc$faname[j]) > 0) next

			patt = FApattern(desc$ltype[j],desc$faname[j],frlow$ilev,selev)
			if (is.null(patt)) next

			ss = desc$symbol[j]

			cat(". get fields matching pattern",patt,"\n")
			stopifnot(is.null(ldata[[j]]))
			n = length(lc)+1
			lc[[n]] = mcparallel(getField(fics[i],patt,ss,frame,frlow))
			if (n == 4 || j == npar) {
				ld = mccollect(lc)
				try(parallel:::mckill(lc,15),silent=TRUE)
				for (k in 1:n) {
					if (is.null(ld[[k]])) stop("--> field not found in FA/GRIB file\n")
					ldata[[j-n+k]] = ld[[k]]
				}

				rm(lc,ld)
				lc = list()
			}
		}

		if (! is.null(framo)) {
			ldatao = lc = list()
			selev = length(framo$ilev) < length(framo$eta)
			for (j in seq(npar)) {
				if (isbin[i] || regexpr("^\\.",desc$faname[j]) > 0) next

				ss = desc$symbol[j]
				patto = FApattern(desc$ltype[j],desc$faname[j],framo$ilev,selev)
				#if (patto != patt) cat("--> new pattern for ref:",patto,"\n")

				cat(". get ref fields matching pattern",patt,"\n")
				if (length(ldatao) >= j) stopifnot(is.null(ldatao[[j]]))
				n = length(lc)+1
				lc[[n]] = mcparallel(getField(ficref,patto,ss,framo,frlow,cache.alt=TRUE,
					mc.cores=8))
				if (n == 4 || j == npar) {
					ld = mccollect(lc)
					try(parallel:::mckill(lc,15),silent=TRUE)
					for (k in 1:n) {
						if (is.null(ld[[k]])) stop("--> field not found in FA/GRIB file\n")
						ldatao[[j-n+k]] = ld[[k]]
					}

					rm(lc,ld)
					lc = list()
				}
			}
		}

		lsave = lsaveo = list()

		for (j in seq(npar)) {
			ss = desc$symbol[j]

			if (isbin[i] && regexpr("^\\.",desc$faname[j]) > 0) {
				cat(". get field as diag",desc$faname[j],"\n")
				if (is.null(ldata[[j]])) {
					cat("--> field not found in .bin file\n")
					next
				}

				f = setDiagData(frame$g4,frlow$eta,frlow$ilev,ldata[[j]])
			} else if (! isbin[i] && regexpr("^\\.",desc$faname[j]) < 0) {
				cat(". get field",desc$faname[j],"\n")
				f = ldata[[j]]
			} else {
				cat("--> no field",desc$faname[j],"in",fics[i],"(OK)\n")
				next
			}

			if (! is.null(cargs$test)) {
				browser()
				interptest(f)
				if (! is.null(framo)) interptesto(f,framo$g4)
			}

			nl = length(f@eta)
			if (desc$conv[j] != "") {
				cat("Convert data by:",desc$conv[j],"\n")
				fconv = eval(parse(text=desc$conv[j]))
				f[] = fconv(f)
			}

			cat("Primary checks\n")
			if (i == 1) {
				stopifnot(f@grid == frlow$g4)
				if (length(f@eta) > 1) stopifnot(all(f@eta == frlow$eta[frlow$ilev]))
				dataok = sapply(seq(along=f@eta),
					function(i) all(! is.infinite(f[,i])) && any(! is.na(f[,i])))
				stopifnot(all(dataok))
			}

			fnx = range(f,na.rm=TRUE)
			if (diff(fnx) == 0) {
				cat("--> constant value of ",fnx[1],"\n")
				next
			}

			cat("Compute stats and graphics for domains\n")
			tt = desc$longname[j]
			tt[2] = sprintf("on %s +%s",base,hh)
			qv = mapexi(f,desc[j,],doms,prefix=prefixp,main=tt,mnx=mapstat,mc.cores=4)
			if (is.null(lstat[[j]])) lstat[[j]] = statarray(qv,dates,dff$file)
			stopifnot(all(dim(qv) == dim(lstat[[j]])[1:3]))
			lstat[[j]][,,,id,i] = qv

			if (nzchar(desc$save[j])) {
				n = length(lsave)
				lsave[[n+1]] = f
				names(lsave)[n+1] = desc$save[j]
			}

			if (length(ldatao) < j || is.null(ldatao[[j]])) {
				rm(qv)
				gc()
				next
			}

			cat(". get ref field",desc$faname[j],"\n")
			if (regexpr("^\\.",desc$faname[j]) > 0) {
				if (! isbin[i]) next
				if (is.null(ldata[[j]])) next

				fo = setDiagData(framo$g4,frlow$eta,frlow$ilev,ldatao[[j]])
			} else {
				if (isbin[i]) next

				fo = ldatao[[j]]
			}

			if (desc$conv[j] != "") {
				cat("Convert data by:",desc$conv[j],"\n")
				fo[] = fconv(fo)
			}

			cat("Primary checks\n")
			if (i == 1) {
				stopifnot(fo@grid == frlow$g4)
				if (length(fo@eta) > 1) stopifnot(all(fo@eta == frlow$eta[frlow$ilev]))
				dataok = sapply(seq(along=fo@eta),
					function(i) all(! is.infinite(fo[,i])) && any(! is.na(fo[,i])))
				stopifnot(all(dataok))
				stopifnot(all(is.na(f) == is.na(fo)))
			}

			if (nzchar(desc$save[j])) {
				n = length(lsaveo)
				lsaveo[[n+1]] = f
				names(lsaveo)[n+1] = desc$save[j]
			}

			fd = setDataPart(f,f-fo)

			cat("Compute stats and graphics for domains\n")
			tto = c(tt[1],strftime(date,"valid: %Y%m%d %Hh"))
			tt[1] = sprintf("Diff of %s",tt[1])
			mapexi(fo,desc[j,],doms,0:1,prefix=prefixo,main=tto,mnx=mapstat,mc.cores=4)
			qv = mapexi(fd,descb[j,],doms,prefix=prefixd,main=tt,mnx=FALSE,mc.cores=4)

			cat("Compute zonal mean\n")
			t1 = mcparallel(zonalmean(f))
			t2 = mcparallel(zonalmean(fo))

			lc = mccollect(list(t1,t2))
			try(parallel:::mckill(list(t1,t2),15),silent=TRUE)
			zm = lc[[1]]
			zo = lc[[2]]
			rm(lc)

			if (is.null(lstatd[[j]])) {
				lstatd[[j]] = statarray(qv,dates,dff$file)
				lzmeand[[j]] = statarray(zm,dates,dff$file)
				lerr[[j]] = statarray(fd,c("rmse","bias","errx","dayx"),dff$file)
			}

			stopifnot(all(dim(qv) == dim(lstatd[[j]])[1:3]))
			stopifnot(all(dim(zm) == dim(lzmeand[[j]])[1:2]))
			stopifnot(all(dim(zo) == dim(lzmeand[[j]])[1:2]))
			lstatd[[j]][,,,id,i] = qv
			lzmeand[[j]][,,id,i] = zm-zo
			stopifnot(! is.null(dimnames(lstatd[[j]])))
			rm(fo,qv,zm,zo)

			if (mapstat) {
				cat("Compute errors\n")
				stopifnot(all(dim(fd) == dim(lerr[[j]])[1:2]))
				lerr[[j]][,,,i] = addrmxd(lerr[[j]][,,,i],fd,nerr[j,i],id,parallel=id>1)
				stopifnot(! is.null(dimnames(lerr[[j]])))
				nerr[j,i] = nerr[j,i]+1
			}

			rm(fd)
			gc()
			first = id == 1 && i == 1 && j == 1
			write.table(gc(),"gc.txt",append=!first,quote=FALSE,col.names=first)
		}

		for (k in seq(along=compound)) {
			cp = compound[[k]]
			if (! all(cp$vars %in% names(lsave))) next

			cat(". compound field:",names(compound)[k],"\n")
			fcomp = setDataPart(f,cp$formula(lsave[cp$vars]))

			j2 = match(cp$symbol,desc2$symbol)
			tt = desc2$longname[j2]
			tt[2] = sprintf("on %s +%s",as.character(fc@base,"%Y%m%d R%H"),hh)
			qv = mapexi(fcomp,desc2[j2,],doms,prefix=prefixp,main=tt)
			if (length(lstat2) < j2 || is.null(lstat2[[j2]])) {
				lstat2[[j2]] = statarray(qv,dates,dff$file)
				lerr[[j]] = statarray(fd,c("rmse","bias","errx","dayx"),dff$file)
			}

			lstat2[[j2]][,,,id,i] = qv

			if (! all(cp$vars %in% names(lsaveo))) next

			fcompo = setDataPart(f,cp$formula(lsaveo[cp$vars]))

			fd = setDataPart(fcomp,fcomp-fcompo)

			tto[1] = tt[1]
			tt[1] = sprintf("Diff of %s",tt[1])
			mapexi(fcompo,desc2[j2,],doms,0:1,prefix=prefixo,main=tto,mnx=mapstat,mc.cores=4)
			qv = mapexi(fd,desc2b[j2,],doms,prefix=prefixd,main=tt,mnx=FALSE,mc.cores=4)
			if (length(lstat2d) < j2 || is.null(lstat2d[[j2]])) {
				lstat2d[[j2]] = statarray(qv,dates,dff$file)
			}

			lstat2d[[j2]][,,,id,i] = qv

			if (mapstat) {
				lerr2[[j2]][,,,i] = addrmxd(lerr2[[j2]][,,,i],fd,nerr2[j2,i],id,parallel=id>1)
				nerr2[j2,i] = nerr2[j2,i]+1
			}

			rm(qv)
			gc()
		}

		gc()
	}
}

for (j in seq(along=lstat)) stopifnot(! is.null(dimnames(lstat[[j]])))

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
		par(c(Gpart,list(mfcol=c(nr,nc))))

		for (id in 1:min(ndom-nc*i,nc)+nc*i) {
			tt[2] = sprintf("domain %s",names(doms)[id])
			il0 = 1

			for (il in indl) {
				if (nl > 1) {
					tt[2] = sprintf("domain %s, lev. [%d,%d]",names(doms)[id],frlow$ilev[il0],
						frlow$ilev[il])
				}
				qfct = matrix(qfc[,il0:il,id,],ncol=length(ht))
				plotb(ht,qfct,main=tt,xlab=tlab,ylab=ss,xaxt="n")
				axis(1,ht6)
				il0 = min(il+1,max(indl))
			}
		}

		dev.off()
	}

	if (nl == 1) next

	indt = which(apply(qfc,4,function(x) any(! is.na(x))))
	if (length(indt) > 3) indt = indt[c(1,length(indt)%/%2,length(indt))]
	if (etahigh > 0) {
		indl = which(etai < etahigh)
		if (length(indl) < 3) indl = which(etai < .5)
		if (length(indl) < 3) stop(sprintf("not enough levels above .5"))
	}

	nc = length(indt)
	if (etahigh > 0) {
		nr = 2
		nj = 1
		tth = tt
		tth[1] = sprintf("..., highest levels",nom)
	} else if (wide) {
		nr = nj = 1
	} else {
		nr = nj = min(2,ndom)
	}

	for (i in seq((ndom-1)%/%nj+1)-1) {
		png(sprintf("%s/statv%d_%s.png",pngd,i,ss))
		par(c(Gpart,list(mfrow=c(nr,nc))))

		for (id in 1:min(ndom-nj*i,nj)+nj*i) {
			if (etahigh > 0) {
				for (it in indt) {
					tth[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],tunit)
					plotbv(qfc[,indl,id,it],etai[indl],main=tth,xlab=ss,ylab="eta",
						ylim=c(etahigh,yeta[2]))
				}
			}

			for (it in indt) {
				tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],ht[it],tunit)
				plotbv(qfc[,,id,it],etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
			}

		}

		dev.off()
	}
}

if (any(sapply(lstatd,is.array))) {
	for (j in seq(along=lstatd)) {
		stopifnot(! is.null(dimnames(lstatd[[j]])))
		stopifnot(! is.null(dimnames(lzmeand[[j]])))
	}

	i2 = sapply(lstat2d,is.array)
	if (any(i2)) {
		cat("Add lstatd2 to the statistics of errors\n")
		lstatd = c(lstatd,lstat2d[i2])
		lerr = c(lerr,lerr2[i2])
		nerr = rbind(nerr,nerr2[i2,])
	}

	cat("Save statistics\n")
	saveStat("diag.RData",frlow,params,names(doms),dates,htime,lstatd,lzmeand)

	if (mapstat) {
		cat("Maps and statistics of error (mean bias, RMSE, max and dayx error)\n")
		for (j in seq(along=lstatd)) stopifnot(! is.null(dimnames(lerr[[j]])))
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
				tt[2] = sprintf("step +%s",formatStep(fc@step,"h"))

				prefix = character()
				if (dff$graph[i]) prefix = formatStep(fc@step)

				cat(". bias, RMSE, errmax and index of errmax for",dff$file[i],"\n")
				bias = array(lerr[[j]][,,2,i]/nerr[j,i],dim(lerr[[j]])[1:2])
				fe = setDataPart(fe,bias)
				mapexi(fe,descb[j,],doms,0:1,prefix=sprintf("bias%s",prefix),main=tt,mnx=FALSE)

				rmse = array(sqrt(lerr[[j]][,,1,i]/nerr[j,i]),dim(lerr[[j]])[1:2])
				fe = setDataPart(fe,rmse)
				mapexi(fe,desc[j,],doms,0:1,prefix=sprintf("rmse%s",prefix),main=tt,mnx=FALSE)

				errx = array(lerr[[j]][,,3,i],dim(lerr[[j]])[1:2])
				fe = setDataPart(fe,errx)
				mapexi(fe,desc[j,],doms,0:1,prefix=sprintf("errx%s",prefix),main=tt,mnx=FALSE)

				dayx = array(lerr[[j]][,,4,i],dim(lerr[[j]])[1:2])
				fe = setDataPart(fe,dayx)
				mapexi(fe,desc[j,],doms,0:1,prefix=sprintf("dayx%s",prefix),main=tt,mnx=FALSE)
			}
		}
	}
}
