Gdiag = "~/util/diag"

library(maps)
library(mapproj)
library(parallel)

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

	doms = lapply(seq(dim(df)[1]),
		function(i) list(xlim=c(df$west[i],df$east[i]),ylim=c(df$south[i],df$north[i])))

	for (i in seq(dim(df)[1])) {
		doms[[i]]$proj = df$proj[i]
		if (df$param[i] != "-") doms[[i]]$param = df$param[i]
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

		grid = list(nfp=nfp,mucen=gem[1],locen=gem[2],rstret=gem[3],nlat=ndglg,nlong=nloeng,
			npdg=sum(nloeng),gauss=TRUE,lam=FALSE)
	} else if (type[2] == 1) {
		ndglg = type[3]
		ndlon = type[4]
		grid = list(nfp=nfp,nlat=ndglg,nlong=rep(ndlon,ndglg),npdg=ndglg*ndlon,
			gauss=FALSE,lam=TRUE)
	} else {
		stop("unknown grid type")
	}

	nx = readBin(con,"integer",1,endian="swap")
	nl1 = nx/nfp
	Ah = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))
	nx = readBin(con,"integer",1,endian="swap")
	Bh = readBin(con,"numeric",nl1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	grid$eta = ((Ah[-1]+Ah[-nl1])/Gvp0+(Bh[-1]+Bh[-nl1]))/2
	grid$nlevel = nl1-1

	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/4 == 2)
	base = readBin(con,"integer",2,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	nx = readBin(con,"integer",1,endian="swap")
	stopifnot(nx/nfp == 1)
	grid$step = readBin(con,"numeric",1,endian="swap")
	stopifnot(nx == readBin(con,"integer",1,endian="swap"))

	grid$base = as.POSIXct(as.Date(as.character(base[1]),"%Y%m%d"))+base[2]

	grid
}

is.equal = function(grid,frame)
{
	if (grid$npdg != frame$npdg || grid$nlevel != frame$nlevel ||
		grid$nlat != frame$nlat) return(FALSE)
	if (! identical(grid$nlong,frame$nlong)) return(FALSE)

	if (grid$gauss) {
		if (! isTRUE(all.equal(grid$mucen,sin(pi/180*latPole(frame)))) ||
			! isTRUE(all.equal(grid$locen,pi/180*lonPole(frame)))) return(FALSE)
	}

	grid$gauss == frame$gauss
}

toFrame = function(grid,frame)
{
	grid = c(grid,frame[c("lat","long","base")])
	if (grid$gauss) grid = c(grid,frame[c("theta","mucen","locen","rstret")])

	grid
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
		stopifnot(nx/grid$nfp == nl*grid$npdg)

		i = match(sub("^\\.?","",varname),cvars)
		if (is.na(i)) {
			seek(con,nl*grid$npdg*grid$nfp+4,"current")
			next
		}

		lvar[[i]] = readBin(con,"numeric",nl*grid$npdg,endian="swap")
		stopifnot(nx == readBin(con,"integer",1,endian="swap"))
		dim(lvar[[i]]) = c(grid$npdg,nl)
		names(lvar)[i] = varname
		if (all(! sapply(lvar,is.null))) break
	}

	lvar
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
		ech = as.integer(gsub(".+\\+0*([0-9]+)","\\1",file))
		if (step != ech*3600) step = ech*3600
	}

	lats = readBin(con,"numeric",ngp)
	longs = readBin(con,"numeric",ngp)

	Ah = readBin(con,"numeric",nl+1)
	Bh = readBin(con,"numeric",nl+1)
	eta = ((Ah[-1]+Ah[-(nl+1)])/Gvp0+(Bh[-1]+Bh[-(nl+1)]))/2

	if (gauss) gem = readBin(con,"numeric",3)

	close(con)

	if (gauss) {
		nlat = length(nlong)
		theta = 90-180*(seq(nlat)-.5)/nlat

		frame = list(nlat=nlat,nwave=dims[2],nlevel=nl,npdg=ngp,nlong=nlong,theta=theta,
			lat=lats,long=longs,eta=eta,base=base,step=step,mucen=gem[1],locen=gem[2],
			rstret=gem[3],gauss=TRUE,lam=FALSE)
	} else {
		nlong = rep(dims[2],dims[1])
		frame = list(nlat=dims[1],nlong=nlong,nwavex=dims[3],nwavey=dims[4],nlevel=nl,
			npdg=ngp,lat=lats,long=longs,eta=eta,base=base,step=step,gauss=FALSE,lam=TRUE)
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

getField = function(fic,param,symbol,frame,frlow=frame)
{
	stopifnot(! is.null(frlow$ilev))

	fsave = sprintf("%s/e%s/%s.RData",dirname(fic),formatStep(frame$step),symbol)

	if (file.exists(fsave) && file.info(fsave)$mtime < file.info(fic)$mtime) {
		cat("--> cache file older than data file, removed\n")
		file.remove(fsave)
	}

	if (file.exists(fsave)) {
		ilev = 0
		fic.save = ""
		load(fsave)

		if (! identical(ilev,frlow$ilev) || fic.save != fic || dim(data)[1] != frlow$npdg) {
			cat("--> cache file does not fit, change for 2nd cache file\n")
			fsave = sprintf("%s/e%s/%s.2.RData",dirname(fic),formatStep(frame$step),symbol)
			if (file.exists(fsave) && file.info(fsave)$mtime < file.info(fic)$mtime) {
				cat("--> 2nd cache file older than data file, removed\n")
				file.remove(fsave)
			}
		}
	}

	if (file.exists(fsave)) {
		load(fsave)

		if (fic.save != fic) {
			cat("--> different origin file",fic.save," read data again\n")
		} else if (! identical(ilev,frlow$ilev)) {
			cat("--> different levels, read data again\n")
		} else if (dim(data)[1] != frlow$npdg) {
			cat("--> different grid, read data again\n")
		} else {
			cat("--> cache file read:",fsave,"\n")
			return(data)
		}
	}

	data = getGPFields(fic,param,frame)

	h = hist(data,plot=FALSE)
	if (length(which(h$counts > 0)) == 2) {
		stopifnot(length(table(data)) == 2)
		meth = "ppp"
	} else {
		meth = "lin"
	}

	if (frame$npdg != frlow$npdg) {
		if (! "ind" %in% names(frlow) || frlow$npdg0 != frame$npdg) {
			cat("--> interpolation from",frame$npdg,"gridpoints - method:",meth,"\n")
			data = interpGauss(data,frame,frlow,meth)
		}

		if ("ind" %in% names(frlow)) {
			cat("--> selection of gridpoints:",frlow$npdg,"/",frame$npdg,"\n")
			data = data[frlow$ind,,drop=FALSE]
		}
	}

	ilev = frlow$ilev

	if (dim(data)[2] > 1 && dim(data)[2] != length(ilev)) {
		if (frame$nlevel != frlow$nlevel) {
			cat("--> interpolation from",frame$nlevel,"levels - method:",meth,"\n")
			if (meth == "ppp") {
				data = interpAB0(data,frame,frlow)
			} else {
				data = interpAB(data,frame,frlow)
			}
		}
	}

	fic.save = fic
	if (! file.exists(dirname(fsave))) dir.create(dirname(fsave))
	save("data","ilev","base","step","fic.save",file=fsave)

	data
}

getGPFields = function(file,field,frame)
{
	out = system(sprintf("epy_dump.py %s -f '%s' -o %s",file,field,Gficbin),
		intern=TRUE)
	out = grep("date|time|step|gridpoint size",out,ignore.case=TRUE,value=TRUE)
	cat(paste(out,collapse=", "),"\n")

	con = file(Gficbin,"rb")
	nl = readBin(con,"integer",1,size=8)
	stopifnot(length(nl) == 1)

	data = matrix(nrow=frame$npdg,ncol=nl,dimnames=list(NULL,seq(nl)))
	for (j in seq(nl)) data[,j] = readBin(con,"numeric",frame$npdg)
	stopifnot(nl == readBin(con,"integer",1,size=8))
	close(con)

	data
}

getSPFields = function(file,field,frame)
{
	system(sprintf("epy_dump.py %s -f %s -o %s",file,field,Gficbin),ignore.stdout=TRUE)
	con = file(Gficbin,"rb")
	nsp2 = readBin(con,"integer",1,size=8)
	data = readBin(con,"numeric",nsp2)
	close(con)

	stopifnot(nsp2 == frame$nwave*(frame$nwave-1))

	matrix(data,nrow=frame$nwave)
}

degrade = function(frame,ilat,nlat,nlon=2)
{
	stopifnot(all(frame$nlong%%nlon == 0))
	frame$npdg0 = frame$npdg

	if (frame$gauss) {
		if (missing(ilat)) {
			stopifnot(frame$nlat%%2 == 0)

			ilatn = seq(1,frame$nlat/2,nlat)
			ilat = c(ilatn,rev(frame$nlat-ilatn+1))
		}

		stopifnot(all(ilat %in% seq(frame$nlat)) && identical(ilat,sort(ilat)))

		clats = c(0,cumsum(frame$nlong))
		lats = datatoGauss(frame$lat,frame)
		longs = datatoGauss(frame$long,frame)

		# save
		nlong = frame$nlong

		# update
		frame$nlat = length(ilat)
		frame$nlong = frame$nlong[ilat]%/%nlon
		frame$npdg = sum(frame$nlong)
		frame$theta = frame$theta[ilat]

		frame$lat = frame$long = numeric(frame$npdg)
		frame$ind = integer(frame$npdg)

		ip = 0
		for (i in ilat) {
			ind = seq(1,nlong[i],nlon)
			np = length(ind)
			frame$lat[ip+1:np] = lats[ind,i]
			frame$long[ip+1:np] = longs[ind,i]
			frame$ind[ip+1:np] = clats[i]+ind
			ip = ip+np
		}
	} else {
		stopifnot(all(duplicated(frame$nlong)[-1]))
		ndlon = frame$nlong[1]

		ilat = seq(1,frame$nlat,nlat)
		ilon = seq(1,ndlon,nlon)
		frame$ind = rep((ilat-1)*ndlon,each=length(ilon))+ilon
		frame$lat = frame$lat[frame$ind]
		frame$long = frame$long[frame$ind]
		frame$nlat = length(ilat)
		frame$nlong = rep(length(ilon),frame$nlat)
		frame$npdg = length(frame$ind)
	}

	frame
}

csLat = function(frame)
{
	stop("not working")

	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	lat = numeric(frame$nlat)

	ind = frame$lat[cumsum(frame$nlong)-frame$nlong/2]
	cslon = pi

	for (ilat in seq(frame$nlat)) {
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		lat[ilat] = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
	}

	lat
}

csLat = function(frame)
{
	d1 = abs(diff(frame$long[frame$nlong[1]/2+0:1]))
	ioff = frame$npdg-frame$nlong[frame$nlat]/2
	d2 = abs(diff(frame$long[ioff+0:1]))
	C = round(sqrt(d2/d1),1)

	mu = sin(frame$theta*pi/180)
	sinla = ((C^2+1)*mu+C^2-1)/((C^2-1)*mu+C^2+1)
	theta2 = asin(sinla)*180/pi
	#Clat = seq(C,1/C,length.out=frame$nlat)

	#schmidt = Clat*tan(acos(sin(frame$theta*pi/180))/2)
	#theta2 = asin(cos(2*atan(schmidt)))*180/pi
	stopifnot(all(theta2 > frame$theta))

	theta2
}

csLong = function(frame)
{
	longs = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (i in seq(frame$nlat)) {
		nl = frame$nlong[i]
		longs[seq(nl),i] = 360*(seq(nl)-1)/nl
	}

	longs
}

datatoGauss = function(data,frame)
{
	clats = c(0,cumsum(frame$nlong))
	datag = matrix(nrow=max(frame$nlong),ncol=frame$nlat)

	for (ilat in seq(frame$nlat)) {
		ip = clats[ilat]
		np = frame$nlong[ilat]
		datag[1:np,ilat] = data[ip+1:np]
	}

	datag
}

latPole = function(frame)
{
	n1 = frame$nlong[1]

	if (all(duplicated(frame$lat[1:n1])[-1])) {
		latcen = 90
	} else if (diff(range(frame$long[1:n1])) < 180) {
		# lat of Pole is so that Pole is out of 1st lat circle
		latcen = mean(frame$lat[c(1,n1/2+1)])
	} else {
		# lat of Pole is close to geo North Pole and is inside 1st lat circle
		# latcen = (latx+2*(90-latx)+latn)/2 = 90+(latn-latx)/2 = 90-(latx-latn)/2
		latcen = 90-abs(diff(frame$lat[c(1,n1/2+1)]))/2
	}

	latcen
}

lonPole = function(frame)
{
	n1 = frame$nlong[1]
	frame$long[n1/2+1]
}

compass = function(frame)
{
	mucen = sin(pi/180*latPole(frame))
	locen = pi/180*lonPole(frame)
	C = 2.4
	c2 = C^2

	sqm2 = sqrt(1-mucen^2)
	sinlocen = sin(locen)
	coslocen = cos(locen)

	l = m = numeric(frame$npdg)

	off = 0
	for (ilat in seq(frame$nlat)) {
		nlon = frame$nlong[ilat]
		cslon = 2*pi*seq(0,nlon-1)/nlon
		mu = sin(pi/180*frame$theta[ilat])
		sqmu2 = sqrt(1-mu^2)
		a = 1/(c2+1+(c2-1)*mu)
		b = c2-1+(c2+1)*mu
		ca1 = a*(2*C*mucen*sqmu2-b*sqm2*cos(cslon))
		gemu = a*(2*C*sqm2*sqmu2*cos(cslon)+b*mucen)
		rcoslat = 1/sqrt(1-gemu^2)
		l[off+1:nlon] = -sqm2*sin(cslon)*rcoslat
		m[off+1:nlon] = ca1*rcoslat
		off = off+nlon
	}

	data.frame(l=l,m=m)
}

toU = function(nord,x,y)
{
	x*nord$m+y*nord$l
}

toV = function(nord,x,y)
{
	y*nord$m-x*nord$l
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

interpGauss = function(datao,framo,frame,meth)
{
	stopifnot(frame$gauss)

	# only for (Gauss) grids nested inside
	stopifnot(length(framo$nlong) >= length(frame$nlong))
	stopifnot(max(framo$nlong) >= max(frame$nlong))

	if (meth == "ppp") interp = interp0

	longs = csLong(frame)
	stopifnot(all(longs[1,] == 0))
	theta2 = csLat(frame)
	thetao2 = csLat(framo)

	data = matrix(nrow=frame$npdg,ncol=dim(datao)[2])
	dimnames(data)[[2]] = dimnames(datao)[[2]]

	clats = c(0,cumsum(frame$nlong))
	clato = c(0,cumsum(framo$nlong))
	dlat = 180/framo$nlat
	dlon = 360/framo$nlong

	for (ilat in seq(frame$nlat)) {
		ilato = max(which(thetao2 >= theta2[ilat]))
		stopifnot(0 < ilato && ilato < length(framo$nlong))

		#e0 = (90-dlat/2-frame$theta[ilat])/dlat
		#ilato = floor(e0)+1
		#e0 = e0-(ilato-1)
		e0 = (thetao2[ilato]-theta2[ilat])/(thetao2[ilato]-thetao2[ilato+1])
		stopifnot(0 <= e0 && e0 < 1)

		ilon = seq(frame$nlong[ilat])
		e1 = longs[ilon,ilat]/dlon[ilato]
		e2 = longs[ilon,ilat]/dlon[ilato+1]
		ilono1 = floor(e1)+1
		ilono2 = floor(e2)+1
		stopifnot(all(0 < ilono1 & ilono1 <= framo$nlong[ilato]))
		stopifnot(all(0 < ilono2 & ilono2 <= framo$nlong[ilato+1]))
		e1 = e1-(ilono1-1)
		e2 = e2-(ilono2-1)
		stopifnot(all(0 <= e1 & e1 < 1))
		stopifnot(all(0 <= e2 & e2 < 1))
		ip = clats[ilat]+ilon
		data[ip,] = interp(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
	}

	data
}

interp0 = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	if (e0 <= .5) {
		ip1 = clato[ilato]+ilono1
		ip = ifelse(e1 <= .5,ip1,ip1%%clato[ilato+1]+1)
	} else {
		ip1 = clato[ilato+1]+ilono2
		ip = ifelse(e2 <= .5,ip1,ip1%%clato[ilato+2]+1)
	}

	datao[ip,]
}

interp = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	ip1 = clato[ilato]+ilono1
	ip2 = ip1%%clato[ilato+1]+1
	d1 = datao[ip1,]+e1*(datao[ip2,]-datao[ip1,])
	ip1 = clato[ilato+1]+ilono2
	ip2 = ip1%%clato[ilato+2]+1
	d2 = datao[ip1,]+e2*(datao[ip2,]-datao[ip1,])
	d1+e0*(d2-d1)
}

interpAB0 = function(datao,framo,frame)
{
	eta = frame$eta[frame$ilev]
	etao = framo$eta
	if (! is.null(framo$ilev)) etao = etao[framo$ilev]
	ind = findInterval(eta,etao)
	e = (eta-etao[ind])/(etao[ind+1]-etao[ind])
	stopifnot(all(0 < ind & ind < length(etao)))
	stopifnot(all(0 <= e & e < 1))

	ind = ifelse(e <= .5,ind,ind+1)
	datao[,ind,drop=FALSE]
}

interpAB = function(datao,framo,frame)
{
	eta = frame$eta[frame$ilev]
	etao = framo$eta
	if (! is.null(framo$ilev)) etao = etao[framo$ilev]
	ind = findInterval(eta,etao)
	e = (eta-etao[ind])/(etao[ind+1]-etao[ind])
	stopifnot(all(0 < ind & ind < length(etao)))
	stopifnot(all(0 <= e & e < 1))

	data = matrix(nrow=dim(datao)[1],ncol=length(eta))
	dimnames(data)[[2]] = eta

	for (i in seq(along=eta)) {
		d1 = datao[,ind[i]]
		if (e[i] < 1.e-3) {
			data[,i] = d1
			next
		}

		d2 = datao[,ind[i]+1]
		if (1-e[i] < 1.e-3) {
			data[,i] = d2
			next
		}

		data[,i] = d1+e[i]*(d2-d1)
	}

	data
}

dilat = function(frame)
{
	stopifnot(frame$gauss)

	clats = c(0,cumsum(frame$nlong))

	ddm = numeric(frame$nlat)

	for (i in seq(frame$nlat)) {
		ip = clats[i]+seq(frame$nlong[i])

		dlon = abs(diff(frame$long[ip[1:2]]))
		a = cos(frame$lat[ip[1]]*pi/180)
		if (frame$mucen < 1) {
			stopifnot(abs(diff(frame$lat[ip[1:2]])) < abs(diff(frame$lat[ip[3:4]])))
		}

		ddm[i] = sqrt(diff(frame$lat[ip[1:2]])^2+(a*pmin(dlon,360-dlon))^2)
	}

	ddm/max(ddm[1:(frame$nlat/2)])
}

lat45 = function(data,frame)
{
	stopifnot(frame$gauss)

	dlon = cos(frame$theta*pi/180)*360/frame$nlong
	gm = dlon[frame$nlat%/%2]/dlon

	clats = c(0,cumsum(frame$nlong))

	i = frame$nlat%/%4
	ip = clats[i]+seq(frame$nlong[i])
	data[ip,,drop=FALSE]
}

ddiffn = function(x,y,n,yprev)
{
	if (n == 0) {
		return(y)
	} else if (n == 1) {
		return(diff(y[1:2])/diff(x[1:2]))
	}

	if (missing(yprev)) {
		yx1 = ddiffn(x[-(n+1)],y[-(n+1)],n-1)
	} else {
		yx1 = yprev
	}

	yx2 = ddiffn(x[-1],y[-1],n-1)
	(yx2-yx1)/(x[n+1]-x[1])
}

udiffn = function(y,n,yprev)
{
	if (n == 0) {
		return(y)
	} else if (n == 1) {
		return(diff(y[1:2]))
	}

	if (missing(yprev)) {
		yx1 = udiffn(y[-(n+1)],n-1)
	} else {
		yx1 = yprev
	}

	yx2 = udiffn(y[-1],n-1)
	(yx2-yx1)/n
}

Pn = function(x,y,n,xh)
{
	yprev = y
	for (i in (1:n)+1) yprev[i] = ddiffn(x[1:i],y[1:i],i-1,yprev[i-1])

	yprev = yprev[-1]
	if (missing(xh)) return(yprev)

	y[1]+sum(sapply(1:n,function(i) yprev[i]*prod(xh-x[1:i])))
}

uPn = function(y,n,h)
{
	yprev = y
	for (i in (1:n)+1) yprev[i] = udiffn(y[1:i],i-1,yprev[i-1])

	yprev = yprev[-1]
	if (missing(h)) return(yprev)

	y[1]+sum(sapply(1:n,function(i) yprev[i]*prod(h-(1:i))))
}

dprodf = function(x,xh)
{
	if (length(x) == 1) return(1)

	sum(sapply(seq(along=x),function(j) prod(xh-x[-j])))
}

dPn = function(x,y,n,xh)
{
	newt = Pn(x,y,n)
	sum(sapply(1:n,function(i) newt[i]*dprodf(x[1:i],xh)))
}

duPn = function(y,n,h)
{
	newt = Pn(y,n)
	sum(sapply(1:n,function(i) newt[i]*dprodf(1:i,h)))
}

gradl = function(data,frame,method="newton")
{
	stopifnot(frame$gauss)
	stopifnot(method %in% c("newton","FD"))

	clats = c(0,cumsum(frame$nlong))

	datal = data

	for (i in seq(frame$nlat)) {
		ip = clats[i]+seq(frame$nlong[i])
		np = length(ip)

		ip1 = c(ip[np],ip[-np])
		ip3 = c(ip[-1],ip[1])
		ip0 = c(ip[np-1:0],ip[0:1-np])
		ip4 = c(ip[-(1:2)],ip[1:2])

		if (method == "newton") {
			m = t(matrix(c(ip0,ip1,ip,ip3,ip4),ncol=5))

			for (j in seq(dim(data)[2])) {
				a = mclapply(seq(dim(m)[2]),function(k) dPn(1:5,data[m[,k],j],4,3),mc.cores=8)
				datal[ip,j] = simplify2array(a)
			}
		} else {
			dlon = abs(diff(frame$long[ip[1:2]]))
			a = cos(frame$lat[ip[1]]*pi/180)
			dl = sqrt(diff(frame$lat[ip[1:2]])^2+(a*pmin(dlon,360-dlon))^2)
			d3 = (data[ip3,]-data[ip1,])/(2*dl)
			d5 = (data[ip3,]-data[ip1,])/(4*dl)
			datal[ip,] = 2*d3-d5
		}
		if (i%%10 == 0) cat(". lat",i,"\n")
	}

	datal
}

vgrad = function(data,frame)
{
	eta = frame$A/Gvp0+frame$B
	deta = diff(eta[frame$ilev])
	data = apply(data,1,function(x) diff(x)/deta)

	t(data)
}

vmean = function(data,m=median,decreasing=TRUE)
{
	np = dim(data)[1]
	nlev = dim(data)[2]

	indf = apply(data,1,order,decreasing=decreasing)
	mm1 = matrix(c(1:np,indf[,1]),ncol=2)
	ddm01 = data[mm1]
	if (nlev < 4) return(ddm01)

	mm2 = matrix(c(1:np,indf[,2]),ncol=2)
	ddm02 = data[mm2]

	dd = data
	ip = logical(np)

	for (im in seq(nlev%/%4)) {
		for (i in seq(nlev-2*im)+im) dd[,i] = apply(dd,1,function(x) m(x[c(i-im,i+im)]))
		#indf = apply(dd,1,order,decreasing=decreasing)
		ddm1 = dd[mm1]
		ddm2 = dd[mm2]

		ip = ip | abs(ddm1-ddm01) < abs(ddm2-ddm02)
		ddm1[ip] = ddm2[ip]
		if (all(ip)) break

		iip = ! ip
		ddm01[iip] = ddm1[iip]
		ddm02[iip] = ddm2[iip]
	}

	ddm1
}

makeDomain = function(frame,ind,dlat,dlon)
{
	stopifnot(length(ind) == 1)
	stopifnot(length(dlat) == 2 && length(dlon) == 2)
	stopifnot(diff(dlat) > 0 && diff(dlon) > 0)

	xlim = frame$long[ind]+dlon
	ylim = frame$lat[ind]+dlat
}

getDomain = function(frame)
{
	xlim = longRange(frame$long)
	ylim = range(frame$lat)

	list(xlim=xlim,ylim=ylim)
}

longRange = function(long)
{
	dlong1 = range(long%%360)
	dlong2 = range((long+180)%%360-180)

	# if long span more than 180, make no choice (take whole globe)
	if (diff(dlong1) > 180 && diff(dlong2) > 180) return(c(-180,180))
	#stop("longitudes span more than 180dg (case not supported)")

	if (diff(dlong1) <= diff(dlong2)) {
		return(dlong1)
	} else {
		return(dlong2)
	}
}

inDomain = function(points,domain)
{
	ilat = domain$ylim[1] <= points$lat & points$lat <= domain$ylim[2]

	# case xlim=[0,0] (whole globe) is in alternative
	xlim = (domain$xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= points$long & points$long <= xlim[2]
	} else {
		ilon = points$long <= xlim[2] | points$long >= xlim[1]
	}

	ilat & ilon
}

area = function(dom1,dom2)
{
	diff(dom1$xlim)*diff(dom1$ylim)/(diff(dom2$xlim)*diff(dom2$ylim))
}

select = function(frame,ind)
{
	list(lat=frame$lat[ind],long=frame$long[ind])
}

zoom = function(data,frame,domain)
{
	mask = inDomain(frame,domain)
	frame$ind = which(mask)

	frame$lam = frame$lam || diff(domain$xlim) < 360
	frame$lat = frame$lat[frame$ind]
	frame$long = frame$long[frame$ind]
	frame$data = data[frame$ind,,drop=FALSE]

	frame
}

section = function(data,frame,long,lat=c(-90,90))
{
	stopifnot(frame$gauss)
	stopifnot(all(-90 <= lat & lat <= 90))

	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	if (length(lat) == 1) {
		dlat = 180/frame$nlat
		ilat = which.max(frame$theta <= lat)
		e0 = (90-dlat/2-frame$theta[ilat])/dlat-ilat+1
		ind = c(ilat,ilat+1)
	} else {
		ind = which(min(lat) <= frame$theta & frame$theta <= max(lat))
	}

	data1 = matrix(nrow=length(ind),ncol=dim(data)[2])
	lats = longs = numeric(length(ind))

	# frame longs belong to [-180,180[
	long = (long+180)%%360-180

	for (i in seq(along=ind)) {
		ilat = ind[i]
		e = long/dlon[ilat]
		ilon = floor(e)+1
		e = e-(ilon-1)
		ip1 = clats[ilat]+ilon
		ip2 = ip1%%clats[ilat+1]+1
		data1[i,] = (1-e)*data[ip1,]+e*data[ip2,]
		lats[i] = (1-e)*frame$lat[ip1]+e*frame$lat[ip2]
		longs[i] = frame$long[ip1]
	}

	if (length(lat) == 1) data1 = data1[1,]+e0*(data1[2,]-data1[1,])

	list(lat=frame$theta[ind],lats=lats,longs=longs,data=data1)
}

sectiongeo = function(data,frame,long=c(0,360),lat=c(-90,90))
{
	stopifnot(all(-90 <= lat & lat <= 90))

	# frame longs belong to [-180,180[
	if (length(long) == 1) {
		stopifnot(diff(lat) > 0)
		xf = frame$lat
		yf = frame$long
		x = lat
		y = (long+180)%%360-180
	} else if (length(lat) == 1) {
		stopifnot(diff(long) > 0)
		xf = frame$long
		yf = frame$lat
		x = (long+180)%%360-180
		y = lat
	} else {
		stop("lat or long must be of length 1")
	}

	clats = c(0,cumsum(frame$nlong))

	data1 = array(dim=c(2,frame$nlat,dim(data)[2]))
	x1 = array(dim=c(2,frame$nlat))
	iy = integer(frame$nlat)
	my = numeric(frame$nlat)

	for (ilat in seq(frame$nlat)) {
		off = clats[ilat]
		nlon = frame$nlong[ilat]
		xfn = xf[off+1:nlon]
		yfn = yf[off+1:nlon]
		yi = y
		if (frame$gauss) {
			# 1 more point in longitude (1st one) since frame is a Gaussian global grid
			xfn = c(xfn,xfn[1])
			yfn = c(yfn,yfn[1])
			if (length(long) == 1) {
				il = which(abs(diff(yfn)) > 180)
				stopifnot(length(il) <= 2)
				for (i in rev(il)) yfn[-(1:i)] = yfn[-(1:i)]-360*sign(diff(yfn)[i])
				if (yi < min(yfn)) yi = yi+360
				stopifnot(all(diff(yfn) > -180))
			} else {
				il = which(abs(diff(xfn)) > 180)
				stopifnot(length(il) <= 2)
				for (i in rev(il)) xfn[-(1:i)] = xfn[-(1:i)]-360*sign(diff(xfn)[i])
				stopifnot(all(diff(xfn) > -180))
			}
		}

		#ind = yfn[-nlon] <= yi & yi < yfn[-1] | yfn[-1] <= yi & yi < yfn[-nlon]
		n = length(yfn)
		iy[ilat] = which.min(abs(yi-yfn))
		indy = iy[ilat]+seq(-2,2)
		indy = indy[0 < indy & indy <= length(yfn)]
		my[ilat] = mean(yfn[indy])
		ind = yfn[-n] <= yi & yi < yfn[-1] | yfn[-1] <= yi & yi < yfn[-n]
		i1 = which(ind)
		if (length(i1) == 0) next

		stopifnot(length(i1) <= 2)
		i2 = i1+1
		i2[i1 == n] = 1
		e = (yi-yfn[i1])/(yfn[i2]-yfn[i1])
		# e can be 1 because of precision:
		# e = (b*(1-eps)-a)/(b-a)=1-b*eps/(b-a)=1-eps/(1-a/b), and e=1 for some cases
		stopifnot(all(0 <= e & e <= 1))

		for (i in seq(along=i1)) {
			data1[i,ilat,] = (1-e[i])*data[off+i1[i],]+e[i]*data[off+i2[i],]
			x1[i,ilat] = (1-e[i])*xfn[i1[i]]+e[i]*xfn[i2[i]]
		}
	}

	x1 = as.vector(x1)
	if (length(na.omit(x1)) < 5) {
		stopifnot(length(lat) == 1)

		ilat = 1+(which(! is.na(x1))-1)%/%2
		cat("--> few lats crossing:",length(ilat),"\n")
		if (length(ilat) == 0) {
			ilat = which.min(my)
		} else if (length(ilat) >= 2) {
			ilat = ilat[which.min(my[ilat])]
		}

		ind = clats[ilat]+seq(frame$nlong[ilat])
		return(list(longs=xf[ind],lat=lat,data=data[ind,]))
	}

	if (length(na.omit(x1)) == 0) stop("no lat/long crossing given parameter")

	if (length(lat) == 1 && diff(x) < 0) {
		xn = x1[x1 < 0]
		xp = x1[x1 >= 0]
		indn = order(xn,na.last=NA)
		indp = order(xp,na.last=NA)
		ind = match(c(xp[indp],xn[indn]),x1)
		ii = which(x[1] <= x1[ind] | x1[ind] <= x[2])
		cat("--> long:",x1[ind],"\n")
	} else {
		ind = order(x1,na.last=NA)
		ii = which(x[1] <= x1[ind] & x1[ind] <= x[2])
	}

	ind = ind[ii]
	stopifnot(all(! is.na(x1[ind])))

	# account for points "beyond pole", ie long=long+180%%360
	if (length(long) == 1) {
		#if (min(lat) <= -90)... 
		#if (max(lat) >= 90)... 
	}

	dim(data1) = c(2*frame$nlat,dim(data)[2])
	stopifnot(all(! is.na(data1[ind,1])))

	if (length(long) == 1) {
		list(lats=x1[ind],long=long,data=data1[ind,])
	} else {
		list(longs=x1[ind],lat=lat,data=data1[ind,])
	}
}

sectionf = function(data,frame,lat,long=c(0,360))
{
	stopifnot(frame$gauss)

	long = unique(long)
	stopifnot(-90 <= lat && lat < 90)
	stopifnot(length(long) == 2)
	stopifnot(all(0 <= long & long <= 360))

	longs = csLong(frame)
	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	ilat = max(which(frame$theta >= lat))

	if (long[1] < long[2]) {
		ind = which(long[1] <= longs[,ilat] & longs[,ilat] <= long[2])
	} else {
		ind = c(which(long[1] <= longs[,ilat] & longs[,ilat] <= 360),
			which(0 <= longs[,ilat] & longs[,ilat] <= long[2]))
	}

	data1 = matrix(nrow=length(ind),ncol=dim(data)[2])

	for (j in seq(dim(data)[2])) {
		datag = datatoGauss(data[,j],frame)
		data1[,j] = datag[ind,ilat]
	}

	lats = longs = numeric(length(ind))

	e = (lat-frame$theta[ilat])/diff(frame$theta[ilat+0:1])

	ip = clats[ilat]+ind

	list(long=longs[ind,ilat],lats=frame$lat[ip],longs=frame$long[ip],data=data1)
}

prettyBreaks = function(x,breaks="Sturges",nmin=8,n=10,crop=FALSE,split=FALSE)
{
	h = hist(x,breaks,plot=FALSE)

	if (is.character(breaks)) {
		stopifnot(0 < nmin && nmin <= n)

		# case for binary values (e.g. 0/1)
		if (length(which(h$counts > 0)) == 2 && all(which(h$counts > 0) == c(1,length(h$counts)))) {
			h = hist(x,2,plot=FALSE)
		}

		if (length(which(h$counts > 0)) >= n) h = hist(x,n,plot=FALSE)

		if (n > 2 && length(h$breaks) > n && length(which(h$counts > 0)) > 2) {
			h = hist(x,n-1,plot=FALSE)
		}

		if (length(h$breaks) < nmin) {
			hh = hist(x,nmin,plot=FALSE)
			# yes, 0s can be fewer with more bins...
			if (length(which(hh$counts == 0)) <= length(which(h$counts == 0))) h = hh
		}
	}

	if (split && length(h$breaks) > 3) {
		indx = order(h$density,decreasing=TRUE)
		dy = diff(h$density[indx])
		ix = which.min(dy)
		if (length(which(h$counts > 0)) > 1 && ix < length(which(h$density > 0)) &&
			(r=h$density[indx[ix]]/h$density[indx[ix+1]]) > 3*ix) {
			# n is 3, 4 or 5 (ie 2, 3 or 4)
			nb = 1+min(4,as.integer(sqrt(r/(3*ix))+1))
			br = seq(h$breaks[indx[1]],h$breaks[indx[1]+1],length.out=nb)
			h = hist(x,unique(sort(c(h$breaks,br))),plot=FALSE)
		}
	}

	if (crop && length(h$breaks) > 3) {
		br = h$breaks
		dxn = diff(br[1:2])
		dxx = diff(br[-(1:(length(br)-2))])
		if (min(x) > br[2]-dxn/5 && max(x) < br[length(br)-1]+dxx/5) {
			h$breaks[1] = br[2]-dxn/5
			h$breaks[length(br)] = br[length(br)-1]+dxx/5
		} else if (min(x) > br[2]-dxn/2 && max(x) < br[length(br)-1]+dxx/2) {
			h$breaks[1] = br[2]-dxn/2
			h$breaks[length(br)] = br[length(br)-1]+dxx/2
		}
	}

	if (length(h$counts) > 1 && h$breaks[length(h$breaks)-1] == max(x,na.rm=TRUE)) {
		h = hist(x,h$breaks[-length(h$breaks)])
	}

	h
}

mapxy = function(dom,...)
{
	if (dom$proj == "-") {
		l = map(xlim=dom$xlim,ylim=dom$ylim,...)
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
		l = map(projection=dom$proj,parameters=dom$param,xlim=dom$xlim,ylim=dom$ylim,...)
	}

	box()

	l
}

mapdom = function(dom,points,data,main=NULL,breaks="Sturges",palette="YlOrRd",pch=20,
	cex=.6,ppi=72,quiet=FALSE,...)
{
	if (par("mar")[1] > 2.5 || par("mar")[4] < 2.5) stop(paste(par("mar"),collapse=" "))

	l = mapxy(dom,mar=par("mar"),new=TRUE)

	u = par("usr")
	f = par("fin")

	if (ppi > 144) stop("ppi > 144")

	if (cex < .1) {
		warning("cex < 0.1, reset to 0.1")
		cex = .1
	}

	npmax = as.integer(min(prod(f*ppi/(3*cex)),.Machine$integer.max))
	if (length(data) < npmax/100) {
		cat("--> very few points, magnify plotting symbol (x2)\n")
		cex = 2*cex
	} else if (length(data) > npmax) {
		cex = max(.2,round(cex*sqrt(npmax/length(data)),3))
		npmax = as.integer(min(prod(f*ppi/(3*cex)),.Machine$integer.max))
	}

	x = points$long
	y = points$lat

	h = prettyBreaks(data,breaks,crop=TRUE)
	b = cut(data,h$breaks)

	if (length(data) > npmax) {
		if (! quiet) {
			cat("--> reducing xy points from",length(data),"to",npmax,"and cex to",cex,"\n")
		}

		if (points$gauss && points$rstret > 1) {
			for (nbin in 2:20) {
				np = length(data)%/%nbin
				ind = unlist(sapply(1:nbin,function(k) (k-1)*np+seq(1,np,by=nbin+1-k)))
				if (length(ind) < npmax) break
			}
		} else {
			ind = round(seq(1,length(data),length.out=npmax))
		}

		b2 = cut(data[ind],h$breaks)

		ii = which(h$counts > 0 & table(b2) == 0)
		if (length(ii) > 0) {
			ind1 = which(b %in% levels(b)[ii])
			ind = c(ind,ind1)
			stopifnot(all(! duplicated(ind)))
			if (! quiet) {
				cat("--> selecting back",length(ind1),"lost points in",length(ii),"data bins\n")
			}
		}

		x = x[ind]
		y = y[ind]
		data = data[ind]
	}

	br = h$breaks

	ind = findInterval(data,br,rightmost.closed=TRUE)
	rev = regexpr("\\+$",palette) < 0
	cols = hcl.colors(length(br),sub("\\+$","",palette),rev=rev)

	tind = table(ind)

	if (dom$proj == "-") {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			points(x[ii],y[ii],col=cols[ind[ii]],pch=pch,cex=cex,...)
		}
	} else {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			mp = mapproject(x,y)
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
	map.grid(c(long,long,dom$ylim[1],dom$ylim[2]),nx=1,labels=FALSE,pretty=FALSE,...)

	# W->E cross section
	map.grid(c(dom$xlim[1],dom$xlim[2],lat,lat),ny=1,labels=FALSE,pretty=FALSE,...)
}

mapdom2 = function(dom,points,zx,zy,main=NULL,breaks="Sturges",colvec=TRUE,length=.05,
	angle=15,ppi=6,quiet=FALSE,...)
{
	l = mapxy(dom,mar=par("mar"),new=TRUE)

	f = par("fin")

	if (ppi > 96) stop("ppi > 96\n")

	npmax = as.integer(min(prod(f*ppi),.Machine$integer.max))

	if (colvec) {
		ffz = sqrt(zx^2+zy^2)
		br = prettyBreaks(ffz,crop=TRUE)$breaks
	}

	x = points$long
	y = points$lat

	if (length(zx) > npmax) {
		if (! quiet) cat("--> reducing xy2 plot from",length(zx),"to",npmax,"points\n")
		if (points$gauss && points$rstret > 1) {
			for (nbin in 2:20) {
				np = length(zx)%/%nbin
				ind = unlist(sapply(1:nbin,function(k) (k-1)*np+seq(1,np,by=nbin+1-k)))
				if (length(ind) < npmax) break
			}
		} else {
			ind = as.integer(seq(1,length(zx),length.out=npmax))
		}

		x = x[ind]
		y = y[ind]
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

		ic = findInterval(ffz,br,rightmost.closed=TRUE)
		palette = rainbow(length(br),start=2/3,end=0)
		cols = palette[ic]
	} else {
		ff1 = ffu*scale(ffz,center=FALSE)
		fx = zx/ffz*ff1
		fy = zy/ffz*ff1
		cols = rep(1,length(ff1))
	}

	x2 = x+fx
	y2 = y+fy

	ffi = sqrt((fx/ux*f[1])^2+(fy/uy*f[2])^2)

	ind = ffi > 2.e-3

	if (dom$proj == "-") {
		arrows(x,y,x2,y2,length,angle,col=cols[ind],...)
	} else {
		mp = mapproject(x[ind],y[ind])
		mp2 = mapproject(x2[ind],y2[ind])
		arrows(mp$x,mp$y,mp2$x,mp2$y,length,angle,col=cols[ind],...)
	}

	if (any(! ind)) {
		if (dom$proj == "-") {
			points(x[! ind],y[! ind],pch=1,col=1)
		} else {
			mp = mapproject(x[! ind],y[! ind])
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

hist.annot = function(data,nmin=15,n=20,split=TRUE,...)
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

prettydiff = function(x,...)
{
	lab = pretty(x,...)
	n = length(lab)
	i = n%/%8
	j = i
	if (n%/%2*2 == n && j > 0) j = j-1
	lab[-(n%/%2+seq(-j,i))]
}

revpos = function(x,at)
{
	# for any xi in [a,b], xi = a+(xi-a), xi' = b-(xi-a) = a+b-xi
	sum(range(x))-at
}

plotv0 = function(x,y,z,ylim=NULL,xaxt="s",yaxt="s",...)
{
	if (! is.null(ylim)) {
		indy = which(ylim[1] <= y & y <= ylim[2])
		y = y[indy]
		indy = rev(indy)
	} else {
		indy = rev(seq(along=y))
	}

	contour(x,y,z[,indy],xaxt="n",yaxt="n",...)

	if (xaxt != "n") {
		lab = pretty(x/30)*30
		axis(1,at=lab)
	}

	if (yaxt != "n") {
		lab = pretty(y)
		axis(2,revpos(y,lab),lab)
	}
}

plotv = function(x,y,z,ylim=rev(range(y,finite=TRUE)),xaxs="i",yaxs="i",palette="YlOrRd",
	long=FALSE,...)
{
	if (par("mar")[1] > 2.5 || par("mar")[4] < 2) stop("mar does not fit")

	br = prettyBreaks(z,crop=TRUE)$breaks
	dbr = diff(range(br))*.Machine$double.eps
	stopifnot(all(min(br)-dbr <= z & z <= max(br)+dbr))
	if (length(br) == 2) br = c(br,br[2]+diff(br))
	rev = regexpr("\\+$",palette) < 0
	cols = hcl.colors(length(br)-1,sub("\\+$","",palette),rev=rev)
	if (long) {
		i = which(diff(x) < 0)
		if  (length(i) > 0) x[-(1:i)] = x[-(1:i)]+360
	}

	image(x,y,z,ylim=ylim,col=cols,breaks=br,xaxs=xaxs,yaxs=yaxs,...)
	lev = sprintf("% .3g",br)
	DOPfill.legend(lev,col=cols)
}

plotf = function(y,z,ylim=NULL,plot.axes=TRUE,...)
{
	if (is.null(ylim)) {
		indy = rev(seq(along=y))
	} else {
		indy = which(ylim[1] <= y & y <= ylim[2])
		y = y[indy]
		indy = rev(indy)
	}

	filled.contour(y=y,z=z[,indy],plot.axes=plot.axes,...)

	if (plot.axes) {
		lab = pretty(y)
		axis(2,revpos(y,lab),lab)
	}
}

plotsectiongeo = function(data,etai,frame,dom,main,...)
{
	longm = mean(dom$xlim)
	latm = mean(dom$ylim)

	yz = sectiongeo(data,frame,longm,dom$ylim)
	stopifnot(all(is.finite(yz$data)))

	main[2] = sprintf("S->N section, longitude %g",round(longm,1))
	plotv(yz$lats,etai,yz$data,main=main,xlab="",ylab="",...)

	xz = sectiongeo(data,frame,dom$xlim,latm)
	stopifnot(all(is.finite(xz$data)))

	main[2] = sprintf("W->E section, latitude %g",round(latm,1))
	plotv(xz$longs,etai,xz$data,main=main,xlab="",ylab="",long=TRUE,...)
}

mapexi = function(data,frame,frmap,desc,doms,prefix=NA,np=4,...)
{
	nl = dim(data)[2]
	qv = array(dim=c(101,nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))
	tt = desc$longname

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (frame$gauss && diff(doms[[i]]$xlim) > 350 && diff(doms[[i]]$ylim) > 170) {
			# no domain, no zoom
			domin = doms[[i]]

			if (is.null(frmap)) {
				xymap = frame
				xymap$data = data
			} else {
				xymap = frmap
				xymap$data = data[frmap$ind,,drop=FALSE]
			}

			xymap$lam = FALSE
		} else {
			pts = select(frame,ind)
			domin = getDomain(pts)
			if (area(domin,doms[[i]]) < .1) next

			xymap = zoom(data,frame,doms[[i]])
		}

		stopifnot(all(is.finite(xymap$data)))
		qv[,,i] = apply(xymap$data,2,quantile,prob=seq(0,100)/100)

		if (is.na(prefix)) next

		nom = names(doms)[i]
		cat(".. domain:",nom,"- points/total:",dim(xymap$data)[1],dim(data)[1],"\n")

		pre = sub("[0-9]+$","",prefix)
		ficpng = sprintf("%s/map%s%s_%s.png",pngd,prefix,nom,desc$symbol)
		cat(sprintf("\t<option>%s</option>\n",ficpng),
			file=sprintf("%s/map%s%s_%s.html",pngd,pre,nom,desc$symbol),append=TRUE)

		pngalt(ficpng)

		if (nl == 1) {
			op = par(c(Gparm,list(cex=.8)))
			mapdom(doms[[i]],xymap,xymap$data[,1],main=tt,palette=desc$palette,...)
		} else {
			op = par(c(Gparm,list(mfrow=c(2,2),cex=.7)))

			# be careful of nl=2...
			indl = seq(1,nl,length.out=min(4,nl))
			if (length(Gindl) == 2) indl = c(1,Gindl,nl)
			for (il in indl) {
				ilev = frame$ilev[il]
				tt[2] = sprintf("level %d (eta %g)",ilev,round(frame$eta[ilev],2))
				mapdom(doms[[i]],xymap,xymap$data[,il],main=tt,palette=desc$palette,
					quiet=il!=indl[1],...)
				mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
			}
		}

		pngoff()

		ficpng = sprintf("%s/hist%s%s_%s.png",pngd,prefix,nom,desc$symbol)
		cat(sprintf("<option>%s</option>\n",ficpng),
			file=sprintf("%s/hist%s%s_%s.html",pngd,pre,nom,desc$symbol),append=TRUE)

		pngalt(ficpng)

		if (nl == 1) {
			op = par(c(Gparh,list(cex=.8)))
			hist.annot(xymap$data,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
		} else {
			op = par(c(Gparm,list(mfcol=c(2,2)),cex=.7))
			plotsectiongeo(data,etai,frame,doms[[i]],main=desc$longname,palette=desc$palette)

			par(mar=Gparh$mar)
			hist.annot(xymap$data,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
			plotbv(qv[,,i],etai,main=tt[1],xlab="",ylab="",ylim=yeta)
			abline(v=0,col="darkgrey",lty=2)
		}

		pngoff()

		if (nl < 4) next

		il = (nl+1)%/%2

		ficpng = sprintf("%s/mapn%s%s_%s.png",pngd,prefix,nom,desc$symbol)
		pngalt(ficpng)

		op = par(c(Gparm,list(mfrow=c(2,2),cex=.7)))
		datan = apply(xymap$data[,1:il],1,min)
		tt[2] = sprintf("min. of lev. [1,%d]",il)
		mapdom(doms[[i]],xymap,datan,main=tt,palette="Heat+",...)
		mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
		par(mar=Gparh$mar)
		hist.annot(datan,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

		par(Gparm)
		datan = apply(xymap$data[,-(1:il)],1,min)
		tt[2] = sprintf("min. of lev. [%d,%d]",il+1,nl)
		mapdom(doms[[i]],xymap,datan,main=tt,palette="Heat+",...)
		mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
		par(mar=Gparh$mar)
		stopifnot(all(par("mar") == Gparh$mar))
		hist.annot(datan,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
		pngoff()

		ficpng = sprintf("%s/mapx%s%s_%s.png",pngd,prefix,nom,desc$symbol)
		pngalt(ficpng)

		op = par(c(Gparm,list(mfrow=c(2,2),cex=.7)))
		datax = apply(xymap$data[,1:il],1,max)
		tt[2] = sprintf("max. of lev. [1,%d]",il)
		mapdom(doms[[i]],xymap,datax,main=tt,palette="Heat",...)
		mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
		par(mar=Gparh$mar)
		hist.annot(datax,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)

		par(Gparm)
		datax = apply(xymap$data[,-(1:il)],1,max)
		tt[2] = sprintf("max. of lev. [%d,%d]",il+1,nl)
		mapdom(doms[[i]],xymap,datax,main=tt,palette="Heat",...)
		mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
		par(mar=Gparh$mar)
		hist.annot(datax,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
		pngoff()
	}

	qv
}

mapex2 = function(u,v,frame,doms,prefix=NA,np=4,palette,...)
{
	nl = dim(u)[2]
	qv = array(dim=c(101,nl,length(doms)),dimnames=list(NULL,seq(nl),names(doms)))

	for (i in seq(along=doms)) {
		ind = inDomain(frame,doms[[i]])
		if (length(which(ind)) < np) next

		if (frame$gauss && diff(doms[[i]]$xlim) > 350 && diff(doms[[i]]$ylim) > 170) {
			# world: no domain, no zoom
			domin = c(doms[[i]],list(lam=FALSE))
			x = list(lat=frame$lat,long=frame$long,data=u)
			y = list(lat=frame$lat,long=frame$long,data=v)
		} else {
			pts = select(frame,ind)
			domin = getDomain(pts)
			if (area(domin,doms[[i]]) < .1) next

			x = zoom(u,frame,doms[[i]])
			y = zoom(v,frame,doms[[i]])
		}

		ffdom = sqrt(x$data^2+y$data^2)
		stopifnot(all(is.finite(ffdom)))
		qv[,,i] = apply(ffdom,2,quantile,prob=seq(0,100)/100)

		if (is.na(prefix)) next

		nom = names(doms)[i]
		tt = "u/v wind"
		cat(".. domain:",nom,"\n")

		pre = sub("[0-9]+$","",prefix)
		ficpng = sprintf("%s/map%s%s_ff.png",pngd,prefix,nom)
		cat(sprintf("\t<option>%s</option>\n",ficpng),
			file=sprintf("%s/map%s%s_ff.html",pngd,pre,nom),append=TRUE)
		pngalt(ficpng)
		if (nl == 1) {
			op = par(c(Gparm,list(cex=.8)))
			mapdom2(doms[[i]],x,x$data,y$data,main=tt,...)
		} else {
			op = par(c(Gparm,list(mfcol=c(2,2),cex=.7)))

			# if nl=2...
			indl = seq(1,nl,length.out=min(4,nl))
			if (length(Gindl) == 2) indl = c(1,Gindl,nl)
			for (il in indl) {
				ilev = frame$ilev[il]
				tt[2] = sprintf("level %d (eta %g)",ilev,round(frame$eta[ilev],2))
				mapdom2(doms[[i]],x,x$data[,il],y$data[,il],main=tt,...)
				mapsegments(domin,lat=mean(domin$ylim),long=mean(domin$xlim),col="darkgrey")
			}
		}

		pngoff()

		ficpng = sprintf("%s/hist%s%s_ff.png",pngd,prefix,nom)
		cat(sprintf("\t<option>%s</option>\n",ficpng),
			file=sprintf("%s/hist%s%s_ff.html",pngd,pre,nom),append=TRUE)
		pngalt(ficpng)
		tt = "wind speed"
		if (nl == 1) {
			op = par(c(Gparh,list(cex=.8)))
			hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab="wind speed")
		} else {
			ff = sqrt(u^2+v^2)
			op = par(c(Gparm,list(mfcol=c(2,2),cex=.7)))
			plotsectiongeo(ff,etai,frame,domin,main=desc$longname,palette=desc$palette)

			par(mar=Gparh$mar)
			hist.annot(ffdom,col="whitesmoke",main=tt[1],xlab=NULL,ylab=NULL)
			plotbv(qv[,,i],etai,main=tt[1],xlab="",ylab="",ylim=yeta)
			abline(v=0,col="darkgrey",lty=2)
		}

		pngoff()
	}

	qv
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
	c(sqrt(mean(x^2,...)),mean(x,...),sd(x,...))
}

plotb = function(x,y,...)
{
	plot(x,apply(y,2,median),type="l",ylim=range(y,na.rm=TRUE),...)
	w = diff(range(x))
	boxplot(y,range=0,col="bisque",lty=1,pars=list(boxwex=.03*w),add=TRUE,at=x,xaxt="n",
		yaxt="n")
}

plotbv = function(x,y,...)
{
	plot(apply(x,2,median),y,type="l",xlim=range(x,na.rm=TRUE),...)
	w = diff(range(y))
	boxplot(x,range=0,col="bisque",lty=1,pars=list(boxwex=.03*w),horizontal=TRUE,add=TRUE,
		at=y,xaxt="n",yaxt="n")
}

matplott = function(x,y,col=1,...)
{
	matplot(x,y,type="o",lty=1:3,col=col,pch="+",...)
	abline(h=0,col="darkgrey")
	legend("topleft",c("RMS","bias","std dev"),col=col,lty=1:3,pch="+",bg="transparent")
}

matplotv = function(x,y,col=1,...)
{
	matplot(x,y,type="o",lty=1:3,col=col,pch="-",...)
	abline(v=0,col="darkgrey")
	legend("topleft",c("RMS","bias","std dev"),col=col,lty=1:3,pch="+",bg="transparent")
}

pngalt = function(...)
{
	if (ask && ! is.null(dev.list())) invisible(readline("Press enter to continue"))

	if (! hasx11) {
		stopifnot(is.null(dev.list()))
		png(...)
	}
}

pngoff = function(op)
{
	if (! hasx11) {
		invisible(dev.off())
	} else if (! missing(op)) {
		par(op)
	}
}

saveStat = function(filename,framep,paramsp,datesp,htimep,lstatdp,lstat2dp)
{
	if (! file.exists(filename)) {
		frlow = framep
		params = paramsp
		dates = datesp
		htime = htimep
		lstatd = lstatdp
		lstat2d = lstat2dp
		save("frlow","params","dates","htime","lstatd","lstat2d",file=filename)
		return()
	}

	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lstatd","lstat2d") %in% noms))

	if (frlow$gauss == frame$gauss &&
		frlow$nlat == frame$nlat && identical(frlow$nlong,frame$nlong) &&
		frlow$nlevel == frame$nlevel && identical(frlow$ilev,framep$ilev) &&
		identical(dates,datesp) && identical(htime,htimep)) {
		params = unique(c(params,paramsp))

		indp = match(paramsp,params)
		ip = is.na(indp)

		if (any(! ip)) lstatd[na.omit(indp)] = lstatdp[which(! ip)]
		if (any(ip)) lstatd = c(lstatd,lstatdp[which(ip)])
	}

	save("frlow","params","dates","htime","lstatd","lstat2d",file=filename)
}

loadStat = function(filename,framep,paramsp,datesp,htimep)
{
	noms = load(filename)
	stopifnot(all(c("frlow","params","dates","htime","lstatd","lstat2d") %in% noms))

	if (frlow$gauss != frame$gauss) return(NULL)
	if (frlow$nlat != frame$nlat || ! identical(frlow$nlong,frame$nlong)) return(NULL)
	if (frlow$nlevel != frame$nlevel || ! identical(frlow$ilev,framep$ilev)) return(NULL)
	if (! identical(sort(dates),sort(datesp))) return(NULL)

	indp = match(paramsp,params)
	if (any(is.na(indp))) return(NULL)

	indt = match(htimep,htime)
	if (any(is.na(indt))) return(NULL)

	lstatd = lapply(lstatd[indp],"[",indt)

	indp = sapply(lstat2d,is.array)
	if (any(indp)) lstatd = c(lstatd,lapply(lstat2d[indp],function(x) x[,,,,indt]))

	lstatd
}

Gparm = list(mar=c(1.8,1.8,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)
Gparh = list(mar=c(1.8,1.8,2.5,1),mgp=c(1.5,.5,0),tcl=-.3)
Gpart = list(mar=c(2.1,2.1,2.5,1),mgp=c(1.5,.5,0),tcl=-.3)
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
	dff = read.table("config/file.txt",header=TRUE)
}

dates = NA

if (file.exists("config/date.txt")) {
	dates = scan("config/date.txt",quiet=TRUE)
	res = as.integer(dates[1])
	dates = dates[-1]
} else if ("date" %in% names(cargs)) {
	cdate = strsplit(cargs$date,split="/")[[1]]
	dates = as.integer(cdate[1])
	res = as.integer(cdate[2])
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

hasx11 = ! "png" %in% names(cargs) && capabilities("X11")
ask = ! "noask" %in% names(cargs) && hasx11 && interactive()

fhtml = dir(pngd,"(map|hist)\\w+_\\w+\\.html",full.names=TRUE)
if (length(fhtml) > 0) invisible(file.remove(fhtml))

dstats = dstatd = array(dim=c(101,length(dates),length(dff$file),npar),
	dimnames=list(NULL,dates,dff$file,params))
lstat2 = lstat2d = lerr2 = vector("list",dim(desc2)[1])
names(lstat2) = names(lstat2d) = desc2$symbol
lstat = lstatd = lerr = vector("list",npar)
htime = integer(length(dff$file))
tstep = rep(NA,length(dff$file))
frame.save = NULL
quiet = FALSE

if (interactive()) browser()

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
			ldata = getVars(con,grid,desc$faname,quiet)
			close(con)
			quiet = TRUE
			stopifnot(exists("frame",mode="list") && is.equal(grid,frame))
			frame = toFrame(grid,frame)
		} else {
			frame = getFrame(fics[i])
		}

		base = as.character(frame$base,"%Y%m%d")
		if (! is.na(dates[id]) && base != dates[id]) {
			message("dates mismatch: ",frame$base," != ",dates[id])
			stop("date mismatch")
		}

		cat("--> Gauss grid:",frame$gauss,"\n")
		s = sprintf("base/step: %s +%s",base,formatStep(frame$step),"h")
		s = paste(s,sprintf("file: %s - graph: %s",fics[i],dff$graph[i]),sep=" - ")

		if (is.na(tstep[i])) {
			cat("-->",s,"\n")
			tstep[i] = s
		} else {
			stopifnot(sub("base.+? \\+(.+?) -.+","",s) == sub("base.+? \\+(.+?) -.+","",tstep[i]))
		}

		gg = id == 1 && dff$graph[i]
		frmap = NULL

		if (gg && frame$gauss) {
			mapf = dilat(frame)
			ilateq = equalize(mapf,offset=.5)
			#if (length(ilateq) < frame$nlat) frmap = degrade(frame,ilateq)

			pngalt(sprintf("%s/stretch%d.png",pngd,i))
			C = sqrt(mapf[length(mapf)]/mapf[1])
			tt = c("Stretching coefficient",sprintf("C = %g",round(C,2)))
			plot(mapf,type="l",main=tt,xlab="Latitude index",ylab="Stretching coef.",log="y")
			#points(ilateq,mapf[ilateq],pch=20,cex=.5)
			grid(col="darkgrey")
			pngoff()
		}

		if (is.null(frame.save)) {
			frame.save = frame
		} else {
			stopifnot(identical(frame$nlong,frame.save$nlong))
			stopifnot(frame$nlevel == frame.save$nlevel)
			stopifnot(frame$gauss == frame.save$gauss)
		}

		if (frame$nlat > nlatmax) {
			cat("--> halving initial resolution:",frame$nlat,max(frame$nlong),"\n")
			frlow = degrade(frame,nlat=2,nlon=2)
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

		date = frame$base+frame$step
		cat("--> validity:",strftime(date,"%Y%m%d %H"),"\n")
		htime[i] = frame$step
		ficref = sprintf("%s/%s/%s",ref,strftime(date,"%H/%Y%m%d"),dff$ref[i])
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
				data = ldata[[j]]
				if (is.null(data)) {
					cat("--> field not found in .bin file\n")
					next
				}

				if (dim(data)[2] > 1 && dim(data)[2] != length(frlow$ilev)) {
					cat("--> change 'bin' levels:",dim(data)[2],"-->",length(frlow$ilev),"\n")
					data = data[,frlow$ilev,drop=FALSE]
				}
			} else if (! isbin[i] && regexpr("^\\.",desc$faname[j]) < 0) {
				cat(". get fields matching pattern",patt,"\n")
				data = getField(fics[i],patt,ss,frame,frlow)
				if (is.null(data)) {
					cat("--> field not found in FA/GRIB file\n")
					next
				}
			} else {
				cat("--> no field",desc$faname[j],"in",fics[i],"(OK)\n")
				next
			}

			nl = dim(data)[2]
			if (desc$conv[j] != "") {
				fconv = eval(parse(text=desc$conv[j]))
				data[] = fconv(data)
			}

			# primary checks
			stopifnot(all(! is.infinite(data)))
			stopifnot(all(apply(data,2,function(x) any(! is.na(x)))))

			#if (! all(apply(data,2,function(x) length(unique(na.omit(x))) > 1))) {
			if (length(unique(na.omit(as.vector(data)))) == 1) {
				cat("--> constant value of ",na.omit(as.vector(data))[1],"\n")
				next
			}

			dstats[,id,i,j] = quantile(data[,nl],seq(0,100)/100,names=FALSE)

			prefix = ifelse(gg,i,NA_character_)
			qv = mapexi(data,frlow,frmap,desc[j,],doms,prefix)
			if (is.null(lstat[[j]])) {
				lstat[[j]] = array(dim=c(dim(qv),length(dates),length(fics)))
			}
			lstat[[j]][,,,id,i] = qv

			if (desc$save[j] == "X") {
				datax = data
			} else if (desc$save[j] == "Y") {
				datay = data
			} else if (desc$save[j] == "Z") {
				dataz = data
			} else if (desc$save[j] == "L") {
				datal = data
			} else if (desc$save[j] == "M") {
				datam = data
			} else if (desc$save[j] == "GPL" && frame$gauss) {
				datagpl = gradl(data,frame)
			}

			if (! file.exists(ficref)) next

			cat(". get ref field",desc$faname[j],"\n")
			if (regexpr("^\\.",desc$faname[j]) > 0) {
				if (! isbin[i]) next
				datao = ldata[[j]]
				if (is.null(datao)) next

				# for some vars (like cty) can be nl+1
				if (dim(datao)[2] > 1 && dim(datao)[2] != length(frlow$ilev)) {
					datao = datao[,frlow$ilev,drop=FALSE]
				}
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

				datao = getField(ficref,patt,ss,framo,frlow)
			}

			if (desc$conv[j] != "") datao = fconv(datao)

			if (desc$save[j] == "X") {
				dataox = datao
			} else if (desc$save[j] == "Y") {
				dataoy = datao
			} else if (desc$save[j] == "Z") {
				dataoz = datao
			}

			ddiff = data-datao
			if (is.null(lerr[[j]])) {
				lerr[[j]] = array(dim=c(dim(ddiff),length(dates),length(fics)))
			}
			lerr[[j]][,,id,i] = ddiff

			dstatd[,id,i,j] = quantile(ddiff[,nl],seq(0,100)/100,names=FALSE)

			prefix = ifelse(gg,sprintf("diff%d",i),NA_character_)
			mapexi(datao,frlow,frmap,desc[j,],doms,prefix)
			qv = mapexi(ddiff,frlow,frmap,desc[j,],doms)
			if (length(lstatd[[j]]) == 0) {
				lstatd[[j]] = array(dim=c(dim(qv),length(dates),length(fics)))
			}
			lstatd[[j]][,,,id,i] = qv

			rm(data,datao)
			gc()
		}

		if (! is.null(datax) && ! is.null(datay)) {
			cat(". compound field wind speed\n")
			ff = sqrt(datax^2+datay^2)
			print(summary(as.vector(ff)))

			if (gg && frame$gauss && latPole(frame) < 90) {
				cat("rotating wind to geo North\n")
				nord = compass(frame)
				if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
				u = ucs(nord,datax,datay)
				v = vcs(nord,datax,datay)

				pole = list(xlim=lonPole(frame)+c(-6,6),ylim=latPole(frame)+c(-5,5),
					proj="stereo")
				x = zoom(u,frlow,pole)
				y = zoom(v,frlow,pole)
				pngalt(sprintf("%s/map%dpole_ff.png",pngd,i))
				mapdom2(pole,x,x$data[,1],y$data[,1],main="u/v wind components")
				pngoff()
			} else {
				u = datax
				v = datay
			}

			prefix = ifelse(gg,i,NA_character_)
			qv = mapex2(u,v,frlow,doms,prefix,palette="PuRd")
			if (is.null(lstat2[[1]])) {
				lstat2[[1]] = array(dim=c(dim(qv),length(dates),length(fics)))
			}
			lstat2[[1]][,,,id,i] = qv
			rm(datax,datay)
		}

		if (! is.null(datal) && ! is.null(datam)) {
			cat(". compound field modulus\n")
			gsp = sqrt(datal^2+datam^2)
			print(summary(as.vector(gsp)))

			prefix = ifelse(gg,i,NA_character_)
			j2 = match("gradsp",descall$symbol)
			qv = mapexi(gsp,frlow,frmap,descall[j2,],doms,prefix)
			if (is.null(lstat2[[2]])) {
				lstat2[[2]] = array(dim=c(dim(qv),length(dates),length(fics)))
			}
			lstat2[[2]][,,,id,i] = qv
			rm(datal,datam)
		}

		if (! is.null(datagpl)) {
			cat(". computed field gradl\n")
			print(summary(as.vector(datagpl)))

			prefix = ifelse(gg,i,NA_character_)
			j2 = match("gradl",descall$symbol)
			qv = mapexi(datagpl,frlow,frmap,descall[j2,],doms,prefix)
			if (is.null(lstat2[[3]])) {
				lstat2[[3]] = array(dim=c(dim(qv),length(dates),length(fics)))
			}
			lstat2[[3]][,,,id,i] = qv
			rm(datagpl)
		}

		if (! is.null(dataox) && ! is.null(dataoy)) {
			cat(". compound field wind speed\n")
			ffo = sqrt(dataox^2+dataoy^2)
			print(summary(as.vector(ffo)))

			if (framo$gauss && latPole(framo) < 90) {
				cat("rotating wind to geo North\n")
				nord = compass(framo)
				if (! is.null(frlow$ind)) nord = nord[frlow$ind,]
				uo = ucs(nord,dataox,dataoy)
				vo = vcs(nord,dataox,dataoy)

				pole = list(xlim=lonPole(frame)+c(-6,6),ylim=latPole(frame)+c(-5,5),
					proj="stereo")
				x = zoom(dataox,frlow,pole)
				y = zoom(dataoy,frlow,pole)
				pngalt(sprintf("%s/mapdiff%dpole_ff.png",pngd,i))
				mapdom2(pole,x,x$data[,nl],y$data[,nl],main="u/v wind components")
				pngoff()
			}

			ffd = ff-ffo
			print(summary(as.vector(ffd)))
			drms = sqrt(mean(ffd^2))
			cat("Global bias/RMSE:",mean(ffd),drms,"\n")

			cat("Global bias/RMSE by level:\n")
			dbias = apply(ffd,2,mean)
			drms = apply(ffd,2,function(x) sqrt(mean(x^2)))
			print(summary(dbias))
			print(summary(drms))

			descff = descall[which(descall$symbol=="ff")[1],]
			prefix = ifelse(gg,sprintf("diff%d",i),NA_character_)
			qv = mapexi(ffd,frlow,frmap,descff,doms,prefix)
			if (is.null(lstat2d[[1]])) {
				lstat2d[[1]] = array(dim=c(dim(qv),length(dates),length(fics)))
				lerr2[[1]] = array(dim=c(dim(ffd),length(dates),length(fics)))
			}
			lstat2d[[1]][,,,id,i] = qv
			lerr2[[1]][,,id,i] = ffd
			rm(dataox,dataoy)
		}

		gc()
	}
}

if (! is.null(cargs$png)) {
	write(tstep,file=sprintf("%s/steps.txt",cargs$png))
	tt = grep(" TRUE",tstep,value=TRUE)
	write(sprintf("<option>%s</option>",tt),file=sprintf("%s/steps.html",cargs$png))
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
htime = htime/unit[iu]
tfreq = 1
if (unit[iu] == 3600 && diff(range(htime)) > 12) tfreq = 3

ht6 = pretty(htime/tfreq,7)*tfreq
if (all((ht6-ht6[1])%%5 == 0)) {
	tfreq = 6
	ht6 = pretty(htime/tfreq,7)*tfreq
}
etai = frlow$eta[frlow$ilev]
yeta = rev(range(frlow$eta))

saveStat("diag.RData",frlow,params,dates,htime,lstatd,lstat2d)

lstatr = NULL
if ("cmp" %in% names(cargs)) {
	ficsave = sprintf("%s/diag.RData",cargs$cmp)
	cat("Load cmp file",ficsave,"\n")
	lstatr = loadStat(ficsave,frlow,params,dates,htime)
}

ndom = length(doms)
nt = dim(lstat[[1]])[5]
if (nt < length(htime)) {
	length(htime) = nt
	length(ht6) = nt
}

i2 = sapply(lstat2,is.array)
if (any(i2)) {
	lstat = c(lstat,lstat2[i2])
	desc = rbind(desc,desc2[i2,])
}

if (length(fics) > 1) {
if (FALSE) {
	cat("Summary statistics of forecasts (whole domain, last (ie lowest) level)\n")
	tt = sprintf("Forecast of %s",desc$longname)

	nf = dim(desc)[1]
	nr = min(3,nf)

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/stats%d.png",pngd,i))
		op = par(c(Gpart,list(mfrow=c(nr,1))))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			ds = apply(dstats[,,,j,drop=FALSE],c(1,3),mean,na=TRUE)
			plotb(htime,ds,main=tt[j],xlab=tlab,ylab=desc$symbol[j],xaxt="n")
			axis(1,ht6)
		}

		pngoff()
	}
}

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
		stopifnot(all(dim(lstat[[j]]) == c(101,nl,ndom,length(dates),length(fics))))

		# dims are [val5,lev,dom,time]
		qv = apply(lstat[[j]],c(2,3,5),quantile,prob=c(0,2,5,8,10)/10,na.rm=TRUE)

		tt = sprintf("Statistics of %s",nom)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = c(1,1+nl%/%2,nl)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/stat%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				tt[2] = sprintf("domain %s",names(doms)[id])

				for (il in indl[1:nr]) {
					if (nl > 1) tt[2] = sprintf("%s, level %d",tt[2],frlow$ilev[il])
					plotb(htime,qv[,il,id,],main=tt,xlab=tlab,ylab=ss,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff()
		}

		if (nl == 1) next

		indt = which(apply(qv,4,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[c(1,length(indt)%/%2,length(indt))]
		nr = min(2,ndom)
		nc = min(3,length(indt))

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/statv%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfrow=c(nr,nc))))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				for (it in indt[1:nc]) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],tunit)
					plotbv(qv[,,id,it],etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
				}
			}

			pngoff()
		}
	}
}

if (any(! is.na(dstatd))) {
if (FALSE) {
	cat("Summary statistics of forecast error (whole domain, last (ie lowest) level)\n")
	tt = sprintf("Error of %s",desc$longname)

	nf = dim(desc)[1]
	nr = min(3,nf)

	for (i in seq((nf-1)%/%nr+1)-1) {
		pngalt(sprintf("%s/dstats%d.png",pngd,i))
		op = par(c(Gpart,list(mfrow=c(nr,2))))

		for (j in 1:min(nf-nr*i,nr)+nr*i) {
			ds = apply(dstatd[,,,j,drop=FALSE],c(1,3),mean,na.rm=TRUE)
			plotb(htime,ds,main=tt[j],xlab=tlab,ylab=ss,xaxt="n")
			abline(h=0,col="grey")
			axis(1,ht6)
		}

		pngoff()
	}
}

	cat("Scores of forecasts over",ndom,"domains\n")
	i2 = sapply(lstat2d,is.array)
	if (any(i2)) {
		lstatd = c(lstatd,lstat2d[i2])
		lerr = c(lerr,lerr2[i2])
	}

	for (j in seq(along=lstatd)) {
		ss = desc$symbol[j]
		nom = desc$longname[j]
		cat(". param",ss,"\n")

		if (is.null(lstatd[[j]])) {
			cat("--> no statistics of errors to plot\n")
			next
		}

		nl = dim(lstatd[[j]])[2]
		stopifnot(all(dim(lstatd[[j]]) == c(101,nl,ndom,length(dates),length(fics))))

		# dims are [val5,lev,dom,time]
		qv = apply(lstatd[[j]],c(2,3,5),quantile,prob=c(0,2,5,8,10)/10,na.rm=TRUE)
		rmsv = apply(lstatd[[j]],c(2,3,5),rms,na.rm=TRUE)
		rmsv = aperm(rmsv,c(2:4,1))
		cols = rep(c("black","royalblue1"),each=3)
		if (! is.null(lstatr)) {
			rmsvr = apply(lstatr[[j]],c(2,3,5),rms,na.rm=TRUE)
			rmsvr = aperm(rmsvr,c(2:4,1))
			rmsv = array(c(rmsv,rmsvr),c(dim(rmsv)[1:3],2*dim(rmsv)[4]))
		}

		tt = sprintf("Error of %s",nom)
		nr = min(3,nl)
		nc = min(2,ndom)
		indl = c(1,1+nl%/%2,nl)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/err%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				for (il in indl[1:nr]) {
					tt[2] = sprintf("domain %s, level %d",names(doms)[id],frlow$ilev[il])
					qvt = matrix(qv[,il,id,],nrow=5)
					plotb(htime,qvt,main=tt,xlab=tlab,ylab=ss,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff()
		}

		tt = sprintf("Score of %s",nom)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			pngalt(sprintf("%s/score%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				tt[2] = sprintf("domain %s",names(doms)[id])

				for (il in indl[1:nr]) {
					if (nl > 1) tt[2] = sprintf("%s, level %d",tt[2],frlow$ilev[il])
					rmst = matrix(rmsv[il,id,,],ncol=3)
					matplott(htime,rmst,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
					axis(1,ht6)
				}
			}

			pngoff()
		}

		indt = which(apply(qv,4,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[c(1,length(indt)%/%2,length(indt))]

		if (length(dates) > 1) {
			dd = as.Date(as.character(dates),format="%Y%m%d")
			indd = order(dd)
			nr = min(3,length(indt))
			rmsd = apply(lstatd[[j]],3:5,rms,na.rm=TRUE)
			rmsd = aperm(rmsd,c(2:4,1))
			if (! is.null(lstatr)) {
				rmsr = apply(lstatr[[j]],3:5,rms,na.rm=TRUE)
				rmsr = aperm(rmsr,c(2:4,1))
				rmsd = array(c(rmsd,rmsr),c(dim(rmsd)[1:3],2*dim(rmsd)[4]))
			}

			for (i in seq((ndom-1)%/%nc+1)-1) {
				pngalt(sprintf("%s/scoret%d_%s.png",pngd,i,ss))
				op = par(c(Gpart,list(mfcol=c(nr,nc))))

				for (id in 1:min(ndom-nc*i,nc)+nc*i) {
					for (it in indt[1:nr]) {
						tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],
							tunit)
						matplott(dd[indd],rmsd[id,indd,it,],xlab="Date",ylab=ss,main=tt,
							col=cols,xaxt="n")
						axis.Date(1,dd,format="%Y%m%d")
					}
				}

				pngoff()
			}

			if (nl == 1) next

			ttd = sprintf("RMSE of %s",nom)
			rmsvd = apply(lstatd[[j]],2:5,rms,na.rm=TRUE)

			for (i in seq((ndom-1)%/%nc+1)-1) {
				pngalt(sprintf("%s/rmsevt%d_%s.png",pngd,i,ss))
				op = par(c(Gparm,list(mfcol=c(nr,nc))))

				for (id in 1:min(ndom-nc*i,nc)+nc*i) {
					for (it in indt[1:nr]) {
						ttd[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],
							tunit)
						plotv(dd[indd],etai,t(rmsvd[1,,id,indd,it]),ylim=yeta,xlab="Date",
							ylab="eta",main=ttd,xaxt="n",las=0)
						axis.Date(1,dd,format="%Y%m%d")
					}
				}

				pngoff()
			}
		}

		if (nl == 1) next

		nr = min(2,ndom)
		nc = min(3,length(indt))

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/scorev%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfrow=c(nr,nc))))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				for (it in indt[1:nc]) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],
						tunit)
					matplotv(rmsv[,id,it,],etai,xlab=ss,ylab="eta",main=tt,ylim=yeta,col=cols)
				}
			}

			pngoff()
		}

		# back to err, after scores
		tt = sprintf("Error of %s",nom)

		for (i in seq((ndom-1)%/%nr+1)-1) {
			pngalt(sprintf("%s/errv%d_%s.png",pngd,i,ss))
			op = par(c(Gpart,list(mfrow=c(nr,nc))))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				for (it in indt[1:nc]) {
					tt[2] = sprintf("domain %s, lead-time %g%s",names(doms)[id],htime[it],
						tunit)
					plotbv(qv[,,id,it],etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
				}
			}

			pngoff()
		}
	}

	cat("Compute statistics (mean bias, RMSE and max error) and graphics\n")
	descb = desc
	descb$palette = "Blue-Red+"

	for (j in seq(along=lstatd)) {
		cat(". param",desc$symbol[j],"\n")

		for (i in seq(along=dff$file)) {
			cat(". mean bias\n")
			# efficient mean
			bias = lerr[[j]][,,1,i]
			for (id in seq(along=dates)[-1]) bias = bias+lerr[[j]][,,id,i]
			bias = bias/length(dates)
			dim(bias) = dim(lerr[[j]])[1:2]
			prefix = ifelse(gg,sprintf("bias%d",i),NA_character_)
			invisible(mapexi(bias,frlow,frmap,descb[j,],doms,prefix))

			cat(". RMSE\n")
			rmse = lerr[[j]][,,1,i]^2
			for (id in seq(along=dates)[-1]) rmse = rmse+lerr[[j]][,,id,i]^2
			rmse = sqrt(rmse/length(dates))
			dim(rmse) = dim(lerr[[j]])[1:2]
			prefix = ifelse(gg,sprintf("rmse%d",i),NA_character_)
			invisible(mapexi(rmse,frlow,frmap,desc[j,],doms,prefix))

			cat(". maximum error\n")
			errx = abs(lerr[[j]][,,1,i])
			dim(errx) = dim(lerr[[j]])[1:2]
			for (id in seq(along=dates)[-1]) errx = pmax(errx,abs(lerr[[j]][,,id,i]))
			prefix = ifelse(gg,sprintf("errx%d",i),NA_character_)
			invisible(mapexi(errx,frlow,frmap,desc[j,],doms,prefix))

			#cat(". index of maximum error\n")
			#jmax = errx
			#for (ip in seq(prod(dim(errx)))) {
				#err = lerr[[j]][,,1,i]
				#jmax[ip] = sapply(lerr,function(j) which.min(abs(errx[ip]-err[ip])))
			#}
			#prefix = ifelse(gg,sprintf("errx%d",i),NA_character_)
			#invisible(mapexi(errx,frlow,frmap,desc[j,],doms,prefix))
		}
	}
}
