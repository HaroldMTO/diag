Gdiag = "~/util/diag"

library(maps)
library(mapproj)
library(parallel)
library(mffield)

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

rms = function(x,...)
{
	sqrt(mean(x^2,...))
}

rmx = function(x,...)
{
	c(rms=rms(x,...),ave=mean(x,...),q90=quantile(abs(x),prob=.9,...))
}

matplott = function(x,y,col=1,pch="+",lty=1:3,x.leg="topleft",legend=c("RMS","bias","Q9"),
	...)
{
	matplot(x,y,type="o",lty=lty,col=col,pch=pch,...)
	abline(h=0,col="darkgrey")
	if (! is.null(legend)) legend(x.leg,legend,col=col,lty=lty,pch=pch,bg="transparent")
}

matplotv = function(x,y,col=1,pch="+",lty=1:3,x.leg="topleft",legend=c("RMS","bias","Q9"),
	...)
{
	matplot(x,y,type="o",lty=lty,col=col,pch=pch,...)
	abline(v=0,col="darkgrey")
	if (! is.null(legend)) legend(x.leg,legend,col=col,lty=lty,pch=pch,bg="transparent")
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
		cat("--> some params missing in ref:",paramsp[is.na(indp)],"\n")
		return(NULL)
	}

	indz = match(domsp,doms)
	if (any(is.na(indz))) {
		cat("--> some domains missing in ref:",domsp[is.na(indz)],"\n")
		return(NULL)
	}

	indd = match(datesp,dates)
	if (any(is.na(indd))) {
		cat("--> some dates missing in ref:",datesp[is.na(indd)],"\n")
		return(NULL)
	}

	indt = match(htimep,htime)
	if (any(is.na(indt))) {
		cat("--> some time-steps missing in ref:",htimep[is.na(indt)],"\n")
		return(NULL)
	}

	lapply(lstatd[indp],function(x) x[,,indz,indd,indt,drop=FALSE])
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
		cat("--> some params missing in ref:",paramsp[is.na(indp)],"\n")
		return(NULL)
	}

	indd = match(datesp,dates)
	if (any(is.na(indd))) {
		cat("--> some dates missing in ref:",datesp[is.na(indd)],"\n")
		return(NULL)
	}

	indt = match(htimep,htime)
	if (any(is.na(indt))) {
		cat("--> some time-steps missing in ref:",htimep[is.na(indt)],"\n")
		return(NULL)
	}

	lapply(lzmeand[indp],function(x) x[,,indd,indt,drop=FALSE])
}

source(sprintf("%s/plot.R",Gdiag))

Gpart = list(mar=c(2.5,2.5,2.5,1),mgp=c(1.5,.5,0),tcl=-.3,cex=.8)
Gparmt = list(mar=c(2.5,2.5,2.5,4.8),mgp=c(1.5,.5,0),tcl=-.3)

doms = readDom(sprintf("%s/config/domain.txt",Gdiag))
if (file.exists("config/domain.txt")) doms = readDom("config/domain.txt",doms)

descall = read.table(sprintf("%s/config/params.txt",Gdiag),header=TRUE)

args = strsplit(commandArgs(trailingOnly=TRUE),split="=")
cargs = lapply(args,function(x) unlist(strsplit(x[-1],split=":")))
names(cargs) = sapply(args,function(x) x[1])

pngd = "."
if ("png" %in% names(cargs)) pngd = cargs$png
if (! "cmp" %in% names(cargs)) cargs$cmp = character()

if ("png" %in% names(cargs) || ! capabilities("X11")) {
	cat("--> sending plots to PNG files in",pngd,"\n")
} else {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
}

wide = FALSE
if ("wide" %in% names(cargs)) wide = as.logical(cargs$wide)

etahigh = 0
if ("etahigh" %in% names(cargs)) etahigh = as.numeric(cargs$etahigh)

if (interactive()) browser()

if (any(duplicated(basename(cargs$cmp)))) stop("duplicated directory names")

lstatr = lzmeanr = list()

for (cmp in c(".",cargs$cmp)) {
	stopifnot(file.exists(cmp) && file.info(cmp)$isdir)
	ficsave = sprintf("%s/diag.RData",cmp)
	cat("Load cmp file",ficsave,"\n")

	n = length(lstatr)+1

	if (n == 1) {
		noms = load(ficsave)
		stopifnot(all(c("frlow","params","dates","htime","lstatd") %in% noms))
		if (is.list(doms)) {
			cat("--> no domains in diag file, use the domains read\n")
			doms = names(doms)
		}

		eta = frlow$eta
		ilev = frlow$ilev
		indt = order(htime)
		if (length(indt) > 1 && any(indt != seq(along=htime))) {
			for (i in seq(along=lstatd)) lstatd[[i]] = lstatd[[i]][,,,,indt,drop=FALSE]
			htime = htime[indt]
		}
	} else {
		lstatd = loadStat(ficsave,frlow,params,doms,dates,htime)
		if (is.null(lstatd)) stop("no stat in common\n")
	}

	lstatr = c(lstatr,list(lstatd))
	names(lstatr)[n] = basename(cmp)

	lzmeand = try(loadZmean(ficsave,frlow,params,dates,htime),silent=TRUE)
	if (is(lzmeand,"try-error")) {
		ficsave = sprintf("%s/diagz.RData",cmp)
		cat("Load cmp file for zonal mean",ficsave,"\n")
		if (! file.exists(ficsave)) {
			cat("--> no file",ficsave,"\n")
			next
		}

		lzmeand = loadZmean(ficsave,frlow,params,dates,htime)
	}

	lzmeanr = c(lzmeanr,list(lzmeand))
	names(lzmeanr)[n] = basename(cmp)
}

#if (length(lstatr) < 2) stop("no stat to compare")
if (any(sapply(lzmeanr,is.null))) lzmeanr = NULL
lstatd = lstatr[[1]]
lzmeand = lzmeanr[[1]]

ind = match(params,descall$faname)
if (any(is.na(ind))) stop("unknown parameters, see config/params.txt")

desc = descall[ind,]

if (length(htime) > 2 && diff(htime)[1] == 0) {
	htime[1] = -diff(htime[2:3])/3
	cat("--> changing 1st time-step for graphics:",htime[1:3],"...\n")
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

etai = eta[ilev]
yeta = rev(range(eta))

ndom = length(doms)
nday = length(dates)
nt = length(htime)

cat("Scores of forecasts over",ndom,"domains\n")
for (j in seq(along=lstatd)) {
	ss = desc$symbol[j]
	nom = desc$longname[j]
	nl = dim(lstatd[[j]])[2]
	cat(". param",ss,"- levels:",nl,"\n")

	stopifnot(all(dim(lstatd[[j]]) == c(101,nl,ndom,nday,nt)))

	# cmp[stat, level, domain, date, step] -> qvalues, level, domain, step
	qe = apply(lstatd[[j]],c(2,3,5),quantile,prob=seq(0,50)/50,na.rm=TRUE)

	tt = sprintf("Error of %s",nom)
	nr = min(3,nl)
	nc = min(2,ndom)
	indl = round(seq(0,nl,length.out=nr+1)[-1])

	# qe err: ll box 3lx2
	for (i in seq((ndom-1)%/%nc+1)-1) {
		png(sprintf("%s/err%d_%s.png",pngd,i,ss))
		par(c(Gpart,list(mfcol=c(nr,nc))))

		for (id in 1:min(ndom-nc*i,nc)+nc*i) {
			il0 = 1
			tt[2] = sprintf("domain %s",doms[id])
			for (il in indl) {
				if (nl > 1) {
					tt[2] = sprintf("domain %s, levels [%d,%d]",doms[id],
						frlow$ilev[il0],frlow$ilev[il])
				}
				qet = matrix(qe[,il0:il,id,],ncol=length(ht))
				plotb(ht,qet,main=tt,xlab=tlab,ylab=ss,xaxt="n")
				axis(1,ht6)
				il0 = min(il+1,max(indl))
			}
		}

		dev.off()
	}

	if (nl > 1) {
		indt = which(apply(lstatd[[j]],5,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[seq(1,length(indt),len=3)]

		nc = length(indt)
		if (etahigh > 0) {
			indl = which(etai < etahigh)
			if (length(indl) < 3) indl = which(etai < .5)
			if (length(indl) < 3) stop(sprintf("not enough levels above .5"))
			nr = 2
			nj = 1
			tth = tt
			tth[1] = sprintf("..., high levels",nom)
		} else if (wide) {
			nr = nj = 1
		} else {
			nr = nj = min(2,ndom)
		}

		# qe errv: box v 2x3t
		for (i in seq((ndom-1)%/%nj+1)-1) {
			png(sprintf("%s/errv%d_%s.png",pngd,i,ss))
			par(c(Gpart,list(mfrow=c(nr,nc))))

			for (id in 1:min(ndom-nj*i,nj)+nj*i) {
				it0 = 1
				if (etahigh > 0) {
					for (it in indt) {
						tth[2] = sprintf("dom. %s, lead-time [%g,%g]%s",doms[id],
							ht[it0],ht[it],tunit)
						qet = qe[,indl,id,it0:it]
						if (length(dim(qet)) == 3) {
							qet = matrix(aperm(qet,c(3,1:2)),ncol=length(indl))
						}
						plotbv(qet,etai[indl],main=tth,xlab=ss,ylab="eta",
							ylim=c(etahigh,yeta[2]))
					}
				}

				for (it in indt) {
					tt[2] = sprintf("dom. %s, lead-time [%g,%g]%s",doms[id],
						ht[it0],ht[it],tunit)
					qet = qe[,,id,it0:it]
					if (length(dim(qet)) == 3) qet = matrix(aperm(qet,c(3,1:2)),ncol=nl)
					plotbv(qet,etai,main=tt,xlab=ss,ylab="eta",ylim=yeta)
					it0 = min(it+1,max(indt))
				}
			}

			dev.off()
		}
	}

	cols = c("black","royalblue1","orangered3","forestgreen","mediumpurple3")
	cols = rep(cols[seq(along=lstatr)],each=3)
	cat("--> colors:",unique(cols),"\n")

	cat("Scores by lead-time, levels in rows\n")
	# cmp[stat, level, domain, dates, step] -> rmx, level, domain, step, cmp
	rmxv = sapply(lstatr,function(statd) apply(statd[[j]],c(2,3,5),rmx,na.rm=TRUE),
		simplify="array")
	# rmx, level, domain, step, cmp -> level, domain, step, rmx, cmp
	rmxv = aperm(rmxv,c(2:4,1,5))

	tt = sprintf("Score of %s",nom)
	nr = min(3,nl)
	nc = min(2,ndom)
	if (wide) nc = 1
	indl = as.integer(seq(1,nl,length.out=nr))

	# rmxv score: ll plot 3lx2
	for (i in seq((ndom-1)%/%nc+1)-1) {
		png(sprintf("%s/score%d_%s.png",pngd,i,ss))
		par(c(Gpart,list(mfcol=c(nr,nc))))

		for (id in 1:min(ndom-nc*i,nc)+nc*i) {
			tt[2] = sprintf("domain %s",doms[id])
			for (il in indl) {
				if (nl > 1) tt[2] = sprintf("domain %s, level %d",doms[id],ilev[il])
				rmxvt = matrix(rmxv[il,id,,,],nrow=nt)
				matplott(ht,rmxvt,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
				axis(1,ht6)
			}
		}

		dev.off()
	}

	if (nl > 1) {
		tt = sprintf("Score vert. mean of %s",nom)
		rmxm = sapply(lstatr,function(statd) apply(statd[[j]],c(3,5),rmx,na.rm=TRUE),
			simplify="array")
		# rmx, domain, step, cmp -> domain, step, rmx, cmp
		rmxm = aperm(rmxm,c(2:3,1,4))

		# rmxm score: plot 2dx2d
		off = (ndom-1)%/%nc+1
		nr = min(2,ndom)
		nc = min(2,1+(ndom-1)%/%2)
		nj = nr*nc
		for (i in seq((ndom-1)%/%nj+1)-1) {
			png(sprintf("%s/score%d_%s.png",pngd,off+i,ss))
			par(c(Gpart,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nj*i,nj)+nj*i) {
				tt[2] = sprintf("domain %s",doms[id])
				rmxt = matrix(rmxm[id,,,],nrow=nt)
				matplott(ht,rmxt,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
				axis(1,ht6)
			}

			dev.off()
		}

		cat("Profile of scores by lead-time\n")
		tt = sprintf("Score of %s",nom)
		indt = which(apply(lstatd[[j]],5,function(x) any(! is.na(x))))
		if (length(indt) > 3) indt = indt[seq(1,length(indt),len=3)]
		nc = length(indt)
		if (etahigh > 0) {
			indl = which(etai < etahigh)
			if (length(indl) < 3) indl = which(etai < .5)
			if (length(indl) < 3) stop(sprintf("not enough levels above .5"))
			nr = 2
			tth = tt
			tth[1] = sprintf("..., high levels",nom)
		} else if (wide) {
			nr = 1
		} else {
			nr = min(2,ndom)
		}

		# rmxv scorev: v plot 2x3t
		for (i in seq((ndom-1)%/%nr+1)-1) {
			png(sprintf("%s/scorev%d_%s.png",pngd,i,ss))
			par(c(Gpart,list(mfrow=c(nr,nc))))

			for (id in 1:min(ndom-nr*i,nr)+nr*i) {
				if (etahigh > 0) {
					for (it in indt) {
						tth[2] = sprintf("dom. %s, lead-time %g%s",doms[id],ht[it],tunit)
						rmxvl = matrix(rmxv[indl,id,it,,],nrow=length(indl))
						matplotv(rmxvl,etai[indl],x.leg="topright",xlab=ss,ylab="eta",main=tth,
							ylim=c(etahigh,yeta[2]))
					}
				}

				for (it in indt) {
					tt[2] = sprintf("dom. %s, lead-time %g%s",doms[id],ht[it],tunit)
					rmxvl = matrix(rmxv[,id,it,,],nrow=nl)
					matplotv(rmxvl,etai,x.leg="topright",xlab=ss,ylab="eta",main=tt,ylim=yeta,
						col=cols)
				}
			}

			dev.off()
		}
	}

	if (length(lzmeand) >= j) {
		cat("Scores of zonal mean\n")
		# cmp[lat, level, dates, step] -> rmx, level, step, cmp
		# dims are [lat,lev,date,time]
		# rmxv score: ll plot 6l
		zrmxv = sapply(lzmeanr,function(zmean) apply(zmean[[j]],c(2,4),rmx,na.rm=TRUE),
			simplify="array")
		zrmxv = aperm(zrmxv,c(2:3,1,4))

		tt = sprintf("Zonal score of %s",nom)
		nr = min(3,nl%/%2)
		nc = min(2,nl%/%2)
		if (nl == 1) nr = nc = 1
		indl = as.integer(seq(1,nl,length.out=nr*nc))

		png(sprintf("%s/scorez1_%s.png",pngd,ss))
		par(c(Gpart,list(mfcol=c(nr,nc))))

		for (il in indl) {
			if (nl > 1) tt[2] = sprintf("level %d",frlow$ilev[il])
			zt = matrix(zrmxv[il,,,],nrow=length(ht))
			matplott(ht,zt,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
			axis(1,ht6)
		}

		dev.off()

		if (nl > 1) {
			cat("Scores of zonal mean - maps\n")
			tt = sprintf("Score zon.+vert. mean, %s",nom)
			zrmx = sapply(lzmeanr,function(zmean) t(apply(zmean[[j]],4,rmx,na.rm=TRUE)),
				simplify="array")
			zrmx = matrix(zrmx,nrow=length(ht))

			png(sprintf("%s/scorez2_%s.png",pngd,ss))
			par(Gpart)

			matplott(ht,zrmx,xlab=tlab,ylab=ss,main=tt,col=cols,xaxt="n")
			axis(1,ht6)

			dev.off()

			tt = sprintf("Zonal score of %s",nom)
			indt = which(apply(lzmeand[[j]],4,function(x) any(! is.na(x))))
			nr = min(2,length(indt)%/%2)
			nc = min(3,length(indt)%/%2)
			if (length(indt) == 1) nr = nc = 1
			if (length(indt) > nr*nc) indt = indt[seq(1,length(indt),len=nr*nc)]

			png(sprintf("%s/scorevz1_%s.png",pngd,ss))
			par(c(Gpart,list(mfrow=c(nr,nc))))

			for (it in indt) {
				tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
				zv = matrix(zrmxv[,it,,],nrow=nl)
				matplotv(zv,etai,x.leg="topright",ylim=yeta,xlab=ss,ylab="eta",main=tt,
					col=cols)
			}

			dev.off()
		}

		tt = sprintf("Zonal bias of %s",nom)
		indt = which(apply(lzmeand[[j]],4,function(x) any(! is.na(x))))
		nr = min(3,length(indt))
		nc = 1
		if (length(indt) == 1) nr = nc = 1
		if (length(indt) > nr*nc) indt = indt[seq(1,length(indt),len=nr*nc)]
		lat = sort(frlow$g4@theta)
		latx = pretty(lat,8)

		png(sprintf("%s/scorevz2_%s.png",pngd,ss))
		par(c(Gparmt,list(mfrow=c(nr,nc))))

		if (nl == 1) {
			#zrmx = sapply(lzmeand,function(zmean) t(apply(zmean[[j]],c(1,4),rmx,na.rm=TRUE)),
			#	simplify="array")
			#zrmx = aperm(zrmx,c(2:3,1,4))
			zrmx = apply(lzmeand[[j]],c(1,4),mean,na.rm=TRUE)

			for (it in indt) {
				tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
				zz = matrix(zrmx[,it],nrow=length(lat))
				matplott(lat,zz,main=tt,xlab="Latitude",ylab="Bias",col=cols,xaxt="n")
				axis(1,latx)
			}
		} else {
			zbiasv = apply(lzmeand[[j]],c(1:2,4),mean,na.rm=TRUE)

			for (it in indt) {
				tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
				plotv(lat,etai,zbiasv[,,it],main=tt,ylim=yeta,xlab="Latitude",ylab="eta",
					palette="BlueRed+",xaxt="n")
				axis(1,latx)
			}
		}

		dev.off()

		if (nl > 1) {
			tt = sprintf("Zonal RMSE of %s",nom)
			zrmsev = apply(lzmeand[[j]],c(1:2,4),rms,na.rm=TRUE)

			png(sprintf("%s/scorevz3_%s.png",pngd,ss))
			par(c(Gparmt,list(mfrow=c(nr,nc))))

			for (it in indt) {
				tt[2] = sprintf("lead-time %g%s",ht[it],tunit)
				plotv(lat,etai,zrmsev[,,it],main=tt,ylim=yeta,xlab="Latitude",ylab="eta",
					xaxt="n")
				axis(1,latx)
			}

			dev.off()
		}
	}

	if (length(dates) < 2) next

	cat("Scores for",length(dates),"days\n")
	dd = as.Date(as.character(dates),format="%Y%m%d")
	indt = which(apply(lstatd[[j]],5,function(x) any(! is.na(x))))
	if (length(indt) > 3) indt = indt[seq(1,length(indt),len=3)]
	indd = order(dd)
	nr = length(indt)
	nc = min(2,ndom)

	# cmp[stat, level, domain, dates, step] -> rmx, domain, dates, step, cmp
	rmxt = sapply(lstatr,function(statd) apply(statd[[j]],3:5,rmx,na.rm=TRUE),
		simplify="array")
	# rmx, domain, dates, step, cmp -> domain, dates, step, rmx, cmp
	rmxt = aperm(rmxt,c(2:4,1,5))

	# rmxt scoret: day t plot 3tx2
	for (i in seq((ndom-1)%/%nc+1)-1) {
		png(sprintf("%s/scoret%d_%s.png",pngd,i,ss))
		par(c(Gpart,list(mfcol=c(nr,nc))))

		for (id in 1:min(ndom-nc*i,nc)+nc*i) {
			for (it in indt) {
				tt[2] = sprintf("domain %s, lead-time %g%s",doms[id],ht[it],tunit)
				rmxtd = matrix(rmxt[id,indd,it,,],nrow=length(indd))
				matplott(dd[indd],rmxtd,main=tt,xlab="Date",ylab=ss,col=cols,xaxt="n",las=3)
				axis.Date(1,dd,format="%Y%m%d")
			}
		}

		dev.off()
	}

	if (nl > 1) {
		cat("Score RMSE for days - profile\n")
		# rmxvd rmsev: day v image 3tx2
		tt = sprintf("RMSE of %s",nom)
		rmxvd = apply(lstatd[[j]],2:5,rms,na.rm=TRUE)

		for (i in seq((ndom-1)%/%nc+1)-1) {
			png(sprintf("%s/rmsevt%d_%s.png",pngd,i,ss))
			par(c(Gparmt,list(mfcol=c(nr,nc))))

			for (id in 1:min(ndom-nc*i,nc)+nc*i) {
				br = prettyBreaks(rmxvd[,id,,indt[1:nr]],crop=TRUE)$breaks

				for (it in indt) {
					tt[2] = sprintf("domain %s, lead-time %g%s",doms[id],ht[it],tunit)
					plotv(dd[indd],etai,t(rmxvd[,id,indd,it]),br,main=tt,ylim=yeta,
						xlab="Date",ylab="eta",xaxt="n",las=3)
					axis.Date(1,dd,format="%Y%m%d")
				}
			}

			dev.off()
		}
	}
}
