readDef = function(nd)
{
	linfo = vector("list",9)

	info = c("statkind","surfaces","types","params","items","instrument","flagfilter")
	names(linfo) = info
	ind = c(2,4:9)
	for (i in seq(along=info)) {
		val = sub(sprintf("^ *%s *= *([^']+) *",info[i]),"\\1",nd[ind[i]])
		linfo[[i]] = as.integer(strsplit(sub(" +\\(.+\\)","",val)," +")[[1]])
	}

	sinfo = c("comment","areaNSEW")
	ind = c(1,3)
	ni = 7
	names(linfo)[-(1:ni)] = sinfo
	for (i in seq(along=sinfo)) {
		val = sub(sprintf("^ *%s *= *'?([^']+)'? *",sinfo[i]),"\\1",nd[ind[i]])
		linfo[[ni+i]] = sub("(.+) +\\(.+\\)","\\1",val)
	}

	stopifnot(nzchar(linfo$comment))
	stopifnot(nzchar(linfo$areaNSEW))

	if (length(nd) == length(linfo)) {
		return(linfo)
	} else if (length(nd) == length(linfo)+1) {
		binfo = "nbbin"
		ind = 10
	} else if (length(nd) == length(linfo)+3) {
		binfo = c("sizebin","refval","nbbin")
		ind = 10:12
	} else if (length(nd) == length(linfo)+4) {
		binfo = c("nbbin","coorditem1","sizebin1?","refval1?")
		ind = 10:13
	} else {
		stop("unknown statdef")
	}

	ni = length(linfo)
	for (i in seq(along=binfo)) {
		val = sub(sprintf("^ *%s *= *'?([^']+)'? *",binfo[i]),"\\1",nd[ind[i]])
		linfo[[ni+i]] = as.numeric(val)
	}

	names(linfo)[-(1:ni)] = binfo
	linfo$nbbin = as.integer(linfo$nbbin)
	linfo
}

readItems = function(nd,ii,items,statkind)
{
	indi = grep("^ *BEGIN STATITEM",nd[ii])+ii[1]
	indie = grep("^ *END STATITEM",nd[ii])+ii[1]-2
	inul = which(indi > indie)
	if (length(inul) > 0) {
		#cat("--> remove",length(inul),"nul items\n")
		items = items[-inul]
	}

	if (length(items) > 6) {
		#cat("--> limiting items to 6 out of",length(items),"\n")
		length(items) = 6
	}

	litems = lapply(seq(along=items),function(j) readItem(nd[indi[j]:indie[j]],statkind))

	litems
}

readItem = function(nd,stat)
{
	item = sub("^ *#item *: *(.+) *","\\1",nd[1])
	item = sub(" +$","",item)

	ind = 1

	if (stat == 2) {
		stopifnot(regexpr("^ *# *pop +min +max",nd[2]) > 0)
		val = as.numeric(strsplit(nd[3]," +")[[1]])
		ind = 1:3
	}

	con = file()
	writeLines(nd[-ind],con)
	df = read.table(con,header=TRUE,comment.char="")
	close(con)
	names(df)[1] = sub("X\\.","",names(df)[1])
	list(item=item,data=df)
}

args = commandArgs(trailingOnly=TRUE)

fobs = args[1]
pngd = NULL
if (length(args) == 2) pngd = args[2]

if (! is.null(pngd) || ! capabilities("X11")) {
	if (is.null(pngd)) pngd = "obstat"
	cat("--> no X11 device, sending plots to PNG files in",pngd,"\n")
} else {
	png = dev.off = function(...) return(invisible(NULL))
	if (interactive()) options(device.ask.default=TRUE)
}

cat("Read file",fobs,"\n")
nd = readLines(fobs)

ind = grep("^ *BEGIN STATDEF",nd)+1
inde = grep("^ *END STATDEF",nd)-1
stopifnot(length(ind) == length(inde))

ii = which(ind > inde)
nstat = length(ind)
n0 = length(ii)
cat("--> parse",nstat-n0,"statdef, with",n0,"nul def\n")

# for next "BEGIN STATDEF"
ind = c(ind,length(nd))

defs = vector("list",length(inde))

if (interactive()) browser()

for (i in seq(along=inde)) {
	def = readDef(nd[ind[i]:inde[i]])
	s = paste(def$params,collapse=":")
	if (length(def$params) != 1)
		cat(". def:",def$comment,def$areaNSEW,"- params/kind:",s,def$statkind,"\n")

	# items
	ii = (inde[i]+1):(ind[i+1]-1)
	def$items = readItems(nd,ii,def$items,def$statkind)

	defs[[i]] = def

	if (length(def$items) == 0) next

	tt = paste(def$comment,def$areaNSEW)
	if (any(def$instrument != 999)) {
		tt[2] = sprintf("instrument(s): %s",paste(def$instrument,collapse=" "))
	}
	s = gsub(" ","-",def$comment)
	if (regexpr("-",s) < 0) s = sprintf("%s-X",s)
	ficpng = sprintf("%s/%s_%s.png",pngd,s,gsub("\\.","",def$areaNSEW))
	#cat(". file",ficpng,"\n")
	png(ficpng)

	nc = max(1,length(def$items)%/%3)
	nr = length(def$items)%/%nc
	par(mfrow=c(nr,nc),mgp=c(2,1,0))

	if (def$statkind == 2) {
		for (j in seq(along=def$items)) {
			item = def$items[[j]]
			barplot(item$data$population,col="grey98",space=0,main=tt,xlab=item$item,
				ylab="population")
			axis(1)
		}
	} else {
		for (j in seq(along=def$items)) {
			item = def$items[[j]]
			ylim = range(item$data[,1])
			if (names(item$data)[1] == "Pressure") ylim = rev(ylim)
			matplot(item$data[,-(1:2)],item$data[,1],type="o",main=tt,xlab=item$item,
				ylab=names(item$data)[1],ylim=ylim,lty=c(0,0,2,1),pch=c("-","+",NA,NA),
				col=1)
			#abline(h=0)
			abline(v=0,col="darkgrey",lty=2)
		}
	}

	dev.off()
}

comm = sapply(defs,"[[","comment")
area = sapply(defs,"[[","areaNSEW")
tt = paste(comm,area)

it = which(duplicated(tt))
if (length(it) > 0) {
	tt1 = unique(tt[it])
	idup = integer()

	for (i in seq(along=tt1)) {
		ind = which(tt == tt1[i])
		items = defs[[ind[1]]]$items
		ii = which(sapply(ind[-1],function(j) identical(items,defs[[j]]$items)))
		if (length(ii) > 0) {
			#cat("--> duplicates for index",ind[1],":",ind[ii+1],"\n")
			idup = c(idup,ind[ii+1])
		}
	}
}

if (length(idup) > 0) {
	cat("--> duplicates:",idup,"\n")
	defs = defs[-idup]
}

i = unlist(sapply(defs,function(x) x$instrument))
s = unlist(sapply(defs,function(x) x$statkind))
p = unlist(sapply(defs,function(x) x$params))
f = sapply(defs,function(x) x$flagfilter)
t = unlist(sapply(defs,function(x) x$types))
a = sapply(defs,function(x) x$areaNSEW)
n = sapply(defs,function(x) length(x$items))
cat("Summary:\n")
cat(". instruments:",sort(unique(i)),"\n")
cat(". params:",sort(unique(p)),"\n")
cat(". types:",sort(unique(t)),"\n")
cat(". flagfilter:\n")
print(table(f))
cat(". statkind:\n")
print(table(s))
cat(". area:\n")
print(table(a))
