#!/usr/bin/env python

import os

def usage():
	print("Syntax:\n\
	epy_dump.py FILE -f frame|field_name -o FILEOUT [-prof] [-h]\n\
\n\
Description:\n\
	Dump information (geometry) or field data from a file read by Epygram\n\
\n\
Arguments:\n\
	FILE: name of the file to be read\n\
	frame: target information on the date and geometry of the fields in the file\n\
	field_name: target access to the data of one field of the file\n\
	FILEOUT: name of the (binary) file to write to\n\
	-prof: activate simple profiling in this script\n\
	-h: print this help and exit\n\
\n\
Details:\n\
	Information is dumped into a FILEOUT, containing data in binary format.\n\
	For grid option, it consists in dimensions and vectors of horizontal and vertical \
grids:\n\
		- dimensions contain numbers (4) of lats, waves (0 if spectral), levels and \
grid-points.\n\
		- vectors (2) are numbers of longs per lat and truncation number per wave number.\n\
	For field_name option, it consists in the data of the field, pre-pended by its size if \
spectral.\n")

tag = ""
fbin = ""
prof = False

args = os.sys.argv[1:]
if len(args) == 0:
	usage()
	exit()

for i,arg in enumerate(args):
	if arg == "-h":
		usage()
		exit()

	if arg == "-f":
		if tag != "": exit("too many options, see usage")

		tag = args[i+1]
	elif arg == "-o":
		if fbin != "": exit("too many options, see usage")

		fbin = args[i+1]
	elif arg == "-prof":
		prof = True
	elif (tag == "" or arg != tag) and (fbin == "" or arg != fbin):
		fname = arg

import epygram
import numpy
import time
import concurrent.futures as cf

def sp2gp(f):
	f.sp2gp()
	return f

epygram.init_env()

print("Reading file "+fname)
a = epygram.formats.resource(fname,"r")

if a.format != "GRIB":
	geom = a.geometry
	grid = geom.grid
	dimg = geom.dimensions

	islam = geom.name in ("lambert","mercator","regular_lonlat")
	if islam:
		print("Geometry:",geom.name,"(Limited Area)")
		nlat = dimg["Y"]
		nlon = dimg["X"]
	else:
		print("Geometry:",geom.name,"(global)")
		if geom.name not in ("gauss","reduced_gauss","rotated_reduced_gauss"):
			exit("error unknown geom name")

		nlon = dimg["lon_number_by_lat"]
		if geom.name == "rotated_reduced_gauss":
			locen = grid["pole_lon"].get("degrees")
			mucen = grid["pole_lat"].get()[1]
			gem = (grid["dilatation_coef"],mucen,locen)
		else:
			gem = (1.,1.,0.)

	lonlat = geom.get_lonlat_grid()
	longs = lonlat[0]
	lats = lonlat[1]
	lats = lats.compressed()
	longs = longs.compressed()

	vv = geom.vcoordinate
	if len(vv.levels) == 1:
		exit("one level only")

	vab = vv.grid["gridlevels"]
	Ai = numpy.empty(len(vab))
	Bi = numpy.empty(len(vab))
	for i,vi in enumerate(vab):
		Ai[i] = vi[1]["Ai"]
		Bi[i] = vi[1]["Bi"]

	spgeom = a.spectral_geometry
	tc = spgeom.truncation
	if spgeom.space == "legendre" and tc["shape"] == "triangular":
		nwave = tc["max_zonal_wavenumber_by_lat"]
		# last 0 is for recognition of global grid vs LAM grid
		dims = len(nlon),len(nwave),len(vv.levels),geom.gridpoints_number,0
		print("Dimensions (gaussian lats, waves, levels, grid-points):",dims[:4])
	elif spgeom.space == "bi-fourier" and tc["shape"] == "elliptic":
		nwavex = tc["in_X"]
		nwavey = tc["in_Y"]
		dims = nlat,nlon,nwavex,nwavey,len(vv.levels)
		print("Dimensions (Lambert lats and longs, X and Y waves, levels):",dims)
	else:
		exit("unknown space of spectral geometry")

con = None
if fbin != "":
	con = open(fbin,"wb")

if tag == "frame":
	base = a.validity.getbasis()
	step = a.validity.term()
	print("Base date/time:",base.strftime("%Y%m%d %X"))
	print("Step:",step.total_seconds()/3600,"h")

	if con != None:
		numpy.array(dims).tofile(con)
		if not islam:
			numpy.array(nlon).tofile(con)
			numpy.array(nwave).tofile(con)

		bb = numpy.array(base.timestamp())
		ss = numpy.array(step.total_seconds())
		print("Base/step:",bb,ss)
		bb.tofile(con)
		ss.tofile(con)
		lats.tofile(con)
		longs.tofile(con)
		Ai.tofile(con)
		Bi.tofile(con)
		if not islam: numpy.array(gem).tofile(con)
elif tag != "":
	print(". read fields "+tag)
	if tag == "list":
		l = a.listfields()
		print(l)
		exit()

	if prof: t = os.times().elapsed
	fields = a.readfields(tag)
	if prof: print("reading time:",round(os.times().elapsed-t,3))
	fnames = fields.listfields("FA")
	fnoms = numpy.array(fnames)

	nf = len(fields)
	if nf == 0:
		exit("Error: no field found")

	print("--> found",len(fields),"field(s)")

	if nf > 1:
		import re
		fheads = fields.listfields()

		if re.search('^S\d+\w+',fheads[0].get("FA")) is None:
			exit("Error: field names do not match level form 'S\d+\w+'")

		for h in fheads:
			if re.search('^S\d+\w+',h.get('FA')) is None:
				exit("Error: mixing field names")

	if con != None:
		numpy.array(nf).tofile(con)

	finfo = True
	ginfo = True
	ffsp = list()

	if prof: t1 = t2 = t3 = 0
	for j in fnoms.argsort():
		print("... field",j,fnoms[j])
		ff = fields[j]

		if a.format == "GRIB":
			geom = ff.geometry
			base = ff.validity.getbasis()
			step = ff.validity.term()
			if ggeom is None:
				ggeom = geom
				gbase = base
			elif geom is not ggeom or base != gbase:
				next

			dims = (geom.dimensions["X"],geom.dimensions["Y"],0,0,1)
			if ginfo: print(dims)

			bb = numpy.array(base.timestamp())
			ss = numpy.array(step.total_seconds())
			if con != None:
				if ginfo:
					dims.tofile(con)
					bb.tofile(con)
				ss.tofile(con)

			ginfo = False

		if ff.spectral:
			#ffsp.append(ff)
			#continue

			if prof: t = os.times().elapsed
			splen = len(ff.data)
			spsize = ff.data.size
			ff.sp2gp()
			if finfo:
				if prof: t1 = t1+(os.times().elapsed-t)
				print("Spectral/gridpoint sizes:",splen,spsize,ff.data.size)
			if prof: t2 = t2+(os.times().elapsed-t)
		else:
			if finfo:
				print("Gridpoint size:",ff.data.size)

		data = ff.data
		if (not islam):
			data = data.compressed()
			if finfo:
				print("Reduced gridpoints:",data.size)
				finfo = False

		if con != None:
			if prof: t = os.times().elapsed
			data.tofile(con)
			if prof: t3 = t3+(os.times().elapsed-t)

	# (debug) mark end of data with one known integer value
	if con != None:
		if len(ffsp) > 0:
			#exe = cf.ThreadPoolExecutor(max_workers=4)
			if prof: t = os.times().elapsed
			for ff in ffsp:
				splen = len(ff.data)
				ff.sp2gp()
				if finfo:
					if prof: t1 = t1+(os.times().elapsed-t)
					print("Spectral/gridpoint sizes:",splen,ff.data.size,"- fields:",len(ffsp))
					finfo = False
				#f = exe.submit(ff.sp2gp)
				#f.result()

			#with exe:
			#	f = exe.submit(sp2gp,ffsp)
			#	f.result()
			if prof: t2 = t2+(os.times().elapsed-t)

			for ff in ffsp:
				data = ff.data
				if (not islam): data = data.compressed()
				if prof: t = os.times().elapsed
				data.tofile(con)
				if prof: t3 = t3+(os.times().elapsed-t)

		numpy.array(nf).tofile(con)

if con != None:
	con.close()

if prof: print("times sp2gp(1)/sp2gp/tofile:",round(t1,3),round(t2,3),round(t3,3))

