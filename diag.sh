#!/bin/sh

diags=~/util/diags

usage()
{
	printf "
Description:
	Produce a list of diagnostics on model forecasts

Usage:
	diag.sh [-conf CONFIG] [-o HTML] [-ref PATH] [-h]

Arguments:
	CONFIG: path to a directory containing files for settings: param.txt, domain.txt and \
file.txt as mandatory files, and optionally level.txt and date.txt (see details)
	HTML: path to output files (graphics and HTML files)
	PATH: path where to find reference files
	-h: print this help and exits normally

Details:
	Default value for CONFIG is the local directory config/. This directory must contain \
files param.txt, domain.txt and file.txt, respectively indicating parameters, domains \
and files to target for graphical and statistical diagnostics. Optionnally, it can also \
contain files level.txt and date.txt, indicating levels to filter (given as indices in \
the range [1-NFLEVG]) and base dates to read for files. If date.txt is present, paths to \
files given in file.txt are prepended by 'HH/YYYYMMDD/', meaning that files are to be \
located in directory hierarchy with time (HH) and date (YYYYMMDD) information.
	PATH is the path to a hierarchy of files and directories following \
'HH/YYYYMMDD/filename', where HH is a reference to some 'base time' in 00-24 format, \
YYYYMMDD is the 'base date' and filename is the name of a file (with its extension, \
if any). The 'base date/time' is the moment from which the product is considered to be \
existing (ie it can be used).

Dependencies:
	R software
"
}

if [ $# -eq 0 ] || echo " $*" | grep -qE ' \-h\>'
then
	usage
	exit
fi

conf="config"
graph=1
html=""

while [ $# -ne 0 ]
do
	case $1 in
		-conf)
			conf=$2
			shift
			;;
		-ref)
			ref=$2
			shift
			;;
		-o)
			html=$2
			shift
			;;
		-nograph)
			graph=0
			;;
		*)
			echo "Warning: unknown option '$1', ignored" >&2
			;;
	esac

	shift
done

if [ -n "$conf" ]
then
	fic=$conf/file.txt
	par=$conf/param.txt
	[ -s $conf/level.txt ] && lev=$conf/level.txt || lev=0
fi

if [ -z "$fic" -o -z "$lev" -o -z "$par" ]
then
	echo "Error: input options missing
conf: '$conf'
lev: '$lev'
par: '$par'
fic: '$fic'" >&2
	exit 1
fi

if [ ! -s $lev ] && ! echo $lev | grep -qE '[0-9]+(:[0-9]+)*'
then
	echo "Error: option '-level' uncorrectly defined" >&2
	exit 1
	#echo $lev | tr ':' '\n' > $loc/levels.txt
fi

set -e

ls $par > /dev/null
if [ -n "$ref" ]
then
	ls -d $ref > /dev/null
	opt="ref=$ref"
fi

if [ -n "$html" ]
then
	loc=$html
else
	loc=$(mktemp -d diagXXX)
fi

echo "--> sending output to $loc"
mkdir -p $loc

if [ $graph -eq 1 ]
then
	type R > /dev/null 2>&1 || module -s load intel R > /dev/null 2>&1
	R --slave -f $diags/diag.R --args fic=$fic params=$par level=$lev png=$loc $opt > \
		$loc/diag.log
fi

[ -n "$html" ] || exit

ficdom=$conf/domain.txt
if [ ! -s $ficdom ]
then
	echo "Error: no domains file '$ficdom'" >&2
	exit 1
fi

doms=$(grep -E '^ *\w+' $ficdom | awk 'NR > 1 {print $1}' | xargs)
params=$(ls -1 $loc | grep -E 'map[1-9].+_\w+\.png$' | sed -re 's:.+_(\w+)\.png:\1:' | \
	sort -u | xargs)
echo "HTML files: $(cat $loc/steps.txt | wc -l) forecast steps
domains: $doms"

for par in $params
do
	ficpar=$(dirname $html)/$par.html
	echo ". creating $ficpar"

	{
	it=0
	while read -a tt
	do
		it=$((it+1))
		echo ${tt[*]} | grep -qE 'graph: TRUE' || continue

		echo "<table>"
		echo "<tr><th>Maps, selection of levels</th><th>Cross-sections and diagrams</th>"
		att="colspan='2'"
		if [ -n "$ref" ]
		then
			echo "<th>Maps for reference</th><th>Cross-sections and diagrams for ref</th>"
			att="colspan='4'"
		fi

		echo "</tr>"

		for dom in $doms
		do
			echo "<tr><th name='title' $att>Param '$par' - ${tt[*]} - Domain $dom</th></tr>"
			echo "<tr>"

			for typ in map hist mapdiff histdiff
			do
				fic=$loc/$typ$it${dom}_$par.png
				[ -s $fic ] || continue

				printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
			done

			echo "</tr>"

			cat $diags/step.html

			for stat in bias rmse errx
			do
				echo "<tr><th colspan='2'>Param '$par' $stat - Domain $dom</th></tr>"
				echo "<tr>"

				for typ in map$stat hist$stat
				do
					fic=$loc/$typ$it${dom}_$par.png
					[ -s $fic ] || continue

					printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
				done

				echo "</tr>"
			done

			echo "<tr><th colspan='2'>Param '$par' min/max - Domain $dom</th></tr>"
			echo "<tr>"

			for typ in mapn mapx
			do
				fic=$loc/$typ$it${dom}_$par.png
				[ -s $fic ] || continue

				printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
			done

			echo "</tr>"
		done > $loc/idom.html

		if grep -qE '<img .+ src=' $loc/idom.html
		then
			cat $loc/idom.html
			echo "</table>"
			break
		fi

		echo "</table>"
	done < $loc/steps.txt

	sed -re 's:TAG NAME:step:' -e "/TAG OPT/r $loc/steps.html" $diags/select.html
	for dom in $doms
	do
		for stat in "" diff bias rmse errx
		do
			for typ in map$stat hist$stat
			do
				fic=$loc/$typ${dom}_$par.html
				[ -s $fic ] || continue

				sed -re 's:TAG NAME:map:' -e "/TAG OPT/r $fic" $diags/select.html
			done
		done
	done

	typd=("stat" "statv" "err" "errv" "score" "scorev" "scoret" "rmsevt")
	titd=("Statistics of forecast" "Statistics of profile" "Statistics of forecast error" \
		"Statistics of profile error" "Scores of forecasts" "Scores of profile" \
		"Daily scores" "Daily RMSE")
	i=0
	while [ $i -lt ${#typd[*]} ]
	do
		echo "<h2>${titd[i]} on domains</h2>
<table>
<tr>"

		for ficp in $(ls -1 $loc | grep -E "${typd[i]}[0-9]+_$par.png" | sort)
		do
			printf "\t<td><img src=\"%s\" alt=\"missing image\"/></td>\n" $loc/$ficp
		done

		echo "</tr>
</table>"
		i=$((i+1))
	done
	} > $loc/$par.html

	sed -re "s:TAG PAR:$par:" -e "/TAG BODY/r $loc/$par.html" $diags/par.html > $ficpar
done
