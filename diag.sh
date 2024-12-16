#!/bin/sh

diag=~/util/diag

usage()
{
	printf "
Description:
	Produce a list of diagnostics on model forecasts

Usage:
	diag.sh [-conf CONFIG] [-o PNG] [-ref PATH] [-h]

Arguments:
	CONFIG: path to a directory containing files for settings: param.txt, domain.txt and \
file.txt as mandatory files, and optionally level.txt and date.txt (see details)
	PNG: path to output files (graphics + some HTML files)
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
mapstat=1
loc=""
ref=.
cmp=""

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
		-cmp)
			cmp=$2
			shift
			;;
		-o)
			loc=$2
			shift
			;;
		-nograph)
			graph=0
			;;
		-nomapstat)
			mapstat=0
			;;
		*)
			echo "Warning: unknown option '$1', ignored" >&2
			;;
	esac

	shift
done

if [ -z "$conf" -o -z "$ref" ]
then
	echo "Error: input options missing
conf: '$conf'
ref: '$ref'" >&2
	exit 1
fi

set -e

if [ -n "$ref" ]
then
	ls -d $ref > /dev/null
fi

opt=""
if [ -n "$cmp" ]
then
	opt="cmp=$cmp"
fi

if [ -n "$loc" ]
then
	mkdir -p $loc
else
	loc=$(mktemp -d diagXXX)
	echo "--> sending output to $loc"
fi

type R > /dev/null 2>&1 || module -s load intel R > /dev/null 2>&1

if [ $graph -eq 1 ]
then
	echo "Production of maps and figures for data and errors"
	R --slave -f $diag/diag.R --args png=$loc ref=$ref mapstat=$mapstat $opt > $loc/diag.log
fi

echo "Production of scores"
R --slave -f $diag/score.R --args png=$loc $opt > $loc/score.log

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

temph=$(mktemp -d htmlXXX)
trap 'rm -r $temph' 0

while read -a tt
do
	echo ${tt[*]} | grep -q " TRUE" && printf "\t<option>%s</option>\n" "${tt[*]}"
done < $loc/steps.txt > $temph/steps.html

for par in $params
do
	ficpar=$(dirname $loc)/$par.html
	echo ". creating $ficpar"

	rm -f $temph/*_$par.html

	echo ".. maps, hists and options (maph.html)"
	it=0
	while read -a tt
	do
		it=$((it+1))
		echo ${tt[*]} | grep -qE ' graph: TRUE' || continue

		base=$(echo ${tt[*]} | sed -re 's:base/step\: ([12][0-9]+) R([0-9]+) (\+[0-9]+\w*) .+:\1\2\3:')
		step=$(echo ${tt[*]} | sed -re 's:base/step\: [12][0-9]+ R[0-9]+ \+([0-9]+\w*) .+:\1:')
		title=$(echo ${tt[*]} | sed -re 's: \- graph.+::')

		for dom in $doms
		do
			echo "<table>"
			echo "<tr><th>Maps, selection of levels</th><th>Cross-sections and diagrams</th>"
			att="colspan='2'"
			if false && [ -n "$ref" ]
			then
				echo "<th>Maps for reference</th><th>Cross-sections and diagrams for ref</th>"
				att="colspan='4'"
			fi

			echo "</tr>"

			# rows with name attribute (fig for images, title for row header)
			echo "<tr><th name='title' $att>Domain $dom - $title</th></tr>"

			for suf in "" ref diff
			do
				echo "<tr>"

				for typ in map hist
				do
					fic=$loc/$typ$suf$base${dom}_$par.png
					[ -s $fic ] || continue

					printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
					printf "\t<option>%s</option>\n" $fic >> $temph/$typ$suf${dom}_$par.html
				done

				echo "</tr>"
			done

			cat $diag/step.html

			[ $mapstat -eq 1 ] || continue

			{
			echo "<tr><th $att>Min/max of upper and lower levels - Domain $dom</th></tr>"

			for suf in "" ref
			do
				echo "<tr>"

				for typ in mapn mapx
				do
					fic=$loc/$typ$suf$base${dom}_$par.png
					[ -s $fic ] || continue

					printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
					printf "\t<option>%s</option>\n" $fic >> $temph/$typ$suf${dom}_$par.html
				done

				echo "</tr>"
			done
			} > $temph/mnx.html

			grep -qE '<img .+ src=' $temph/mnx.html && cat $temph/mnx.html

			{
			for stat in bias rmse errx dayx
			do
				echo "<tr><th colspan='2'>Param '$par' $stat - Domain $dom</th></tr>"
				echo "<tr>"

				for typ in map$stat hist$stat
				do
					fic=$loc/$typ$step${dom}_$par.png
					[ -s $fic ] || continue

					printf "\t<td><img name='fig' src='%s' alt='missing image'/></td>\n" $fic
					printf "\t<option>%s</option>\n" $fic >> $temph/$typ${dom}_$par.html
				done

				echo "</tr>"
			done
			} > $temph/stat.html

			grep -qE '<img .+ src=' $temph/stat.html && cat $temph/stat.html
		done > $temph/idom.html

		if grep -qE '<img .+ src=' $temph/idom.html && [ ! -s $temph/idom.html.save ]
		then
			cat $temph/idom.html
			cp $temph/idom.html $temph/idom.html.save
			echo "</table>"
		fi
	done < $loc/steps.txt > $temph/maph.html

	rm -f $temph/idom.html.save
	echo ".. HTML select (stat.html)" # with name attributes step and map
	{
	sed -re 's:TAG NAME:step:' -e "/TAG OPT/r $temph/steps.html" $diag/select.html
	for dom in $doms
	do
		# same order as table images
		for typ in map hist mapref histref mapdiff histdiff mapn mapx mapnref mapxref \
			mapbias histbias maprmse histrmse maperrx histerrx mapdayx histdayx
		do
			fic=$temph/$typ${dom}_$par.html
			[ -s $fic ] || continue

			sed -re 's:TAG NAME:map:' -e "/TAG OPT/r $fic" $diag/select.html
		done
	done
	} > $temph/stat.html

	echo ".. Stats and scores (score.html)"
	typd=(stat statv err errv score scorev scorez scorevz scoret rmsevt)
	titd=("Statistics of forecast" "Statistics of profile" "Statistics of forecast error" \
		"Statistics of profile error" "Scores of forecasts" "Scores of profile" \
		"Scores of zonal mean" "Score of zonal mean - profile" "Daily scores" "Daily RMSE")
	i=0
	while [ $i -lt ${#typd[*]} ]
	do
		echo ${typd[i]} | grep -q score && th="<p>Ref data (blue): $cmp</p>" || th=""
		[ "$cmp" ] || th=""
		echo "<h2>${titd[i]} on domains</h2>
$th
<table>
<tr>"
		n=0
		for ficp in $(ls -1 $loc | grep -E "${typd[i]}[0-9]+_$par.png" | sort)
		do
			if [ $n -eq 2 ]
			then
				echo -e "</tr>\n<tr>"
				n=0
			fi

			printf "\t<td><img src=\"%s\" alt=\"missing image\"/></td>\n" $loc/$ficp
			n=$((n+1))
		done

		echo "</tr>
</table>"
		i=$((i+1))
	done > $temph/score.html

	echo ".. fill par template"
	sed -re "s:TAG PAR:$par:" -e "s:TAG REF:$ref:" -e "/TAG MAP/r $temph/maph.html" \
		-e "/TAG STAT/r $temph/stat.html" -e "/TAG SCORE/r $temph/score.html" \
		$diag/par.html > $ficpar
done

echo ". creating params.html"
typd=(score scorez)
titd=("Scores of forecasts" "Scores of zonal mean")
i=0
while [ $i -lt ${#typd[*]} ]
do
	echo "<h2>${titd[i]} on domains</h2>
<table>"

	for par in $params
	do
		echo "<tr>"
		for ficp in $(ls -1 $loc | grep -E "${typd[i]}[0-9]+_$par.png" | sort)
		do
			printf "\t<td><img src=\"%s\" alt=\"missing image\"/></td>\n" $loc/$ficp
		done

		echo "</tr>"
   done

	echo "</table>"
	i=$((i+1))
done > $temph/score.html

sed -re "s:TAG PAR:multi parameter:" -e "/TAG SCORE/r $temph/score.html" $diag/par.html \
	> params.html
