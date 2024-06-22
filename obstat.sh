#!/bin/sh

diags=~/util/diag

set -e

usage()
{
	echo "
Description:
	Produce HTML files gathereing a few plots from an Obstat file

Usage:
	obsta.sh FOBSTAT [-o PNG] [-h]

Options:
	FOBSTAT: input Obstat text file
	PNG: local directory name (no '/') where output PNG files (+ text files) are produced
	-h: print this help and exit normally

Details:
	Input file is an Obstat file, produced by some NWP application.
	HTML output files are produced, containing text and images built from the Obstat \
file. Files are named after [param].html' with param an Obstat parameter.
	Plots are placed into directory [PNG] (created if needed) aside (ie local) to the \
HTML file."
}

if [ $# -eq 0 ] || echo " $*" | grep -qE '\-h'
then
	usage
	exit
fi

fin=""
png=""

while [ $# -ne 0 ]
do
	case $1 in
	-o)
		png=$2
		shift
		;;
	*)
		if [ -n "$fin" ]
		then
			echo "Error: usage" >&2
			exit 1
		fi

		fin=$1
		;;
	esac

	shift
done

if [ -z "$fin" -o -z "$png" ]
then
	echo "Error: mandatory arguments missing
fin: '$fin'
png: '$png'" >&2
	exit 1
fi

set -e

ls -L $fin > /dev/null
file -L $fin | grep -q text
if ! grep -qi 'BEGIN STATDEF' $fin
then
	echo "Error: $fin is not an Obstat file (no pattern 'BEGIN STATDEF')" >&2
	exit 1
fi

echo "Create graphics from obstat file $fin in $png"
mkdir -p $png

type R >/dev/null 2>&1 || module -s load intel R >/dev/null 2>&1

R --slave -f $diags/obstat.R --args $fin $png

groups=$(ls -1 $png | sed -re 's:(\w+)-.+:\1:' | sort -u | xargs)
temp=$(mktemp -d -t obsXXX)
trap 'rm -r $temp' EXIT

for grp in $groups
do
	echo ". group $grp.html"

	re="s:$grp-(.+)_.+:\1:"
	for par in $(ls -1 $png | grep "^$grp-" | sed -re "$re" | sort -u | xargs)
	do
		echo -e "<h3>Param $par\n<table><tr>"

		np=0
		for fic in $(find $png -name $grp-${par}_*.png | sort)
		do
			[ $np -gt 0 -a $((np%2)) -eq 0 ] && echo -e "</tr>\n<tr>"
			np=$((np+1))
			printf "\t<td><img src='%s' alt='missing image'/></td>\n" $fic
		done

		echo "</tr></table>"
	done > $temp/$grp.html

	sed -re "s:TAG GROUP:$grp:" -e "/TAG IMG/r $temp/$grp.html" $diags/obstat.html > \
		$grp.html
done
