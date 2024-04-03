MAKEFLAGS += --no-print-directory

# ne pas mettre ~ pour P : il faut un chemin absolu
P = $(HOME)/proc/diags
B = ~/bin

.PHONY: build install diags

build:
	# rien

install:
	! git status --porcelain 2>/dev/null | grep -qvE "^\?\? "
	make diags
	make $B/epy_dump.py
	make $B/diag.sh
	make $B/obstat.sh
	if git status >/dev/null 2>&1; then \
		grep -q $(shell git log -1 --pretty=format:%h 2>/dev/null) $P/version || \
			git log -1 --oneline >> $P/version; \
	fi

diags:
	mkdir -p $P
	cp -pruv diag.R obstat.R img.html obstat.html $P

$B/epy_dump.py: epy_dump.py
	cp -uv epy_dump.py $B

$B/diag.sh: diag.sh
	sed -re "s:diags=.+:diags=$P:" diag.sh > $B/diag.sh
	chmod a+x $B/diag.sh

$B/obstat.sh: obstat.sh
	sed -re "s:diags=.+:diags=$P:" obstat.sh > $B/obstat.sh
	chmod a+x $B/obstat.sh
