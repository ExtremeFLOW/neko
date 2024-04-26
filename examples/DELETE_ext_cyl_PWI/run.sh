#!/usr/bin/bash

[ -f neko ] && neko=./neko || neko=neko

if [ ! -x $neko ]; then
    echo -e "Neko not found." >&2
    echo -e "Please ensure Neko is installed and in your PATH." >&2
    echo -e "Alternatively, set the NEKO_DIR environment variable." >&2
    exit 1
fi

mpirun -n 8 $neko ext_cyl.case >ext_cyl.log
