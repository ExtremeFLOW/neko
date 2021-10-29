#!/bin/bash
IN=$(basename $1)
NAME=${IN%.cl}

echo "const char *"${NAME}" =" > $1.h
sed -e 's/\\/\\\\/g;s/"/\\"/g;s/^/"/;s/$/\\n"/' \
    $1 >> $1.h
echo ";" >>$1.h
