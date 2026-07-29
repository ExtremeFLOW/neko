#!/bin/sh

echo "Updating configuration..."
echo "Running libtoolize"
if command -v libtoolize >/dev/null 2>&1; then
    libtoolize --install --force
elif command -v glibtoolize >/dev/null 2>&1; then
    glibtoolize --install --force
else
    echo "No libtoolize or glibtoolize found on your system" >&2
    exit 1
fi

rm -rf autom4te.cache
rm -f aclocal.m4

echo "Running aclocal"
if which aclocal > /dev/null 2>&1; then
    aclocal -I m4 --install
else
    echo "No aclocal found on your system" >&2
    exit 1
fi

echo "Running autoconf"
autoconf

echo "Running automake"
automake --add-missing

echo "Deleting autom4te.cache directory"
rm -rf autom4te.cache