#!/bin/sh -e

URL='https://cfd.ku.edu/hiocfd/data.tgz'

wget "$URL" || curl -O "$URL"
tar -xzf data.tgz spectral_Re1600_512.gdiag
