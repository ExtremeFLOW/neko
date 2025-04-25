#!/bin/bash

# Extract time and enstrophy from neko output from argument $1

echo "# Time Enstrophy"
awk '/enst:/ {print($3, $NF)}' $1
