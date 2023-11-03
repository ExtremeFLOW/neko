#!/bin/bash

# Extract time and enstrophy from neko output on stdin

echo "# Time Enstrophy"
awk '/enst:/ {print($3, $NF)}' $1
