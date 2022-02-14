#!/bin/bash

# Extract time and enstrophy from neko output on stdin

echo "# Time Enstrophy"
grep -a Enstrophy | cut -d" " -f3,8
