#!/bin/bash
echo data process begin...
echo

rm slice*.png

OF=movie-$(date +%M%H-%d%m%y).mp4

R CMD BATCH --slave myplot.r /dev/tty

echo
echo

ffmpeg -y -qscale 5 -r 5 -b 9600 -loglevel 0 -i slice%5d.png $OF

echo
echo
echo $OF is generated.
echo done!
echo
