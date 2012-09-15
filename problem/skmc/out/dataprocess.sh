#!/bin/bash
echo data process begin...

rm slice*.png

OF=movie-$(date +%d%m%y-%k).mp4
echo $OF

R CMD BATCH myplot.r /dev/tty
ffmpeg -y -qscale 1 -r 1 -b 9600 -i slice%5d.png $OF

echo ... data process end