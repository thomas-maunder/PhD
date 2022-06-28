#!/bin/bash

# ffmpeg -r:v 1/2 -i "ea_lists.txt" -r:v 30 -codec:v libx264  -preset veryslow -pix_fmt yuv420p -crf 28 -an "test.mp4"

touch ea_lists.txt
rm ea_lists.txt
touch ea_lists.txt

for FILE in *; do
  if [ ${FILE: -3:3} == "png" ]
  then
    echo file $FILE >> ea_lists.txt
    echo duration 0.01 >> ea_lists.txt
  fi
done

echo ea_lists

ffmpeg -y -r:v 1 -f concat -safe 0 -i "ea_lists.txt" -r:v 30 -c:v libx264 -preset veryslow -vb 30M -vf "setpts=1*PTS" "spectral_evolution.mp4" -hide_banner
