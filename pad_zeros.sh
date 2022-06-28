#!/bin/bash

str0=0
str00=00

for FILE in *; do
  # echo $FILE
  size=${#FILE}
  if [ ${FILE: -3:3} == "txt" ] || [ ${FILE: -3:3} == "pdf" ]
  then
    if [[ "$FILE" == "plotspecemission_"* ]]
    then
      if [ $size == 26 ]
      then
        mv $FILE ${FILE:0:17}$str00${FILE:17:3}$str00${FILE:20:6}
      fi
      if [ $size == 28 ]
      then
        mv $FILE ${FILE:0:17}$str0${FILE:17:4}$str0${FILE:21:7}
      fi
      if [ $size == 29 ]
      then
        mv $FILE ${FILE:0:17}$str0${FILE:17:12}
      fi
    fi
    if [[ "$FILE" == "plotspec_"* ]]
    then
      if [ $size == 18 ]
      then
        mv $FILE ${FILE:0:9}$str00${FILE:9:3}$str00${FILE:12:6}
      fi
      if [ $size == 20 ]
      then
        mv $FILE ${FILE:0:9}$str0${FILE:9:4}$str0${FILE:13:7}
      fi
      if [ $size == 21 ]
      then
        mv $FILE ${FILE:0:9}$str0${FILE:9:12}
      fi
    fi
  fi
done


# for i in {1..9}; do
#   mv $str1$i$str2$i$str3 $str1$str00$i$str2$i$str3
#   mv $str1$i$str2$i$str4 $str1$str00$i$str2$i$str4
# done
#
# for i in {10..99}; do
#   mv $str1$i$str2$i$str3 $str1$str0$i$str2$i$str3
#   mv $str1$i$str2$i$str4 $str1$str0$i$str2$i$str4
# done
