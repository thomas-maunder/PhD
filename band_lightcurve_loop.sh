#!/bin/bash

plotartislightcurve --magnitude -filter U B V R I zs

xdg-open plotlightcurves.pdf

echo "Do you need to adjust any ranges? (y/n)"

read answer

# echo $answer
if [[ $answer == y ]]
then
  var=y
  echo $var
  while [[ $var == y ]]
  do
    echo "xmax? (y/n)"
    read answer
    if [[ $answer == y ]]
    then
      echo enter value
      read xmax
    else
      xmax=100
    fi
    echo "xmin? (y/n)"
    read answer
    if [[ $answer == y ]]
    then
      echo enter value
      read xmin
    else
      xmin=0
    fi
    echo "ymax? (y/n)"
    read answer
    if [[ $answer == y ]]
    then
      echo enter value
      read ymax
    else
      ymax=-20
    fi
    echo "ymin? (y/n)"
    read answer
    if [[ $answer == y ]]
    then
      echo enter value
      read ymin
    else
      ymin=-14
    fi

  plotartislightcurve --magnitude -filter U B V R I zs -xmax $xmax -ymax $ymax -xmin $xmin -ymin $ymin

  echo "Do you need to adjust any ranges again? (y/n)"
  read var
  done

fi
