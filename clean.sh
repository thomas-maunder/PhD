#!/bin/bash

if ls plot*txt >/dev/null 2>&1
then
  rm plot*txt
fi

if ls plot*pdf >/dev/null 2>&1
then
  rm plot*pdf
fi

if ls peak*txt >/dev/null 2>&1
then
  rm peak*txt
fi

[ -f lightcurve.txt ] && rm lightcurve.txt

[ -d "./emission_absorption" ] && rm -dr emission_absorption
[ -d "./spectra" ] && rm -dr spectra
