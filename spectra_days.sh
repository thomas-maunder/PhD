#!/bin/bash

for i in {1..100}; do
	plotartisspectrum -timedays $i --write_data
done


if ls spectrum_ts*.txt >/dev/null 2>&1
then
  mv spectrum_ts*.txt spectra/total_spectrum
fi

if ls plotspec_*.pdf >/dev/null 2>&1
then
  mv plotspec_*.pdf spectra/total_spectrum
fi

if ls plotspec_*.txt >/dev/null 2>&1
then
  mv plotspec_*.txt spectra/total_spectrum
fi

[ -f "time_list.txt" ] && mv time_list.txt spectra/total_spectrum
[ -f "spectra_list.txt" ] && mv spectra_list.txt spectra/total_spectrum


for i in {1..100}; do
	plotartisspectrum -timedays $i -plotviewingangle 0 50 90 --write_data
done

[ -f "time_list.txt" ] && mv time_list.txt spectra/total_spectrum
[ -f "spectra_list.txt" ] && mv spectra_list.txt spectra/total_spectrum

if ls spectrum_ts*.txt >/dev/null 2>&1
then
  mv spectrum_ts*.txt spectra/viewing_angle/bins_0_50_90
fi

if ls plotspec_*.pdf >/dev/null 2>&1
then
  mv plotspec_*.pdf spectra/viewing_angle/bins_0_50_90
fi

if ls plotspec_*.txt >/dev/null 2>&1
then
  mv plotspec_*.txt spectra/viewing_angle/bins_0_50_90
fi
