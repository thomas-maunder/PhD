#!/bin/bash

echo Starting analysis of ARTIS run "$PWD"...


# Remove files that shouldn't be there yet if they exist
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

echo Generating directory structure...
[ ! -d "./emission_absorption" ] && mkdir emission_absorption || echo "Directory emission_absorption already exists!"
[ ! -d "./emission_absorption/default" ] && mkdir emission_absorption/default || echo "Directory emission_absorption/default already exists!"
[ ! -d "./emission_absorption/fixed_ion_list" ] && mkdir emission_absorption/fixed_ion_list || echo "Directory emission_absorption/fixed_ion_list already exists!"
[ ! -d "./emission_absorption/fixed_ion_list/list1" ] && mkdir emission_absorption/fixed_ion_list/list1 || echo "Directory emission_absorption/fixed_ion_list/list1 already exists!"
[ ! -d "./emission_absorption/fixed_ion_list/list1/images" ] && mkdir emission_absorption/fixed_ion_list/list1/images || echo "Directory emission_absorption/fixed_ion_list/list1/images already exists!"
[ ! -d "./spectra" ] && mkdir spectra || echo "Directory spectra already exists!"
[ ! -d "./spectra/viewing_angle" ] && mkdir spectra/viewing_angle || echo "Directory spectra/viewing_angle already exists!"
[ ! -d "./spectra/viewing_angle/bins_0_50_90" ] && mkdir spectra/viewing_angle/bins_0_50_90 || echo "Directory spectra/viewing_angle/bins_0_50_90 already exists!"
[ ! -d "./spectra/total_spectrum" ] && mkdir spectra/total_spectrum || echo "Directory spectra/total_spectrum already exists!"

echo Copying analysis files into relevant directories...
[ ! -f "emission_absorption/fixed_ion_list/list1/pad_zeros.sh" ] && cp pad_zeros.sh emission_absorption/fixed_ion_list/list1 || echo "File pad_zeros.sh is already in emission_absorption/fixed_ion_list/list1"
[ ! -f "emission_absorption/default/pad_zeros.sh" ] && cp pad_zeros.sh emission_absorption/default || echo "File pad_zeros.sh is already in emission_absorption/default"
[ ! -f "emission_absorption/fixed_ion_list/list1/images/create_videos.sh" ] && cp create_videos.sh emission_absorption/fixed_ion_list/list1/images || echo "File pad_zeros.sh is already in emission_absorption/fixed_ion_list/list1/images"
[ ! -f "spectra/total_spectrum/pad_zeros.sh" ] && cp pad_zeros.sh spectra/total_spectrum || echo "File pad_zeros.sh is already in spectra/total_spectrum"
[ ! -f "spectra/viewing_angle/bins_0_50_90/pad_zeros.sh" ] && cp pad_zeros.sh spectra/viewing_angle/bins_0_50_90 || echo "File pad_zeros.sh is already in spectra/viewing_angle/bins_0_50_90"

echo Folder structure created and analysis ready to begin!

echo Plotting light curves and determining peak light...

plotartislightcurve
plotartislightcurve --print_data > lightcurve.txt
python peak_light.py > peak_light.txt

echo Plotting band magnitudes...

# USE THIS IF YOU KNOW THE CORRECT PLOT BOUNDS YOU WANT - CHANGE xmax, xmin,
# ymax, ymin AS REQUIRED, OR USE THE NEXT BASH COMMAND TO CHECK AND CHANGE.

plotartislightcurve --magnitude -filter U B V R I zs -ymax -16.5 -ymin -10

# OPTIONAL IF YOU DON'T KNOW WHAT THE PLOT BOUNDS SHOULD BE - ONCE YOU KNOW WHAT
# THEY SHOULD BE YOU SHOULD JUST MANUALLY ADD THEM TO THIS SCRIPT.

# bash band_lightcurve_loop.sh

echo Plotting colour evolution...

plotartislightcurve -colour_evolution B-V
plotartislightcurve -colour_evolution V-R

echo Light curves completed.

echo Beginning spectra analysis...

bash spectra_days.sh
bash emission_absorption_days.sh

# Pad zeros into file names
bash spectra/total_spectrum/pad_zeros.sh
bash spectra/viewing_angle/bins_0_50_90/pad_zeros.sh
bash emission_absorption/default/pad_zeros.sh
bash emission_absorption/fixed_ion_list/list1/pad_zeros.sh

echo Generating video of spectral evolution for a fixed ion list1...

python plot_emission_absorption_auto.py
mv emission*png emission_absorption/fixed_ion_list/list1/images
mv eission*pdf emission_absorption/fixed_ion_list/list1
cd emission_absorption/fixed_ion_list/list1/images
bash create_videos.sh

echo Video created and saved in emission_absorption/fixed_ion_list/list1/images as spectral_evolution.mp4

cd ../../../../

echo FINISHED ARTIS ANALYSIS!
