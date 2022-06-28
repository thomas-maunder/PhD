#!/bin/bash

# file_ext1=d.txt
# file_ext2=d.pdf
# mid_str=d_
start_str=plotspecemission_
# buffer_str=0
# bound=10

rm em_ab_default_output.txt
touch em_ab_default_output.txt

for i in {1..100}; do
	plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype --notitle >> em_ab_default_output.txt
	# if [ "$i" -lt "$bound" ]
	# then
		# echo $start_str$i$mid_str$i$file_ext2
		# echo $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext2
		# mv $start_str$i$mid_str$i$file_ext2 $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext2
		# mv $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext2 emission_absorption/default
		# mv $start_str*$file_ext2 emission_absorption/default
	# fi
done
mv $start_str* emission_absorption/default
mv em_ab_default_output.txt emission_absorption/default

rm em_ab_fixed_output.txt
touch em_ab_fixed_output.txt

for i in {1..100}; do
	artistools-spectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'free-free' 'Fe II' 'V I' 'Co II' 'Ti II' 'V II' 'Mg II' 'Cr II' 'C I' 'Co I' 'Ca II' 'Ti I' 'Ca II' 'Mg I' 'Other' --write_data >> em_ab_fixed_output.txt
	# if [ "$i" -lt "$bound" ]
	# then
		# echo $start_str$i$mid_str$i$file_ext1
		# echo $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext1
		# mv $start_str$i$mid_str$i$file_ext1 $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext1
		# mv $start_str$buffer_str$i$mid_str$buffer_str$i$file_ext1 emission_absorption/fixed_ion_list/list1
		# mv $start_str*$file_ext1 emission_absorption/fixed_ion_list/list1
	# fi
done
mv $start_str* emission_absorption/fixed_ion_list/list1
mv em_ab_fixed_output.txt emission_absorption/fixed_ion_list/list1
