#!/bin/bash

[ ! -d './emission_absorption/ions' ] && mkdir emission_absorption/ions

ion='He_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'He I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='He_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'He II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='He_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'He III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='C_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'C I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='C_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'C II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='C_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'C III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='N_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'N I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='N_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'N II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='N_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'N III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='O_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'O I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='O_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'O II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='O_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'O III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='F_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'F I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='F_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'F II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='F_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'F III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ne_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ne I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ne_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ne II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ne_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ne III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Na_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Na I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Na_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Na II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Na_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Na III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mg_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mg I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mg_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mg II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mg_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mg III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Al_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Al I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Al_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Al II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Al_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Al III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Si_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Si I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Si_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Si II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Si_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Si III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='P_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'P I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='P_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'P II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='P_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'P III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='S_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'S I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='S_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'S II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='S_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'S III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cl_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cl I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cl_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cl II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cl_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cl III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ar_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ar I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ar_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ar II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ar_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ar III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='K_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'K I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='K_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'K II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='K_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'K III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ca_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ca I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ca_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ca II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ca_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ca III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Sc_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Sc I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Sc_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Sc II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Sc_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Sc III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ti_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ti I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ti_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ti II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ti_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ti III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='V_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'V I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='V_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'V II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='V_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'V III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cr_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cr I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cr_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cr II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cr_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cr III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mn_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mn I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mn_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mn II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Mn_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Mn III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Fe_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Fe I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Fe_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Fe II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Fe_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Fe III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Co_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Co I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Co_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Co II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Co_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Co III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ni_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ni I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ni_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ni II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Ni_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Ni III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cu_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cu I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cu_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cu II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Cu_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Cu III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Zn_I'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Zn I' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Zn_II'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Zn II' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../

ion='Zn_III'
for i in {2..100}; do
  echo $ion $i
  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist 'Zn III' --write_data
  rm plotspecemission_*pdf
  target=${ion/' '/'_'}
  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target
  mv plotspecemission_*txt emission_absorption/ions/$target
done
cp pad_zeros.sh emission_absorption/ions/$target
cd emission_absorption/ions/$target
bash pad_zeros.sh
cd ../../../
