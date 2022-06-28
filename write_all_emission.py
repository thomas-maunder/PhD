ions = ['He_I', 'He_II', 'He_III', 'C_I', 'C_II', 'C_III', 'N_I', 'N_II', 'N_III', 'O_I', 'O_II', 'O_III', 'F_I', 'F_II', 'F_III', 'Ne_I', 'Ne_II', 'Ne_III', 'Na_I', 'Na_II', 'Na_III', 'Mg_I', 'Mg_II', 'Mg_III', 'Al_I', 'Al_II', 'Al_III', 'Si_I', 'Si_II', 'Si_III', 'P_I', 'P_II', 'P_III', 'S_I', 'S_II', 'S_III', 'Cl_I', 'Cl_II', 'Cl_III', 'Ar_I', 'Ar_II', 'Ar_III', 'K_I', 'K_II', 'K_III', 'Ca_I', 'Ca_II', 'Ca_III', 'Sc_I', 'Sc_II', 'Sc_III', 'Ti_I', 'Ti_II', 'Ti_III', 'V_I', 'V_II', 'V_III', 'Cr_I', 'Cr_II', 'Cr_III', 'Mn_I', 'Mn_II', 'Mn_III', 'Fe_I', 'Fe_II', 'Fe_III', 'Co_I', 'Co_II', 'Co_III', 'Ni_I', 'Ni_II', 'Ni_III', 'Cu_I', 'Cu_II', 'Cu_III', 'Zn_I', 'Zn_II', 'Zn_III']

with open("all_emission.sh", "w") as file:
    file.write("#!/bin/bash \n")
    file.write('\n')
    file.write("[ ! -d './emission_absorption/ions' ] && mkdir emission_absorption/ions \n")
    file.write('\n')
    for i in ions:
        file.write("ion='" + i + "'\n")
        file.write('for i in {2..100}; do \n')
        file.write('  echo $ion $i \n')
        file.write("  plotartisspectrum -t $i --emissionabsorption --averagespecpolres --use_lastemissiontype -fixedionlist '" + i.replace('_', ' ') + "' --write_data \n")
        file.write("  rm plotspecemission_*pdf \n")
        file.write("  target=${ion/' '/'_'} \n")
        file.write("  [ ! -d './emission_absorption/ions/'$target ] && mkdir emission_absorption/ions/$target \n")
        file.write("  mv plotspecemission_*txt emission_absorption/ions/$target \n")
        file.write("done \n")
        file.write("cp pad_zeros.sh emission_absorption/ions/$target \n")
        file.write("cd emission_absorption/ions/$target \n")
        file.write("bash pad_zeros.sh \n")
        file.write("cd ../../../ \n")
        file.write("\n")
