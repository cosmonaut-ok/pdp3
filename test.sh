#!/bin/sh

test -z $TESTDIR && TESTDIR=testingdir
test -z $TESTING_DATA_ARCHIVE && TESTING_DATA_ARCHIVE=test_true_data.tar.gz

mkdir -p ${TESTDIR}
cp pdp3 ${TESTDIR}

cat<<EOF>${TESTDIR}/parameters.xml
<?xml version="1.0" encoding="UTF-8"?>
<!-- edited with XMLSpy v2011 rel. 2 -->
<Initial_parameters>

  <geometry>
    <r_size>   0.25  </r_size>
    <z_size>  2.0  </z_size>
    <n_grid_r>  255  </n_grid_r>
    <n_grid_z> 2047</n_grid_z>
  </geometry>

  <PML>
    <comparative_l_1>0</comparative_l_1>
    <comparative_l_2> 0 </comparative_l_2>
    <comparative_l_3> 0 </comparative_l_3>
    <sigma_1> 0.00001 </sigma_1>
    <sigma_2> 0.07 </sigma_2>
  </PML>

  <Time>
    <start_time> 0 </start_time>
    <relaxation_time> 0 </relaxation_time>
    <current_time> 0 </current_time>
    <end_time>2.2e-11</end_time>
    <delta_t>1e-12</delta_t>
  </Time>

  <Number_of_part_kinds>  2 </Number_of_part_kinds>

  <Particle_Name>  Electrons </Particle_Name>
  <Electrons>
    <charge> -1 </charge>
    <mass> 1 </mass>
    <number> 1e6</number>
    <left_density> 2e14</left_density>
    <right_density>2.05e14  </right_density>
    <temperature> 1.0 </temperature>
  </Electrons>

  <Particle_Name> ions </Particle_Name>
  <ions>
    <charge> 1 </charge>
    <mass> 1836 </mass>
    <number> 1e6 </number>
    <left_density> 2e14 </left_density>
    <right_density> 2.05e14 </right_density>
    <temperature> 0.1 </temperature>
  </ions>


  <Inject_Particles>  Electrons_bunch </Inject_Particles>
  <Electrons_bunch>
    <charge> -1 </charge>
    <mass> 1 </mass>
    <number> 1e6</number>
    <duration>1e-8</duration>
    <Radius>0.02</Radius>
    <density> 1e13 </density>
    <initial_velocity> 3e7 </initial_velocity>

  </Electrons_bunch>


  <Boundary_Maxwell_conditions>
    <E_fi_upper> 0  </E_fi_upper>
    <E_fi_left> 0 </E_fi_left>
    <E_fi_right> 0 </E_fi_right>
  </Boundary_Maxwell_conditions>

  <Boundary_conditions type = "dirichlet">0</Boundary_conditions>


  <PathtoResult>pdp3_result/</PathtoResult>
  <PathtoSaveState>pdp3_result/Dump/</PathtoSaveState>

  <file_save_paramters>
    <time_step>100</time_step>
    <number>100</number>
  </file_save_paramters>
  <out_folder>pdp3_result/</out_folder>
  <e1>true</e1>
  <e2>true</e2>
  <e3>false</e3>
  <h1>true</h1>
  <h2>true</h2>
  <h3>false</h3>
  <rho_plasma>true</rho_plasma>
  <rho_beam>true</rho_beam>





</Initial_parameters>

EOF

cd ${TESTDIR}
mkdir -p pdp3_files pdp3_result/Dump
./pdp3

## compare with testing data
mkdir -p true_data
cd true_data
tar -xf ../../${TESTING_DATA_ARCHIVE}

success="true"
for i in `find . -type f`; do
    bn=`basename $i`
    test -z "$(diff ${bn} ../pdp3_result/${bn} 2>&1)" || success="false"
done

if [ "$success" == "false" ]; then
    echo
    echo "Test failed. Generated data is not match control data"
    exit 1
else
    echo
    echo "Test passed."
    exit 0
fi
