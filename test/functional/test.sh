#!/bin/bash

set -e

test -z $TESTDIR && TESTDIR=testingdir
test -z $TRUE_MD5 && TRUE_MD5=test/functional/true_md5sums
test -z $TRUE_EXT_MD5 && TRUE_EXT_MD5=test/functional/true_ext_md5sums

if [ -n "$1" ] && [ "$1" == "extended" ]; then
    END_TIME=1.05e-10
    TRUE_MD5=$TRUE_EXT_MD5
else
    END_TIME=2.2e-11
fi

mkdir -p ${TESTDIR}
cp pdp3 ${TESTDIR}

cat<<EOF>${TESTDIR}/parameters.xml
<?xml version="1.0" encoding="UTF-8"?>
<!-- edited with XMLSpy v2011 rel. 2 -->
<initial_parameters>

  <geometry>
    <r_size>   0.25  </r_size>
    <z_size>  2.0  </z_size>
    <n_grid_r>  255  </n_grid_r>
    <n_grid_z> 2047</n_grid_z>
  </geometry>

  <pml>
    <comparative_l_1>0</comparative_l_1>
    <comparative_l_2> 0 </comparative_l_2>
    <comparative_l_3> 0 </comparative_l_3>
    <sigma_1> 0.00001 </sigma_1>
    <sigma_2> 0.07 </sigma_2>
  </pml>

  <time>
    <start_time> 0 </start_time>
    <relaxation_time> 0 </relaxation_time>
    <current_time> 0 </current_time>
    <end_time>${END_TIME}</end_time>
    <delta_t>1e-12</delta_t>
  </time>

  <number_of_part_kinds>  2 </number_of_part_kinds>


    <particles>
    <particle_kind>
      <name>Electrons</name>
      <charge> -1 </charge>
      <mass> 1 </mass>
      <number> 1e6</number>
      <left_density> 2e14</left_density>
      <right_density>2.05e14  </right_density>
      <temperature> 1.0 </temperature>
    </particle_kind>

    <particle_kind>
      <name>Ions</name>
      <charge> 1 </charge>
      <mass> 1836 </mass>
      <number> 1e6 </number>
      <left_density> 2e14 </left_density>
      <right_density> 2.05e14 </right_density>
      <temperature> 0.1 </temperature>
    </particle_kind>
  </particles>


  <particles_bunch>
    <name> Electrons </name>
    <charge> -1 </charge>
    <mass> 1 </mass>
    <number> 1e6</number>
    <duration>1e-8</duration>
    <radius>0.02</radius>
    <density> 1e13 </density>
    <initial_velocity> 3e7 </initial_velocity>
  </particles_bunch>


  <boundary_maxwell_conditions>
    <e_fi_upper> 0  </e_fi_upper>
    <e_fi_left> 0 </e_fi_left>
    <e_fi_right> 0 </e_fi_right>
  </boundary_maxwell_conditions>

  <boundary_conditions type = "dirichlet">0</boundary_conditions>

  <file_save_parameters>
    <path_to_result>pdp3_result/</path_to_result>
    <path_to_save_state>pdp3_result/dump/</path_to_save_state>
    <data_dump_interval>5</data_dump_interval>
    <frames_per_file>100</frames_per_file>
    <system_state_dump_interval>1000</system_state_dump_interval>
    <dump_data>
      <e1>true</e1>
      <e2>false</e2>
      <e3>true</e3>
      <h1>false</h1>
      <h2>false</h2>
      <h3>false</h3>
      <rho_beam>true</rho_beam>
    </dump_data>
  </file_save_parameters>

</initial_parameters>

EOF

cd ${TESTDIR}
mkdir -p pdp3_result/Dump
time ./pdp3

cd pdp3_result

## compare with testing data
success="true"
for i in `find . -type f`; do
    bn=`basename $i`
    true_md5sum="$(grep ${bn} ../../${TRUE_MD5} | cut -d';' -f1)"
    actual_md5sum="$(md5sum ${bn} | cut -d' ' -f1)"
    test "${true_md5sum}" == "${actual_md5sum}" || success="false"
done

if [ "$success" == "false" ]; then
    echo
    echo "Test failed. Generated data is not match control data"
    exit 1
else
    echo
    echo "Test passed."
    echo
    exit 0
fi
