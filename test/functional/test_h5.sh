#!/bin/bash

set -e

test -z $TESTDIR && TESTDIR=testingdir_h5
test -z $TRUE_MD5 && TRUE_MD5=test/functional/true_md5sums
test -z $TRUE_EXT_MD5 && TRUE_EXT_MD5=test/functional/true_ext_md5sums

## Set some colors
RED='\e[31m'
GREEN='\e[32m'
NC='\e[0m' # No Color

if [ -n "$1" ] && [ "$1" == "extended" ]; then
    END_TIME=1.05e-10
    TRUE_MD5=$TRUE_EXT_MD5
else
    END_TIME=2.2e-11
fi

test -f ${TESTDIR} && rm -rf ${TESTDIR}
mkdir -p ${TESTDIR}
cp pdp3 ${TESTDIR}
cp -r tools ${TESTDIR}

cat<<EOF>${TESTDIR}/parameters.xml
<?xml version="1.0" encoding="UTF-8"?>
<!-- edited with XMLSpy v2011 rel. 2 -->
<initial_parameters>

  <debug>false</debug>
  <use_hdf5>true</use_hdf5>

  <geometry>
    <r_size>0.25</r_size>
    <z_size>2.0</z_size>
    <debug_n_grid_r>15</debug_n_grid_r>
    <debug_n_grid_z>127</debug_n_grid_z>
    <n_grid_r>255</n_grid_r>
    <n_grid_z>2047</n_grid_z>
  </geometry>

  <pml>
    <comparative_l_1>0</comparative_l_1>
    <comparative_l_2>0</comparative_l_2>
    <comparative_l_3>0</comparative_l_3>
    <sigma_1>1e-5</sigma_1>
    <sigma_2>7e-2</sigma_2>
  </pml>

  <time>
    <start_time> 0 </start_time>
    <relaxation_time> 0 </relaxation_time>
    <current_time> 0 </current_time>
    <end_time>${END_TIME}</end_time>
    <step_interval>1e-12</step_interval>
  </time>

  <particles>
    <particle_specie name="Electrons">
      <charge>-1</charge>
      <mass>1</mass>
      <number_macro>1e5</number_macro>
      <debug_number_macro>1e5</debug_number_macro>
      <left_density>2e16</left_density>
      <right_density>2.05e16</right_density>
      <temperature>1.0</temperature>
    </particle_specie>

    <particle_specie name="Ions">
      <charge>1</charge>
      <mass>1836</mass>
      <number_macro>1e5</number_macro>
      <debug_number_macro>1e5</debug_number_macro>
      <left_density>2e16</left_density>
      <right_density>2.05e16</right_density>
      <temperature>0.1</temperature>
    </particle_specie>
  </particles>

  <particle_beam name="Beam">
    <charge>-1</charge>
    <mass>1</mass>
    <initial_velocity>2e8</initial_velocity>
    <number_bunches>1</number_bunches> <!-- NOT IMPLEMENTED -->
    <bunches_distance>0.02</bunches_distance> <!-- NOT IMPLEMENTED -->

    <bunch_number_macro>1e6</bunch_number_macro> <!-- number of macroparticles in one bunch -->
    <debug_bunch_number_macro>1e4</debug_bunch_number_macro>
    <!-- <duration>1e-10</duration> -->

    <bunch_length>2</bunch_length>
    <bunch_radius>0.02</bunch_radius>
    <bunch_density>1e15</bunch_density>
  </particle_beam>

  <boundary_maxwell_conditions>
    <e_fi_upper>0</e_fi_upper>
    <e_fi_left>0</e_fi_left>
    <e_fi_right>0</e_fi_right>
  </boundary_maxwell_conditions>

  <boundary_conditions type = "dirichlet">0</boundary_conditions>

  <file_save_parameters>
    <result_path>pdp3_result/</result_path>
    <save_state_path>pdp3_result/dump/</save_state_path>
    <debug_data_dump_interval>5</debug_data_dump_interval>
    <data_dump_interval>5</data_dump_interval>
    <frames_per_file>1</frames_per_file>
    <system_state_dump_interval>2</system_state_dump_interval>
    <dump_data>
      <E_r>true</E_r>
      <E_phi>false</E_phi>
      <E_z>true</E_z>
      <H_r>false</H_r>
      <H_phi>false</H_phi>
      <H_z>false</H_z>
      <J_r>true</J_r>
      <J_phi>false</J_phi>
      <J_z>true</J_z>

      <position>false</position>

      <velocity>true</velocity>

      <rho_bunch>true</rho_bunch>
    </dump_data>
  </file_save_parameters>

  <plot>
    <video>
      <codec>mjpeg</codec>
      <fps>20</fps>
      <dpi>100</dpi>
      <bitrate>32000</bitrate>
    </video>
  </plot>

</initial_parameters>

EOF

cd ${TESTDIR}
mkdir -p pdp3_result/dump
time ./pdp3

## dump data from hdf5 database
echo -n "Dumping state data from database to plaintext for test..."
for i in `h5ls -r pdp3_result/data.h5/pdp3/result | grep Dataset | awk '{print $1}'`; do
    h5dump -d /pdp3/result/$i -y -r -o pdp3_result/$(basename $(dirname $i))$(basename $i)_raw pdp3_result/data.h5 > /dev/null
done
echo "done"

echo -n "Dumping data from database to plaintext for test..."
for i in `h5ls -r pdp3_result/data.h5/pdp3/dump | grep Dataset | awk '{print $1}'`; do
    h5dump -d /pdp3/dump/$i -y -r -o pdp3_result/$(basename $i)_raw pdp3_result/data.h5 > /dev/null
done
echo "done"

cd pdp3_result

echo -n "Converting RAW data to format for md5sums comparation..."
for i in `find . -name '*_raw'`; do
    echo $(cat ${i} | tr ',' ' ') | tr ' ' '\n' > $(basename $i _raw)n
    rm -f ${i}
done
echo "done"

## compare with testing data
success="true"
for i in `find . -type f -regex '^.*[0-9]n$'`; do
    bn=`basename $i`
    true_md5sum="$(egrep ${bn}$ ../../${TRUE_MD5} | cut -d';' -f1)"
    actual_md5sum="$(md5sum ${i} | cut -d' ' -f1)"
    if test "${true_md5sum}" == "${actual_md5sum}"; then
	printf "%-25b %-10b\n" "file $bn" "${GREEN}match${NC}"
    else
	printf "%-25b %-10b\n" "file $bn" "${RED}not match${NC}"
	success="false"
    fi
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