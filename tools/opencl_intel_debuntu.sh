#!/bin/sh
set e
# source: https://gist.github.com/rmcgibbo/6314452

OPENCL_RUNTIME="opencl_runtime_14.2_x64_4.5.0.8.tgz"
TARGET_URL="http://registrationcenter.intel.com/irc_nas/4181/"$OPENCL_RUNTIME


# download and extract opencl runtime
wget $TARGET_URL
tar xzvf $OPENCL_RUNTIME

# install rpm converter
sudo apt-get install -y \
    rpm \
    alien \
    libnuma1 \
    ocl-icd-libopencl1 \
    opencl-headers

# convert opencl rpms into dpkgs and install them
cd pset_opencl_runtime_14.1_x64_4.5.0.8/rpm

for f in *.rpm; do
    fakeroot alien --to-deb $f
done
for f in *.deb; do
    sudo dpkg -i $f
done

cd ~


# CREATE ICD SYMLINKS
# Install the so-called icd-file, which registers this OpenCL implementation,
# so that it's available in paralell to any other, this isn't done by default.
# The icd files live at /etc/OpenCL/vendors/*.icd and these files tell the ICD
# loader what OpenCL implementations (ICDs) are installed on your machine.
# There's one for each ICD.  Each file is a one-line text file containing the
# name of the dynamic library (aka shared object, aka ".so" file) containing
# the implementation. The single line may either be the full absolute path or
# just the file name, in which case the dynamic linker must be able to find
# that file--perhaps with the help of setting the LD_LIBRARY_PATH environment
# variable. The names of the .icd files themselves are arbitrary, but they must
# have a file extension of .icd. To install the icd file

sudo mkdir -p /etc/OpenCL/vendors
sudo ln -s /opt/intel/opencl-1.2-4.5.0.8/etc/intel64.icd /etc/OpenCL/vendors/intel64.icd


# CREATE LIBRARY SYMLINKS
# If this is the only OpenCL implementation on your machine, you should install a
# symlink to libOpenCL.so into /usr/lib, so that things can be linked up
# easily. If you already have the NVIDIA OpenCL platform (for your GPU) then
# this is not necessary -- installing the icd file into the registry is
# enough to tell the system about your new OpenCL platform.

sudo ln -s /opt/intel/opencl-1.2-4.5.0.8/lib64/libOpenCL.so /usr/lib/libOpenCL.so
sudo ldconfig
