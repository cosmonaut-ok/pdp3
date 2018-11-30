CHANGELOG
---------

CURRENT
=========

- Add enable-legacy option

- Implement frame processing in tools

- Add classes to process probes and particles

- Implement data structure for reading and creating probes (writers)

- Add classes for new data format

- Add directoryExists and makeDirectory functions

release-20181107
=======

- Add stft jupyter notebook

- add frame-step option to skip some frames in tools/movie_3_E_RHObeam_r_z.py

- Add weiting schemes to documentation

- Add imagemagick and doxygen as package dependency

- Add cache to h5_reader

- Add plain to hdf5 converter

- Add compression support to tools

- Add batch mode to nbrun.sh

- Improve nbrun.sh to accept notebook variable names as arguments

- Update test tools

- Decrease number of thread synchronization points, optimize openmp parallelization (speedup)

- Add SSE sqrt to speedup calculations

- Improve autoconf options passing

- Add pdp3 general diagram

- Add caching to h5_reader

- Add plain to HDF5 converter

- Add optional compression to HDF5 format output

- Update tools/nbrun.sh to accept redefinition of any variables, defined in jupyter notebook

- Decrease number of thread synchronization points

- Improve speedup with some SSE support for intel/amd CPUs

- Improve multibunch support

- Migrate to autotools from simple Makefile

- Improve multibunch beam support

- Add profiler support

- Improve messages

- Fix bug with initial rotation coordinate

- Improve parameters passing to application in autoconf/autoheader

release-20180919
=======

- Add nbrun to launch jupyter/ipython notebooks from CLI

- Add full HDF5 support to pdp3 and tools

- Replace tools with completely new code

- Add tools for FFT transform plots

- Add caching to tools

- Add timestamps to 3E\rho_beam view/image/movie tools graphs

- Drop clang support

- Improve quick calculator

- Add tool to build 2.5D color map plot for E_r and E_z in r/t pane

- Fix bug with E_r and E_z labels

- Add with-grid option to set ticks grid

- Improve calculation steps algorithm

- Update tinyxml2 version

- Fix incorrect singlethread behavior

- Add tools to create plots for E_z(t) in fixed point (r, z)"

- tools for 2.5D plots for E_z(z, t) and E_r(z, t)

- Add compiler flag OPENMP_DYNAMIC_THREADS for omp dynamic threading

- Add automatic color limits estimation in tools for plot building

- Minor bugfixes

- Minor improvements

- Improve documentation

release-20180712
=======

- add possibility to inject series of bunches in one particle beam

- naming convenctions refactoring

- bugfixes in python tools

- improve python tools interface

- split python tools to images_3_component.py, movie_3_component.py and view_3_component.py to simplify CLI interface

- implement Parameters class to simplify parameters.xml reading

- optimize Makefile

- add performance debugging options (temporary solution. Profiler-friendly options required)

- add openmp threading with write locks to weighting code (experimental)

- update default CFLAGS set

- add more data dump parameters

- fix fields, positions and velocities dump

- Performance optimization

- move old coordinates definition to Particles class

- implement triple vectors as plain 3-component arrays as tinyvec3d library

- fix particles spatial distribution calculation

- fix particle relativistic calculations in pusher

- improve quick_parameters calculator

- rename directory python to tools

- add quick analytic calculator from parameters.xml

- fix particles distribution loading

- fix cylindrical coordinates pusher

- update python notebook for wake wave calculation

- move math functions to separate namespace

- move random_reverse function to lib.cpp as common for all code

- update defaults of parameters.xml to more realistic

- Add jupyter notebook for fast quick wavelength calculation

- Add debug option to parameters.xml

- remove half_step_coord from bunch as duplicated in parent class

- back to old step_v algorithm

- add some doxy-documentation

- add make target improve target

- minor fixes

release-20180314
================

- add documentation generation

- add make target 'dist' to automatically prepare pdp3 run environment

- improve velocity calculaction algorhythm

- add image set generation python script (#24)

- add data files set range feature for building video

- fix early gcc/clang build error, update tinyxml

- add automatic documentation generation with doxygen

- move common used functions to lib

- improve performance

- Minor build fix

- remove matlab files as unused

- update naming

- remove unneded calculation cycle operations

release-20180223
================

- Update data file names to human readable

- fix MSVS project paths

- change project structure

- update coding style guidelines, update naming, according to new coding standards

- update debug options in makefile

- fix msvs compatibility

- refactoring

- remove unused code

- Add unit testing

release-20180212
================

- Fix Makefile

- fix typo

- update tests naming

- remove most of matlab functions as outdated (replaced by python version)

- update messages

- updated time calculation

- minor bugfixes and updates

- Make python animation builder view mode more interactively, changed defaults

release-20180129
================

- fix MSVS

- Refactoring

- add experimental python/matplotlib support as alternative to matlab plot generation

- update travis. Test on gcc and clang

- add option to set configfile path

- add optimization options

- bugfixing

- add experimental pgi compiler support

- fix requirements

- add experimental intel c++ compiler support

- added clang support

release-20180122
================

- Add OPENMP support

- update documentation/README

- fix get_gamma with exception

- update Load_init_param

- Workaround for very fast speed behavior

- Fix gamma calculation

- Compiler and output optimisations

- fix segfault when calculate gamma

- improve output format

- bugfixing

- updated notification messages

- change defaults

- improve matlab plot generation code

- add function rho_movie_create_light3_from_parameters to get data from parameters.xml file, instead of set it as arguments

- improve parameters.xml format

- update parameters. unhardcode saving parameters

- add octave support to rho_movie_create_light3

release-20171029
================

- make few matlab function arguments optional

- move documentation to doc dir

- update matlab movie/graph generation script

- added basic script to matlab automation

- add script to install intel CPU/GPU openCL runtime

- optimize default values

- add extended tests (#11)

- update documentation/README/TODO

- Updated to compile with MSVS2013

- Update xml loading

- Fix loading params

- update tinyxml to latest version

- remove unused CUDA code

- add TODO

- archive unused matlab

- refactoring

- add experimental octave-based plot generation

- add automatic build and testing with travisCI

- add matlab files

- add functional tests (end-to-end)

- change dos to unix newline

- add gcc and `make` support

- remove unneded old files

base
====

- Fix include naming

- Improve comments

- Add x64 platform
