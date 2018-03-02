CURRENT
=======
- Alexander Vynnyk - Add image set generation python script (#24)

- Alexander Vynnyk - Add data files set range feature for building video

- Alexander Vynnyk - Merge branch 'master' of github.com:cosmonaut-ok/pdp3

- Alexander Vynnyk - Fix early gcc/clang build error, update tinyxml

- Alexander Vynnyk - documentation updates

- Alexander Vynnyk - Use README.md as main page

- Alexander Vynnyk - Initial for adding doxygen documentation (#23)

- Alexander Vynnyk - Merge branch 'master' of github.com:cosmonaut-ok/pdp3

- Alexander Vynnyk - Add data files set range feature for building video

- Alexander Vynnyk - Fix early gcc/clang build error, update tinyxml

- Alexander Vynnyk - documentation updates

- Alexander Vynnyk - Use README.md as main page

- Alexander Vynnyk - Initial for adding doxygen documentation (#23)

- Alexander Vynnyk - move common used functions to lib, optimize performance

- Alexander Vynnyk - Remove unneded small multiplications (#22)

- Alexander Vynnyk - Optimize main calculation cycle

- Alexander Vynnyk - Minor build fix

- Alexander Vynnyk - remove matlab files as unused

- Alexander Vynnyk - Remove unneded weighting

- Alexander Vynnyk - update naming beam->bunch

- Alexander Vynnyk - remove unneded calculation cycle operations, update naming


release-20180223
================
- Alexander Vynnyk - update naming

- Alexander Vynnyk - Update data file names to human readable

- Alexander Vynnyk - Fix python build script

- Alexander Vynnyk - fix error output

- Alexander Vynnyk - fix asterisk

- Alexander Vynnyk - clear MSVS project

- Alexander Vynnyk - Fix MSVS project paths

- Alexander Vynnyk - Unify field classes (#21)

    * Add class field as parent of eField hField
- Alexander Vynnyk - Remove particlesStruct as unused

- Alexander Vynnyk - Changing project structure

- Alexander Vynnyk - update fourier and hField according to coding standards

    * update loadInitParam according to coding standards
    * update Particles according to coding standards
    * update particlesList according to coding standards
    * update particlesStruct according to coding standards
    * update particlesStruct according to coding standards
    * update pdp3Time according to coding standards
    * update Poisson* according to coding standards
    * update Triple* according to coding standards
    * final update according to coding standards
    * update MSVS project according to coding standard
    * Update parameters.xml format
    * update methods naming according to code conventions
    * update tests according to code conventions
    * Syntax fixes
- Alexander Vynnyk - remove field as unused

- Alexander Vynnyk - Update coding style guidelines, update naming, according to new coding standards for geometry

    * Update boundaryMaxwellConditions according to coding standards
    * update bunch according to coding standards
    * update chargeDensity according to coding standards
    * update Constant according to coding standards
    * update eField and hField according to coding standards
- Alexander Vynnyk - Update debug options in Makefile

- Alexander Vynnyk - fix syntax and remove some MSVS old garbage

- Alexander Vynnyk - remove some garbage

- Alexander Vynnyk - Fix MSVS compatibility

- Alexander Vynnyk - Add initial tests for Particles

- Alexander Vynnyk - Remove Particles->continuity_equation as unused

- Alexander Vynnyk - Moved push_back inside of Particles to outside

- Alexander Vynnyk - Removed Boundary_Maxwell_conditions::radiation_source and Boundary_Maxwell_conditions::probe_mode_exitation as unused

- Alexander Vynnyk - removed boundary_conditions, boundary_dirichlet and boundary_neumann as unused

- Alexander Vynnyk - remove dummy unit tests

- Alexander Vynnyk - Remove pml class 2 (#20)

    * Add unit test for Geometry
    * minor E_field updates, refactoring
- Alexander Vynnyk - remove unused variant of E_field->calc_field method

- Alexander Vynnyk - Add unit testing (#18)

    * Add unit testing abilities
    * Add unit tests for Triple, Constant, PML, Geometry (basic)
    * Add dummy unit test files for all classes
    * Improve Geometry constructor

release-20180212

================
- Alexander Vynnyk - Fix Makefile

- Alexander Vynnyk - fix typo

- Alexander Vynnyk - update title name

- Alexander Vynnyk - improve average step execution time rounding

- Alexander Vynnyk - Include end_time to calculations

- Alexander Vynnyk - update tests naming

- Alexander Vynnyk - Remove most of matlab functions as outdated (replaced by python version)

- Alexander Vynnyk - Update messages

- Alexander Vynnyk - Updated time calculation, minor fixes and updates

- Alexander Vynnyk - Make python animation builder view mode more interactively, changed default


release-20180129
================
- Alexander Vynnyk - fix MSVS

- Alexander Vynnyk - Refactoring 7 (#17)
    * remove Bunch::bunch_inject_calc_E as unused
    * remove H_field::magnetostatic_equation as unused
    * remove wrapper.cpp and wrapper.h as unused
    * remove some opencl bindings as unused
    * remove some comments as unused or not informative
    * remove Fourier::fastcosinetransform_old as unused
    * add more parallelization
    * rename Particles::velocity_distribution_v2 to Particles::velocity_distribution and remove old version Particles::velocity_distribution as unused
    * remove commendted out Particles::get_cell_numbers_jr
- Alexander Vynnyk - Update Makefile

- Alexander Vynnyk - Merge branch 'master' of github.com:cosmonaut-ok/pdp3

- Alexander Vynnyk - Fix makefile

- Alexander Vynnyk - Matlab2python (#16)

- Alexander Vynnyk - Add type conversion
    * Fix travis
- Alexander Vynnyk - update test.sh work better with travis

- Alexander Vynnyk - update travis

- Alexander Vynnyk - reducing compiler variants number

- Alexander Vynnyk - update travis

- Alexander Vynnyk - remove gcc5

- Alexander Vynnyk - travis: fix openmp library for clang

- Alexander Vynnyk - update travis. Test on gcc and clang

- Alexander Vynnyk - update README, migrate matlab files to unix EOL

- Alexander Vynnyk - add option set configfile path

- Alexander Vynnyk - add optimization options

- Alexander Vynnyk - Fix typo

- Alexander Vynnyk - Add experimental PGI compiler support

- Alexander Vynnyk - Fix requirements

- Alexander Vynnyk - Add experimental intel C++ compiler support

- Alexander Vynnyk - Added clang support

release-20180122
================
- Alexander Vynnyk - Add OPENMP support (#15)

- Alexander Vynnyk - Fix labels, update README

- Alexander Vynnyk - Add params to save data, update notifications

- Alexander Vynnyk - fix get_gamma with exception

- Alexander Vynnyk - update Load_init_param

- Alexander Vynnyk - Workaround for very fast speed behavior

- Alexander Vynnyk - Fix style

- Alexander Vynnyk - Fix include

- Alexander Vynnyk - Fix gamma calculation

- Alexander Vynnyk - Minor fix

- Alexander Vynnyk - Aviod several exceptions

- Alexander Vynnyk - update build options

- Alexander Vynnyk - Minor bugfix

- Alexander Vynnyk - Compiler and output optimisations

- Alexander Vynnyk - Merge branch 'master' of github.com:cosmonaut-ok/pdp3

- Alexander Vynnyk - updated Makefile

- Alexander Vynnyk - avoid segfault when calculate gamma

- Alexander Vynnyk - Avoid unneded constant definition

- Alexander Vynnyk - Avoid some compilator warnings

- Alexander Vynnyk - Update app output format

- Alexander Vynnyk - update Makefile

- Alexander Vynnyk - Fix cell calculation, Makefile

- Alexander Vynnyk - Updated notification messages

- Alexander Vynnyk - change defaults

- Alexander Vynnyk - Update colorbar ticks

- Alexander Vynnyk - Convert E_bunch to electrons density

- Alexander Vynnyk - fix property name

- Alexander Vynnyk - Fix correct exiting from rho_movie_create_light3_from_parameters

- Alexander Vynnyk - Refactoring of rho_movie_create_light3_from_parameters

- Alexander Vynnyk - Improve matlab rho_movie_create_light3_from_parameters script

- Alexander Vynnyk - Add function rho_movie_create_light3_xml to get data from parameters.xml file, instead of set it as arguments

- Alexander Vynnyk - Fix some bugs, add rho_movie_create_light3_xml which reads data from properties.xml

- Alexander Vynnyk - beautify parameters.xml

- Alexander Vynnyk - Update parameters. unhardcode saving parameters

- Alexander Vynnyk - Add restricted octave support to rho_movie_create_light3 (#14)
	* Updated add octave support to rho_movie_create_light3

release-20171029
=======
- Alexander Vynnyk - Minor fix

- Alexander Vynnyk - Fix PDP3 project

- Alexander Vynnyk - Fix matlab graph function

- Alexander Vynnyk - Make few arguments optional

- Alexander Vynnyk - Move documentation to doc dir

- Alexander Vynnyk - Refactoring 6 (#13)
	* update README
	* Remove unused code
- Alexander Vynnyk - Improve travis builds

- Alexander Vynnyk - Cleaning

- Alexander Vynnyk - Minor update

- Alexander Vynnyk - Update movie generation script

- Alexander Vynnyk - Added basic script to matlab automation

- Alexander Vynnyk - Update matlab function

- Alexander Vynnyk - Add TODO

- Alexander Vynnyk - Add script to install intel CPU/GPU openCL runtime

- Alexander Vynnyk - update README

- Alexander Vynnyk - Removed unused true_data archive

- Alexander Vynnyk - Refactoring 5 (#12)
	* optimize default values

- Alexander Vynnyk - Add extended tests (#11)

- Alexander Vynnyk - Update documentation

- Alexander Vynnyk - Refactoring 4 (#10)
	* Remove not used Beam inclusions
	* rename Time.cpp to pdp3_time.cpp

- Alexander Vynnyk - Merge pull request #9 from cosmonaut-ok/refactoring-3
	* Remove unused code

- Alexander Vynnyk - Remove unused code

- Alexander Vynnyk - minor changes

- Alexander Vynnyk - Minor fixes

- Alexander Vynnyk - Updated to compile with MSVS2013

- Alexander Vynnyk - Update VS project

- Alexander Vynnyk - Merge pull request #8 from cosmonaut-ok/refactoring-xml-parser
	* Refactoring xml parser

- Alexander Vynnyk - Small format fix

- Alexander Vynnyk - Update xml loading

- Alexander Vynnyk - Merge pull request #7 from cosmonaut-ok/fix-loading-params
	* Fix loading params

- Alexander Vynnyk - Added some instructions to README

- Alexander Vynnyk - Fix travis

- Alexander Vynnyk - Moved to modern tinyxml version

- Alexander Vynnyk - Merge pull request #6 from cosmonaut-ok/remove-unused-cuda
	* Remove unused CUDA code

- Alexander Vynnyk - Remove unused CUDA code

- Alexander Vynnyk - Add TODO

- Alexander Vynnyk - Merge pull request #5 from cosmonaut-ok/archive-unused-matlab
	* Archive unused matlab

- Alexander Vynnyk - Archive unused matlab

- Alexander Vynnyk - Merge pull request #4 from cosmonaut-ok/refactoring-2
	* refactoring continue

- Alexander Vynnyk - Refactored till Particles.cpp
	* Fix true data
	* Refactored
	* Removed another local constant
	* Update test_true data

- Alexander Vynnyk - refactoring continue

- Alexander Vynnyk - Add getframe emulation for octave

- Alexander Vynnyk - update travis

- Alexander Vynnyk - Add travis.yml

- Alexander Vynnyk - add matlab files

- Alexander Vynnyk - Clearing, add tests

- Alexander Vynnyk - Merge pull request #1 from cosmonaut-ok/moving-to-makefile
	* Moving to makefile

- Alexander Vynnyk - change dos to unix newline

- Alexander Vynnyk - Moving to gcc and Makefile

- Alexander Vynnyk - Remove unneded old files

base
====
- Alexander Vynnyk - Fix include naming

- Alexander Vynnyk - Add some explanations

- Alexander Vynnyk - Add x64 platform

