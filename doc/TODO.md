# TODO

## RoadMap

- Refactoring:
  - move to modern tinyxml2 [ DONE ]
  - remove CUDA [ DONE ]
  - algorhythm simplification and optimization [ IN PROGRESS ]
- paralellize (openMP) [ DONE ]
- add openCL support
- migrate image processing from matlab to python+scipy/numpy/etc.

## Technologies to research and (maybe) implement

- C++ actor framework (CLF)
- C++11 (as an alternative to CLF)
- openCL (as part of CLF?)
- openMP (for cluster execution) (
- [PVM](https://en.wikipedia.org/wiki/Parallel_Virtual_Machine) (for cluster execution)
- FFT (with openMP, openCL)
- move to CMake (as universal platform)
- Boost (as multiplatform layer)
  - xml https://akrzemi1.wordpress.com/2011/07/13/parsing-xml-with-boost/ (as replacement to tinyXML2)
  - Boost.Thread (as an alternative to CLF)
  - Boost.Filesystem: files and directories autocreation
