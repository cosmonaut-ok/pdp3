# PDP3

Renew PDP3 project

[![Travis Build Status](https://api.travis-ci.org/cosmonaut-ok/pdp3.svg?branch=master)](https://travis-ci.org/cosmonaut-ok/pdp3)

### System And Software Requirements

- Common:
  - C++ compiler: gcc 4.9+ or LLVM/clang 3.9+ or Intel C++ compiler 18.0+ (experimental, tests fails, but still working) or PGI compiler (experimental) or MSVS 2013 (single thread only)
  - OpenMP spec. version 3.0+ (see compiler requirements. As usual, openMP is a part of standard compiler libraries)
  - git (to get sources)
- Linux:
  - 'make' util
  - libomp5 or libgomp (for clang or gcc)
  - python+matplotlib+numpy (or anaconda - python scientific environment)
- Windows:
  - MS Visual Studio 2017 or Cygwin with 'make' util (for gcc/clang/icc)
  - anaconda (python scientific environment)

### Terms and Legend

* `<REQUIRED_VALUE>` - required CLI value (ex. application's parameter, required to launch with)
* `[OPTIONAL_VALUE]` - optionsl CLI value (ex. application's optional launch parameter)
* `VAR=val` - set unix shell environment variable
* `user@host$` - shell terminal prompt for user <user>
* `root@host#` - shell terminal prompt for superuser (root). It can be reached by commands `su`, or `sudo -i`, or `sudo <command with all arguments>` from your user
* `# phrase` - comment in code (in shell-like syntax)

### HOWTO (linux, debian/ubuntu example)

#### 0. INSTALL PREREQUIRED SOFTWARE

``` shell
root@host# apt-get install build-essential git
```
> NOTE: if you are going to use LLVM/clang, you should install different packages

``` shell
root@host# apt-get install clang-<your faforite version> make git libomp5
```

#### 1. **CLONE PROJECT**

```shell
user@host$ git clone https://github.com/cosmonaut-ok/pdp3.git
user@host$ cd pdp3
user@host$ git submodule update --init # require to enable external libraries
```

#### 2. **COMPILE**

* **GCC (recommended)**

```shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make [COMPILE FLAGS] # see below about COMPILE FLAGS
```

* **LLVM/Clang**

```shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make CXX=clang++-<your faforite version> CFLAGS_OPENMP=-fopenmp=libiomp5 [OTHER COMPILE FLAGS] # see below about COMPILE FLAGS
```

* **PGI (Nvidia) C compiler**

```shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make CXX=/path/to/binary/pgc++ CFLAGS="-mp [other pgi-specific compile flags]"
```

#### Compile flags

**Usage:**

```shell
user@host$ make [action] FLAG_1=value1 FLAG_2=value2
# or
user@host$ FLAG=value make [action]
```

**List of built-in compile flags, accepted by PDP3 project:**

- `DEBUG=yes/no`- Compile binary with debug symbols, prepared to use with GDB
- `SPEEDUP=yes/no` - Increase speed up to 30%, by using unsafe math operations (gcc and clang only!). WARNING! it decreases calculations accuracy and can cause incorrect program working
- `SINGLETHREAD=yes/no` - Compile binary without multithreading support. Disables all parallelization features
- `CXX=/foo/bar++` - Use custom c++ compiler (see supported c++ compilers list)
- `CXXFLAGS="foo bar"` - custom C++ compile flags, used by compiler (overrides all autogeneration flags)
- `CFLAGS="foo bar"` - same as CXXFLAGS
  - `CFLAGS_NO_OPENMP` - used to customize flags, when openmp disabled (used for CXXFLAGS autogeneration)
  - `CFLAGS_OPENMP` -  - used to customize flags, when openmp enabled (used for CXXFLAGS autogeneration)
  - `CFLAGS_DEBUG` - used to customize flags, used when debug enabled (used for CXXFLAGS autogeneration)
  - `CFLAGS_SPEEDUP` - used to customize flags, used when fast-math and other speedup options enabled (used for CXXFLAGS autogeneration)
- `LDFLAGS` -  custom C++ linker flags
- `TARGETDIR` - used with `make dist`, which prepares all, required to start modeling in separate directory

#### 3. **TEST (optional)**

**Functional (end-to-end) testing:**

```shell
user@host$ make test # or test-ext for extended testing (require more time)
```

**Unit testing:**

```shell
user@host$ make test-unit
```

#### 4. **INSTALLATION (optional)**

You can install already compiled pdp3 with it's configfile (aka properties.xml) and required subdirs to separate directory. Just run

``` shell
user@host$ make dist [TARGETDIR=/path/to/some/target/directory]
```

Than, you can go to this directory and run pdp3

#### 5. **RUN**

After compilation finished, you just need binary file `pdp3` and `parameters.xml`. You can copy this files to somewhere, edit `parameters.xml` and run pdp3

```shell
user@host$ mkdir pdp_result # or where you defined in parameters.xml. PDP3 does not use smth. like BOOST::filesystem to operate with directories
# NOTE: if you preformed `make dist`, you just need `cd /path/to/target/dir`, edit `parameters.xml` (optional) and run `./pdp3`
user@host$ /path/to/pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # used parameters.xml from current directory, if calling without any options
# or (much better), launch it as high-priority process
user@host$ nice -20 /path/to/pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # give some power to pdp3!
```
> NOTE: it can take several days or weeks, and several hundred gigabytes of diskspace (yep, it is science, my deer friend).

#### 6. **VISUALIZATION (generate images or animation)**

After your application finished modeling, you can build some visual model from generated data. Use python with matplotlib and numpy. [Anaconda](https://www.anaconda.com/download/#linux) as python distribution is recommended.

**...but... if you don't want to use anaconda...:**

``` shell
root@host# apt-get install python3 python3-matplotlib python3-numpy python3-scipy
```

**Pre-required software (for animation):**

``` shell
root@host# apt-get install ffmpeg
```

**Animation generation:**

``` shell
user@host$ /path/to/repository/with/pdp3/tools/movie_3_component.py /path/to/parameters.xml
# USE: /path/to/repository/with/pdp3/tools/movie_3_component.py -h to see list of all available options
```
> NOTE: you may use smth. like `python3 /path/to/repository/with/pdp3/tools/movie_3_component.py /path/to/parameters.xml`, if you don't use anaconda.

**Images generation:**

``` shell
user@host$ /path/to/repository/with/pdp3/tools/images_3_component.py /path/to/parameters.xml --data_set_range=1:2 --frame_range=3:4
# USE: /path/to/repository/with/pdp3/tools/images_3_component.py -h to see list of all available options
```

WAT is *data_set_range* and *frame_range* ? PDP3 saves every modeling frame (step) to file. Number of such frames in one file can be defined in in parameters.xml. So, you getting output as set of files with set of frames in each file. When you generate images, you can define range of files from which you going to generate images and range of frames in each file, from which that images will be generated.

#### 7. **DOCUMENTATION GENERATION** (optional)

**Pre-required software:**

``` shell
root@host# doxygen texlive-latex-base texlive-latex-extra imagemagick
```
**Generate:**

``` shell
user@host$ make doc
```
Find documentation in:

- PDF
  `doc/app/latex/refman.pdf` - application
  `doc/vis/latex/refman.pdf` - visualization
- HTML
  `doc/app/html/index.xhtml` - application
  `doc/vis/html/index.xhtml` - visualization

### Hacking

#### Coding Style

- Indentication: 2 spaces
- Braces: BSD style
- Naming:
  - classes - CamelCase
  - files - camelCase, same as class name, begining from lower case
  - functions - snake_case, lower case
- Spacing:
  - after/around operators
  - no space after unary operators
  - no space before the postfix increment & decrement unary operators
  - no space around the . and -> structure member operators
- Declaring pointer data or a function that returns a pointer type: preferred to use of * is adjacent to the data name or function name and not adjacent to the type name
- Comments: '//' at line begining
- File naming: lowercase

#### GIT

##### Before working

1. Make sure, that you have clean repository and there are no unstaged untracked and uncommited files.

``` shell
user@host$ git status
```
> NOTE: If you have such files, you commit, remove or move out of your working directory (with pdp3), than reset repository
> Also, you can just reset repository
``` shell
user@host$ git reset --hard
```

2. Sync with upstream

``` shell
user@host$ git fetch origin --tags
```

##### Working

1. Checkout to new branch

``` shell
user@host$ git checkout -b <your-branch-name> # use only `-` as word delimiters for better compatability
```

2. Make changes, write new code, add new features, experiment etc.

3. Commit your changes (do it as often, as you press `Ctrl+s` ;) ) and push to upstream

``` shell
user@host$ git add . # or git add your/selected/files
user@host$ git commit -m "<your comment>"
user@host$ git push origin <your-branch-name>
```

##### After working

When you finish some logical step of your work, you should merge your changes to master branch (as stable code). You can do it, using "Pull request" in github

1. Go to https://github.com/cosmonaut-ok/pdp3/pulls

2. Press "New pull request" and choose "master" as base branch and "your-branch-name" as compare. Scroll to view changes against master branch

3. When all is ok, press "Create pull request" and confirm.

4. Wait for travis testing ends and (if they passed), merge your changes to master branch (squash and merge).

#### Debug

**Linux:**

1. build project with DEBUG option
``` shell
user@host$ make DEBUG=yes
```

2. Set options `debug` to `true`,

3. Check out options
   - `particles->particles_kind->debug_number` (at least `1e4` is recommended)
   - `particles_bunch->debug_number` (at least `1e4` is recommended)
   - `geometry->debug_n_grid_r`
   - `geometry->debug_n_grid_z`
   - `file_save_parameters->debug_data_dump_interval`
in `parameters.xml` file before project run.

4. run with gdb

``` shell
user@host$ gdb ./pdp3
(gdb) run ## or perform some modifications first than run, e.g. set breakpoints
```

#### Tools

Quick analytic calculator of plasma (aka Langmur) frequency, wake wavelength, Debye length etc. from parameters.xml file
``` shell
user@host$ ./tools/quick_calculator.py <path/to/parameters.xml>
```