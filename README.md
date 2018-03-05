# PDP3

Renew PDP3 project

[![Travis Build Status](https://api.travis-ci.org/cosmonaut-ok/pdp3.svg?branch=master)](https://travis-ci.org/cosmonaut-ok/pdp3)

## System And Software Requirements

- Common:
  - C++ compiler: gcc 4.9+ or LLVM/clang 3.9+ or Intel C++ compiler 18.0+ (experimental, tests fails, but still working) or PGI compiler (experimental) or MSVS 2013 (single thread only)
  - OpenMP spec. version 3.0+ (see compiler requirements. As usual, openMP is a part of standard compiler libraries)
  - git (to get sources)
- Linux:
  - 'make' util
  - libgomp or libgomp (for clang or gcc)
  - python+matplotlib+numpy (or anaconda - python scientific environment)
- Windows:
  - MS Visual Studio 2017 or Cygwin with 'make' util (for gcc/clang/icc)
  - anaconda (python scientific environment)

## Terms and Legend

* <REQUIRED_VALUE> - required CLI value (ex. application's parameter, required to launch with)
* [OPTIONAL_VALUE] - optionsl CLI value (ex. application's optional launch parameter)
* VAR=val - set unix shell environment variable
* user@host$ - shell terminal prompt for user <user>
* root@host# - shell terminal prompt for superuser (root). You can reach it by command `su`, or `sudo -i`, or `sudo <command with all arguments>` from your user
* `# phrase` - comment in code (in shell-like syntax)

## HOWTO (Linux)

### 0. Install required software (debian/ubuntu example)

``` shell
root@host# apt-get install build-essential git doxygen texlive-latex-base
```
NOTE: if you are going to use LLVM/clang, you should install different packages

``` shell
root@host# apt-get install clang-<your faforite version> make git doxygen texlive-latex-base libomp5
```

### 1. **CLONE PROJECT**

```shell
user@host$ git clone https://github.com/cosmonaut-ok/pdp3.git
user@host$ cd pdp3
user@host$ git submodule update --init # require to enable external libraries
```

### 2. **COMPILE**

* **GCC**

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

* **Intel C compiler**

```shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make CXX=/path/to/binary/icc
```

* **PGI (Nvidia) C compiler**

```shell
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make CXX=/path/to/binary/pgc++ CFLAGS="-mp [other pgi-specific compile flags]"
```

#### Built-in compile flags

**Usage:**

```shell
user@host$ make [action] FLAG_1=value1 FLAG_2=value2
# or
user@host$ FLAG=value make [action]
```

**List of compile flags, related to PDP3 project:**

- `CXX=/foo/bar++` - Use custom c++ compiler (see supported c++ compilers list)
- `DEBUG=yes/no`- Compile binary with debug symbols, prepared to use with GDB
- `SPEEDUP=yes/no` - Increase speed up to 30%, by using unsafe math operations (gcc and clang only!). WARNING! it decreases calculations accuracy and can cause incorrect program working
- `SINGLETHREAD=yes/no` - Compile binary without multithreading support. Disables all parallelization features
- `CFLAGS="foo bar"` - list of custom CFLAGS (and CXXFLAGS), used by compiler


3. **TEST (optional)**

```bash
user@host$ make test # or test-ext for extended testing (require more time)
```

4. **RUN**

You need just file `pdp3` and `parameters.xml`. You can copy this files to somewhere, edit `parameters.xml` and run pdp3

NOTE: pgi: OMP_NUM_THREADS=N # !!!

```bash
user@host$ mkdir pdp_result # or where you defined in configfile
user@host$ ./pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # used parameters.xml from current directory, if calling without any options
## or (much better)
user@host$ nice -20 ./pdp3 [ -h ] | [ -f /path/to/parameters.xml ] # give some power to pdp3!
```

NOTE: it can take several days or weeks, and several hundreds gigabytes of diskspace (yep, it is science, my deer friend).

5. **VISUALIZATION**

After your application finished modeling, you can build some visual model from generated data. Use Mathlab for it (I know, that it is not good and matlab is not free. I plan to migrate visualization code to python+matplotlib+scipy. I will do it. I promise :) ).

``` shell
user@host$ cd /path/to/pdp3/matlab
user@host$ matlab -nodesktop -nosplash # f*ck that GUI sh*it!
...
>> rho_movie_create_light3_from_parameters('/path/to/parameters.xml') # it can take several hours
>> exit() # to exit after visualization finished
```

After visualization finished, please, look into directory with `parameters.xml` and find vieo file `field_movie.avi` there. Also, you can find other matlab scripts in `./matlab` pdp3 subdir.


## Hacking

### Coding Style

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

### Working with GIT

### Before work

1. Make sure, that you have clean repository and there are no unstaged untracked and uncommited files.

```bash
user@host$ git status
```
  1) If you have such files, you should remove them or move out of your working directory (with pdp3), than reset repository
  ```bash
  user@host$ git reset --hard
  ```
2. Sync with upstream
```bash
user@host$ git checkout master ## switch to master branch
user@host$ git pull origin master ## get latest stable code
```

### Begining work

1. Checkout to new branch
```bash
user@host$ git checkout -b your-branch-name # use only `-` as delimiters for better compatability
```

2. Make changes, write new code, add new features, experiment etc.

3. Commit your changes (do it as often, as save word document ;) ) and push to upstream
```bash
user@host$ git commit -a -m "your comment"
user@host$ git push origin your-branch-name
```

### Ending work

When you end some logical step of your working, you should merge your changes to master branch to stable code. You can do it, using "Pull request" in github

1. Go to https://github.com/cosmonaut-ok/pdp3/pulls

2. Press "New pull request" and choose "master" as base branch and "your-branch-name" as compare. Scroll to view changes against master branch

3. When all ok, press "Create pull request" and confirm.

4. Wait for travis testing ends and (if they passed), merge your changes to master branch

### Debug

Linux:
1. build project with DEBUG option

```bash
user@host$ DEBUG=yes make
```
or
```bash
user@host$ make DEBUG=yes
```

2. run with gdb

```bash
user@host$ gdb ./pdp3
(gdb) run ## or perform some modifications first than run, e.g. set breakpoints
```
