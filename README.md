# pdp3
Renew PDP3 project

[![Travis Build Status](https://api.travis-ci.org/cosmonaut-ok/pdp3.svg?branch=master)](https://travis-ci.org/cosmonaut-ok/pdp3)

RoadMap:

- Refactoring:
  - move to modern tinyxml2 [ DONE ]
  - remove CUDA
  - algorhythm simplification and optimization
- paralellize (CAF?)
- add openCL support
- wrappers for octave (CLI) to build images

https://github.com/leethomason/tinyxml2

## HOWTO

### Linux

1. clone project
```bash
user@host$ git clone https://github.com/cosmonaut-ok/pdp3.git
user@host$ cd pdp3
user@host$ git submodule update --init # require to enable external libraries  
```

2. compile

```bash
# change your current directory to project's root directory
user@host$ cd /path/to/pdp3/root/directory
user@host$ make
```
You need just file `pdp3` and `parameters.xml`. You can copy this files to somewhere, edit `parameters.xml` and run pdp3

3. test

```bash
user@host$ make test
```

## Hacking

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
