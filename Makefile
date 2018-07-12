MKDIR                 = mkdir
ROOTDIR               = ./
SRCDIR                = $(ROOTDIR)src
OBJDIR                = $(ROOTDIR)target/obj
TARGETDIR             ?= $(ROOTDIR)target
INCLUDEDIR            = $(ROOTDIR)include
SUBDIRS               =
TESTSUBDIRS           = test/unit
EXES                  = pdp3
TESTDIR               = testdir
### Common settings
CEXTRA                =
CXXEXTRA              =
DEFINES               =
INCLUDE_PATH          = -I$(INCLUDEDIR)
LIBRARY_PATH          =
LIBRARIES             =
DOXYGEN               = doxygen
DOXYGEN_CONFIGS       = doc/app.conf doc/vis.conf
DOXYGEN_FORMATS       ?= latex html rtf
DOXYGEN_DIRS          = doc/app doc/vis
RELEASE               = CURRENT

### pdp3 sources and settings
pdp3_MODULE           = pdp3
pdp3_C_SRCS           =

pdp3_CXX_SRCS := $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(SRCDIR)/math/*.cpp)

pdp3_LDFLAGS          =
pdp3_ARFLAGS          =
pdp3_LIBRARY_PATH     =
pdp3_LIBRARIES        =
pdp3_RESULT_DIR       = pdp3_result

### tinyxml2 sources and settings
tinyxml2_SUBDIR := $(ROOTDIR)lib/tinyxml2
INCLUDE_PATH += -I$(tinyxml2_SUBDIR)
SUBDIRS += $(tinyxml2_SUBDIR)
pdp3_LIBRARIES += tinyxml2
LDFLAGS += -L$(tinyxml2_SUBDIR)

### tinyvec3d sources and settings
tinyvec3d_SUBDIR := $(ROOTDIR)lib/tinyvec3d
INCLUDE_PATH += -I$(tinyvec3d_SUBDIR)
SUBDIRS += $(tinyvec3d_SUBDIR)
pdp3_LIBRARIES += tinyvec3d
LDFLAGS += -L$(tinyvec3d_SUBDIR)

pdp3_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(pdp3_CXX_SRCS))

LDFLAGS := $(LDFLAGS) $(pdp3_LDFLAGS)

### Global source lists
C_SRCS                = $(pdp3_C_SRCS)
CXX_SRCS              = $(pdp3_CXX_SRCS)

### Tools
CC ?= gcc
CXX ?= g++

RC = wrc
AR = ar

CFLAGS_DEFAULT = -m64 -mcmodel=medium -std=c++1z
CFLAGS_SPEEDUP = -ffast-math -O3 -funroll-loops -Wno-write-strings -freciprocal-math -fcx-fortran-rules
CFLAGS_NO_SPEEDUP = -O2
CFLAGS_OPENMP ?= -fopenmp
CFLAGS_NO_OPENMP ?= -Wno-unknown-pragmas
CFLAGS_DEBUG ?= -O0 -Wall -Wextra -ggdb3 -fvar-tracking -ggnu-pubnames -pedantic -time -ftree-vectorizer-verbose=7 -DDEBUG -DPERF_DEBUG

CFLAGS += $(CFLAGS_DEFAULT)

## set debug options
ifeq ($(DEBUG), yes)
CFLAGS += $(CFLAGS_DEBUG)

## set speedup options
else ifeq ($(SPEEDUP), yes)
CFLAGS += $(CFLAGS_SPEEDUP)
else
CFLAGS += $(CFLAGS_NO_SPEEDUP)
endif

## set single thread options
ifeq ($(SINGLETHREAD), yes)
CFLAGS += $(CFLAGS_NO_OPENMP)
else
CFLAGS += $(CFLAGS_OPENMP)
endif

CXXFLAGS = ${CFLAGS} #  -fopenmp

### Generic targets
all: prepare $(SUBDIRS) $(LIBS) $(EXES)

### Build rules
.PHONY: all clean dummy check-syntax prepare doxygen

$(SUBDIRS): dummy
	cd $@ && $(MAKE)

$(DOXYGEN_FORMATS): doxygen
	@for i in $(DOXYGEN_DIRS); do test -f $(ROOTDIR)$$i/$@/Makefile && (cd $(ROOTDIR)$$i/$@ && $(MAKE)) || true; done

# Implicit rules
.SUFFIXES: .cpp .cxx
DEFINCL = $(INCLUDE_PATH) $(DEFINES) $(OPTIONS)

# obj/%.c.o:
# 	$(CC) -c $(CFLAGS) $(CEXTRA) $(DEFINCL) -o $@ $<

# obj/%.cpp.o:
# 	$(CXX) -c $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) -o $@ $<

# obj/%.cxx.o:
# 	$(CXX) -c $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) -o $@ $<

# Rules for cleaning
CLEAN_FILES  = y.tab.c y.tab.h lex.yy.c core *.orig *.rej \
               \\\#*\\\# *~ *% .\\\#*

BUILD_DIRS = pdp3_result pdp3_result/Dump $(OBJDIR) $(OBJDIR)/math
space :=
space +=
comma :=,

clean: $(SUBDIRS:%=%/__clean__) $(EXTRASUBDIRS:%=%/__clean__) $(TESTSUBDIRS:%=%/__clean__) $(DOXYGEN_FORMATS:%=%/__clean__)
	$(RM) $(pdp3_OBJS) $(CLEAN_FILES)
	$(RM) $(LIBS) $(EXES) $(EXES:%=%.so)
	$(RM) -r $(BUILD_DIRS) $(ROOTDIR)/tools/__pycache__/
	$(RM) -r $(DOXYGEN_DIRS)

$(SUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(MAKE) clean

$(DOXYGEN_FORMATS:%=%/__clean__): dummy
	-cd `dirname doc/$@` && test -f Makefile && $(MAKE) clean || return 0


$(EXTRASUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(RM) $(CLEAN_FILES)


$(TESTSUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(MAKE) clean

$(pdp3_MODULE): $(pdp3_OBJS)
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) $(LDFLAGS) -o $@ $(pdp3_OBJS) $(pdp3_LIBRARY_PATH) $(pdp3_LIBRARIES:%=-l%)

mrproper: clean
	$(RM) -r $(TESTDIR) *.avi $(TARGETDIR) $(RELEASE) $(RELEASE).zip

run: bootstrap
	./pdp3

$(TESTSUBDIRS:%=%/__test__): dummy
	cd `dirname $@` && $(MAKE) test

test-unit: $(TESTSUBDIRS:%=%/__test__)

test: mrproper all
	TESTDIR=$(TESTDIR) /bin/bash ./test/functional/test.sh

test-ext: mrproper all
	TESTDIR=$(TESTDIR) /bin/bash ./test/functional/test.sh extended

check-syntax:
	$(CXX) $(LIBRARY_PATH) $(INCLUDE_PATH) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

prepare:
	$(MKDIR) -p $(BUILD_DIRS)

doxygen:
	for i in $(shell ls doc/*.svg); do convert $$i $$(dirname $$i)/$$(basename $$i .svg).png; done
	for i in $(DOXYGEN_CONFIGS); do $(DOXYGEN) $$i; done

doc: doxygen $(DOXYGEN_FORMATS)

dist: all doc
	$(MKDIR) -p $(RELEASE)/$(pdp3_RESULT_DIR)
	cp pdp3 $(RELEASE)
	sed "s/<path_to_result>.*<\/path_to_result>/<path_to_result>${pdp3_RESULT_DIR}\/<\/path_to_result>/g" parameters.xml > $(RELEASE)/parameters.xml
	cp -r CHANGELOG.md $(RELEASE)
	cp -r tools $(RELEASE)/tools
	mkdir -p $(RELEASE)/doc
	cp doc/app/latex/refman.pdf $(RELEASE)/doc/pdp3.pdf
	cp doc/vis/latex/refman.pdf $(RELEASE)/doc/visualization.pdf
	zip -r `basename $(RELEASE)`.zip $(RELEASE)
