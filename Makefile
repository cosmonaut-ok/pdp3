MKDIR                 = mkdir
ROOTDIR               = .
SRCDIR                = $(ROOTDIR)/src
OBJDIR                = $(ROOTDIR)/obj
INCLUDEDIR            = $(ROOTDIR)/include
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

### pdp3 sources and settings
pdp3_MODULE           = pdp3
pdp3_C_SRCS           =

pdp3_CXX_SRCS := $(wildcard $(SRCDIR)/*.cpp)

pdp3_LDFLAGS          =
pdp3_ARFLAGS          =
pdp3_LIBRARY_PATH     =
pdp3_LIBRARIES        =

### tinyxml2 sources and settings
tinyxml2_SUBDIR := $(ROOTDIR)/lib/tinyxml2
INCLUDE_PATH += -I$(tinyxml2_SUBDIR)
SUBDIRS += $(tinyxml2_SUBDIR)
pdp3_LIBRARIES += tinyxml2

LDFLAGS += -L$(tinyxml2_SUBDIR)

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

CFLAGS ?= -m64 -mcmodel=medium
CFLAGS_SPEEDUP = -ffast-math
CFLAGS_NO_OPENMP ?= -Wno-unknown-pragmas
CFLAGS_OPENMP ?= -fopenmp
CFLAGS_DEBUG ?= -O0 -Wall -Wextra -g -ggdb -fvar-tracking -ggnu-pubnames -pedantic

## set debug options
ifeq ($(DEBUG), yes)
CFLAGS += $(CFLAGS_DEBUG)
else
CFLAGS += -O2
endif

## set speedup options
ifeq ($(SPEEDUP), yes)
CFLAGS += $(CFLAGS_SPEEDUP)
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
.PHONY: all clean dummy check-syntax prepare

$(SUBDIRS): dummy
	@cd $@ && $(MAKE)

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

BUILD_DIRS = pdp3_files pdp3_result/Dump $(OBJDIR) pdp3_result

clean: $(SUBDIRS:%=%/__clean__) $(EXTRASUBDIRS:%=%/__clean__) $(TESTSUBDIRS:%=%/__clean__)
	$(RM) $(pdp3_OBJS) $(CLEAN_FILES)
	$(RM) $(LIBS) $(EXES) $(EXES:%=%.so)
	$(RM) -r $(BUILD_DIRS) $(ROOTDIR)/python/__pycache__/

$(SUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean


$(EXTRASUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(RM) $(CLEAN_FILES)


$(TESTSUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean

$(pdp3_MODULE): $(pdp3_OBJS)
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) $(LDFLAGS) -o $@ $(pdp3_OBJS) $(pdp3_LIBRARY_PATH) $(pdp3_LIBRARIES:%=-l%)

mrproper: clean
	$(RM) -r $(TESTDIR) *.avi

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
