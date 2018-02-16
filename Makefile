SRCDIR                = .
SUBDIRS               =
TESTSUBDIRS           = test/unit
DLLS                  =
LIBS                  =
EXES                  = pdp3
TESTDIR               = testdir
### Common settings
CEXTRA                =
CXXEXTRA              =
RCEXTRA               =
DEFINES               =
INCLUDE_PATH          = -I$(SRCDIR)
DLL_PATH              =
DLL_IMPORTS           =
LIBRARY_PATH          =
LIBRARIES             =

### tinyxml2 sources and settings
tinyxml2_SUBDIR := $(SRCDIR)/lib/tinyxml2/
# LDFLAGS := -L$(tinyxml2_SUBDIR)
INCLUDE_PATH += -I$(tinyxml2_SUBDIR)
SUBDIRS += $(tinyxml2_SUBDIR)
DLL_IMPORTS := tinyxml2
LDFLAGS += -L$(tinyxml2_SUBDIR)

### pdp3 sources and settings
pdp3_MODULE           = pdp3
pdp3_C_SRCS           =
pdp3_CXX_SRCS         = Boundary_Maxwell_conditions.cpp \
			Bunch.cpp \
			E_field.cpp \
			Fourier.cpp \
			Geometry.cpp \
			H_field.cpp \
			Load_init_param.cpp \
			Main.cpp \
			Particles.cpp \
			Poisson.cpp \
			Poisson_dirichlet.cpp \
			Poisson_neumann.cpp \
			pdp3_time.cpp \
			Triple.cpp \
			charge_density.cpp \
			current.cpp \
			field.cpp \
			input_output_class.cpp \
			particles_list.cpp \
			particles_struct.cpp

pdp3_RC_SRCS          = pdp3.rc \
			pdp31.rc
pdp3_LDFLAGS          =
pdp3_ARFLAGS          =
pdp3_DLL_PATH         =
pdp3_DLLS             =
pdp3_LIBRARY_PATH     =
pdp3_LIBRARIES        =

pdp3_OBJS             = $(pdp3_C_SRCS:.c=.o) \
			$(pdp3_CXX_SRCS:.cpp=.o)

LDFLAGS := $(LDFLAGS) $(pdp3_LDFLAGS)

### Global source lists
C_SRCS                = $(pdp3_C_SRCS)
CXX_SRCS              = $(pdp3_CXX_SRCS)
RC_SRCS               = $(pdp3_RC_SRCS)

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
all: $(SUBDIRS) $(DLLS:%=%.so) $(LIBS) $(EXES)

### Build rules
.PHONY: all clean dummy check-syntax

$(SUBDIRS): dummy
	@cd $@ && $(MAKE)

# Implicit rules
.SUFFIXES: .cpp .cxx .rc .res
DEFINCL = $(INCLUDE_PATH) $(DEFINES) $(OPTIONS)

.c.o:
	$(CC) -c $(CFLAGS) $(CEXTRA) $(DEFINCL) -o $@ $<

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) -o $@ $<

# Rules for cleaning
CLEAN_FILES     = y.tab.c y.tab.h lex.yy.c core *.orig *.rej \
                  \\\#*\\\# *~ *% .\\\#*

clean:: $(SUBDIRS:%=%/__clean__) $(EXTRASUBDIRS:%=%/__clean__) $(TESTSUBDIRS:%=%/__clean__)
	$(RM) $(CLEAN_FILES) $(C_SRCS:.c=.o) $(CXX_SRCS:.cpp=.o)
	$(RM) $(DLLS:%=%.so) $(LIBS) $(EXES) $(EXES:%=%.so)
	$(RM) -r python/__pycache__/

$(SUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean


$(EXTRASUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(RM) $(CLEAN_FILES)


$(TESTSUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean

### Target specific build rules
DEFLIB = $(LIBRARY_PATH) $(LIBRARIES) $(DLL_PATH) $(DLL_IMPORTS:%=-l%)

$(pdp3_MODULE): $(pdp3_OBJS)
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) $(DEFINCL) $(LDFLAGS) -o $@ $(pdp3_OBJS) $(pdp3_LIBRARY_PATH) $(pdp3_DLL_PATH) $(DEFLIB) $(pdp3_DLLS:%=-l%) $(pdp3_LIBRARIES:%=-l%)

bootstrap: all
	mkdir -p pdp3_files pdp3_result/Dump

mrproper: clean
	$(RM) -r pdp3_files pdp3_result $(TESTDIR) *.avi

run: bootstrap
	./pdp3

$(TESTSUBDIRS:%=%/__test__): dummy
	cd `dirname $@` && $(MAKE) test

unittest: $(TESTSUBDIRS:%=%/__test__)

test: mrproper all
	TESTDIR=$(TESTDIR) /bin/bash ./test/functional/test.sh

test-ext: mrproper all
	TESTDIR=$(TESTDIR) /bin/bash ./test/functional/test.sh extended

check-syntax:
	$(CXX) $(LIBRARY_PATH) $(INCLUDE_PATH) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)
