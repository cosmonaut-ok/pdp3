SRCDIR                = .
SUBDIRS               =
DLLS                  =
LIBS                  =
EXES                  = pdp3
TESTDIR               = testdir
### Common settings
CEXTRA                =
CXXEXTRA              =
RCEXTRA               =
DEFINES               =
INCLUDE_PATH          = -I.
DLL_PATH              =
DLL_IMPORTS           =
LIBRARY_PATH          =
LIBRARIES             =

### pdp3 sources and settings
pdp3_MODULE           = pdp3
pdp3_C_SRCS           =
pdp3_CXX_SRCS         = Beam.cpp \
			Boundary_Maxwell_conditions.cpp \
			Boundary_conditions.cpp \
			Bunch.cpp \
			E_field.cpp \
			Fourier.cpp \
			Geometry.cpp \
			H_field.cpp \
			Load_init_param.cpp \
			Main.cpp \
			PML.cpp \
			Particles.cpp \
			Poisson.cpp \
			Poisson_dirichlet.cpp \
			Poisson_neumann.cpp \
			Test_class.cpp \
			Time.cpp \
			Triple.cpp \
			boundary_dirichlet.cpp \
			boundary_neumann.cpp \
			charge_density.cpp \
			current.cpp \
			field.cpp \
			input_output_class.cpp \
			particles_list.cpp \
			particles_struct.cpp \
			tinyxml2.cpp
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

### Global source lists
C_SRCS                = $(pdp3_C_SRCS)
CXX_SRCS              = $(pdp3_CXX_SRCS)
RC_SRCS               = $(pdp3_RC_SRCS)

### Tools
CC = gcc
CXX = g++
# CC = winegcc
# CXX = wineg++
RC = wrc
AR = ar
CFLAGS = -O2 -Wall
CXXFLAGS = ${CFLAGS}

### Generic targets
all: $(SUBDIRS) $(DLLS:%=%.so) $(LIBS) $(EXES) bootstrap

### Build rules
.PHONY: all clean dummy

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

clean:: $(SUBDIRS:%=%/__clean__) $(EXTRASUBDIRS:%=%/__clean__)
	$(RM) $(CLEAN_FILES) $(C_SRCS:.c=.o) $(CXX_SRCS:.cpp=.o)
	$(RM) $(DLLS:%=%.so) $(LIBS) $(EXES) $(EXES:%=%.so)

$(SUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean

$(EXTRASUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(RM) $(CLEAN_FILES)

### Target specific build rules
DEFLIB = $(LIBRARY_PATH) $(LIBRARIES) $(DLL_PATH) $(DLL_IMPORTS:%=-l%)

$(pdp3_MODULE): $(pdp3_OBJS)
	$(CXX) $(pdp3_LDFLAGS) -o $@ $(pdp3_OBJS) $(pdp3_LIBRARY_PATH) $(pdp3_DLL_PATH) $(DEFLIB) $(pdp3_DLLS:%=-l%) $(pdp3_LIBRARIES:%=-l%)

bootstrap:
	mkdir -p pdp3_files pdp3_result/Dump

mrproper: clean
	rm -rf pdp3_files pdp3_result $(TESTDIR)

run: bootstrap
	./pdp3

test: mrproper all
	TESTDIR=$(TESTDIR) /bin/bash ./test.sh

.PHONY: check-syntax
check-syntax:
	$(CXX) $(LIBRARY_PATH) $(INCLUDE_PATH) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)
