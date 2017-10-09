### Generated by Winemaker 0.8.4
###
### Invocation command line was
### /usr/bin/winemaker-stable .


SRCDIR                = .
SUBDIRS               =
DLLS                  =
LIBS                  =
EXES                  = pdp3



### Common settings

CEXTRA                = -mno-cygwin
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
			kern_accessor.cpp \
			particles_list.cpp \
			particles_struct.cpp \
			tinystr.cpp \
			tinyxml.cpp \
			tinyxmlerror.cpp \
			tinyxmlparser.cpp \
			wrapper.cpp
pdp3_RC_SRCS          = pdp3.rc \
			pdp31.rc
pdp3_LDFLAGS          = -mwindows \
			-mno-cygwin
pdp3_ARFLAGS          =
pdp3_DLL_PATH         =
pdp3_DLLS             = odbc32 \
			ole32 \
			oleaut32 \
			winspool \
			odbccp32
pdp3_LIBRARY_PATH     =
pdp3_LIBRARIES        = uuid

pdp3_OBJS             = $(pdp3_C_SRCS:.c=.o) \
			$(pdp3_CXX_SRCS:.cpp=.o) \
			$(pdp3_RC_SRCS:.rc=.res)



### Global source lists

C_SRCS                = $(pdp3_C_SRCS)
CXX_SRCS              = $(pdp3_CXX_SRCS)
RC_SRCS               = $(pdp3_RC_SRCS)


### Tools

CC = winegcc
CXX = wineg++
RC = wrc
AR = ar


### Generic targets

all: $(SUBDIRS) $(DLLS:%=%.so) $(LIBS) $(EXES)

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

.rc.res:
	$(RC) $(RCFLAGS) $(RCEXTRA) $(DEFINCL) -fo$@ $<

# Rules for cleaning

CLEAN_FILES     = y.tab.c y.tab.h lex.yy.c core *.orig *.rej \
                  \\\#*\\\# *~ *% .\\\#*

clean:: $(SUBDIRS:%=%/__clean__) $(EXTRASUBDIRS:%=%/__clean__)
	$(RM) $(CLEAN_FILES) $(RC_SRCS:.rc=.res) $(C_SRCS:.c=.o) $(CXX_SRCS:.cpp=.o)
	$(RM) $(DLLS:%=%.so) $(LIBS) $(EXES) $(EXES:%=%.so)

$(SUBDIRS:%=%/__clean__): dummy
	cd `dirname $@` && $(MAKE) clean

$(EXTRASUBDIRS:%=%/__clean__): dummy
	-cd `dirname $@` && $(RM) $(CLEAN_FILES)

### Target specific build rules
DEFLIB = $(LIBRARY_PATH) $(LIBRARIES) $(DLL_PATH) $(DLL_IMPORTS:%=-l%)

$(pdp3_MODULE): $(pdp3_OBJS)
	$(CXX) $(pdp3_LDFLAGS) -o $@ $(pdp3_OBJS) $(pdp3_LIBRARY_PATH) $(pdp3_DLL_PATH) $(DEFLIB) $(pdp3_DLLS:%=-l%) $(pdp3_LIBRARIES:%=-l%)


