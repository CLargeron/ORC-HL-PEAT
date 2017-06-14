#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Local parameters for "libioipsl"
#---------------------------------------------------------------------
SHELL = /bin/ksh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
P_P =
MODEL_LIB = $(LIBDIR)/libioipsl.a
#-
#- $Id: AA_make.gdef 2039 2013-04-04 08:43:44Z aclsce $
#-
#- Validate the correlation between the target and the environment
#-
UTIL_DIR = ../../../util
WW_h_t = $(shell cat $(UTIL_DIR)/.host_target)
WW_h_w = $(shell $(UTIL_DIR)/w_i_h)
WW_t_e = $(shell $(UTIL_DIR)/w_i_e $(WW_h_t) $(WW_h_w))
ifeq "$(WW_t_e)" "NO"
 $(error )
endif
#-
#- Global definitions for Mac OSX, gfortran
LIB_MPI = MPI1
LIB_MPI_BIS = MPI1
PRISM_ARCH = Darwin
PRISM_NAME = Darwin
FCM_ARCH = gfortran
M_K = make
P_C = cpp
P_O = -fpreprocessed -P -C -traditional $(P_P)
F_C = gfortran -c -cpp
F_D =
F_P = -fdefault-real-8
w_w = -O5  -funroll-all-loops $(F_D) $(F_P) -I$(MODDIR)
F_O = $(w_w) -J$(MODDIR)
F_L = gfortran
M_M = 0
L_X = 0
L_O =
A_C = ar -rs
A_G = ar -x
C_C = cc -c
C_O =
C_L = cc
FC_C = cc -c
FC_O =
FC_L = cc
#-
###-Q- Darwin   NCDF_INC = /sw/opt/netcdf7/include
###-Q- Darwin   NCDF_LIB = -L/sw/opt/netcdf7/lib -lnetcdf
NCDF_INC = /usr/local/include
NCDF_LIB = -L/usr/local/lib -lnetcdf
#-
######-Q- ada      F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
######-Q- curie  F_O = -DCPP_PARA -p -g -traceback -fp-stack-check -ftrapuv -check bounds $(F_D) $(F_P) -I$(MODDIR) -module $(MODDIR)
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
#- $Id: AA_make 105 2007-07-26 09:46:30Z mafoipsl $
#-
MAKE_NAM = $(MAKE)
ifneq ($(MAKE_NAM),$(M_K))
 @$(error You must invoke this Makefile with the $(M_K) command)
endif
USER_DIR = $(shell /bin/pwd)
MAKE_DIR = /Users/clargeron/ORCHIDEE_MICT3_Cconso_CH4/modeles/IOIPSL/src
ifneq ($(USER_DIR),$(MAKE_DIR))
 @$(error You must invoke this Makefile from its directory)
endif
#-
#---------------------------------------------------------------------
#- Create "libioipsl"
#---------------------------------------------------------------------

.PRECIOUS : $(MODEL_LIB)
SXMODEL_LIB = $(MODEL_LIB)
ifeq ($(L_X),1)
SXMODEL_LIB = $(LIBDIR)/libsxioipsl.a
.PRECIOUS : $(SXMODEL_LIB)
endif

MODS1 = defprec.f90   \
        errioipsl.f90 \
        stringop.f90  \
        mathelp.f90   \
        getincom.f90  \
        calendar.f90  \
        fliocom.f90   \
        flincom.f90   \
        histcom.f90   \
        restcom.f90   \
        ioipsl.f90
OBJSMODS1 = $(MODS1:.f90=.o)

all: $(MODEL_LIB)($(OBJSMODS1))
	@echo IOIPSL is OK

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
ifeq ($(L_X),1)
	$(A_X) $(SXMODEL_LIB) $*.o
endif
ifeq ($(M_M),1)
	mv $*.mod $(MODDIR)
endif
	$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)
ifeq ($(L_X),1)
	$(RM) $(SXMODEL_LIB)
endif
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90
	$(RM) $(MODDIR)/*.mod $(MODDIR)/*.M *.M

#- Specific dependencies

$(MODEL_LIB)(errioipsl.o):   \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(stringop.o):    \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(mathelp.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(getincom.o):    \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(calendar.o):    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(fliocom.o):     \
  $(MODEL_LIB)(defprec.o)    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(flincom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(histcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)   \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(restcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)   \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(ioipsl.o):      \
  $(MODEL_LIB)(getincom.o)   \
  $(MODEL_LIB)(fliocom.o)    \
  $(MODEL_LIB)(flincom.o)    \
  $(MODEL_LIB)(histcom.o)    \
  $(MODEL_LIB)(restcom.o)

