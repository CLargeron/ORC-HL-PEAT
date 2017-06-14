#-
#- $Id: AA_make.ldef 12 2010-11-05 15:42:13Z mmaipsl $
#-
#---------------------------------------------------------------------
#- Creation des elements relatifs a PARAMETER
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
MODEL_LIB = $(LIBDIR)/libparameters.a
SXMODEL_LIB = $(MODEL_LIB)
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
#- $Id: AA_make 956 2012-07-26 14:16:59Z didier.solyga $
#-
IOIPSL_LIB = $(LIBDIR)/libioipsl.a
SXIOIPSL_LIB = $(IOIPSL_LIB)
#-
PARALLEL_LIB = $(LIBDIR)/libparallel.a
SXPARALLEL_LIB = $(PARALLEL_LIB)
#-
MODS1 = constantes.f90 \
	constantes_soil.f90   \
	constantes_mtc.f90   \
	pft_parameters.f90 
#
OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)
#-
all:
	$(M_K) libioipsl
	$(M_K) libparallel
	$(M_K) m_all
	@echo parameter is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libioipsl:
	(cd ../../IOIPSL/src; $(M_K) -f Makefile)

libparallel:
	(cd ../src_parallel; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(RM) $*.o

config:
	$(BINDIR)/Fparser -name PARAMETERS  $(MODS1)
	echo 'Configuration of PARAMETERS done'
 
clean:
	$(RM) $(MODEL_LIB)
 
$(MODEL_LIB)(pft_parameters.o): \
  $(MODEL_LIB)(constantes_mtc.o)
$(MODEL_LIB)(constantes_soil.o): \
  $(MODEL_LIB)(constantes.o)
$(MODEL_LIB)(constantes.o): \
  $(PARALLEL_LIB) \
  $(IOIPSL_LIB)

