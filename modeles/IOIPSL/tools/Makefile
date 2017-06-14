#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Local parameters related to the "tools" of IOIPSL (rebuild, config)
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
#- $Id: AA_make 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Create "tools" related to IOIPSL
#---------------------------------------------------------------------

SXMODEL_LIB = $(MODEL_LIB)
ifeq ($(L_X),1)
SXMODEL_LIB = $(LIBDIR)/libsxioipsl.a
endif

rbld: libioipsl $(BINDIR)/rebuild $(BINDIR)/flio_rbld
	@echo rebuild is OK

ncregular: $(BINDIR)/ncregular 
	@echo ncregular is OK

config: libioipsl $(BINDIR)/Fparser $(BINDIR)/tkparse
	$(BINDIR)/Fparser -name IOIPSL *.f90
	echo 'Configuration of IOIPSL done'

$(BINDIR)/rebuild: ./rebuild
	cp -p ./rebuild $(BINDIR)/rebuild
	chmod a+rx $(BINDIR)/rebuild

$(BINDIR)/flio_rbld: $(SXMODEL_LIB) flio_rbld.f90
	$(F_C) $(F_O) -I$(NCDF_INC) flio_rbld.f90
	$(F_L) flio_rbld.o -o $(BINDIR)/flio_rbld \
               $(L_O) $(SXMODEL_LIB) $(NCDF_LIB)
	chmod a+rx $(BINDIR)/flio_rbld
	$(RM) flio_rbld.o

$(BINDIR)/ncregular: ncregular.f90
	$(F_C) $(F_O) -I$(NCDF_INC) ncregular.f90
	$(F_L) -o $(BINDIR)/ncregular ncregular.o $(L_O) $(NCDF_LIB)
	$(RM) ncregular.o

$(BINDIR)/Fparser: Fparser.o
	$(F_L) -o $(BINDIR)/Fparser Fparser.o -L$(LIBDIR) -lioipsl
	$(RM) Fparser.o

Fparser.o: Fparser.f90
	$(F_C) $(F_O) Fparser.f90

$(BINDIR)/tkparse: tkparse.o tkcond.o tkgen.o
	$(C_L) -o $(BINDIR)/tkparse tkparse.o tkcond.o tkgen.o
	$(RM) tkparse.o tkcond.o tkgen.o

tkparse.o: tkparse.c tkparse.h
	$(C_C) $(C_O) -o tkparse.o tkparse.c

tkcond.o: tkcond.c tkparse.h
	$(C_C) $(C_O) -o tkcond.o tkcond.c

tkgen.o: tkgen.c tkparse.h
	$(C_C) $(C_O) -o tkgen.o tkgen.c

libioipsl:
	(cd ../src; $(M_K) -f Makefile)

clean:
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90 *.M
	$(RM) $(BINDIR)/rebuild $(BINDIR)/flio_rbld
	$(RM) $(BINDIR)/ncregular $(BINDIR)/Fparser $(BINDIR)/tkparse

