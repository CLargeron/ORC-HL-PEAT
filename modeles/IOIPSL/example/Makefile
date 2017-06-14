#-
#- $Id: AA_make.ldef 53 2007-06-06 15:27:26Z bellier $
#-
#---------------------------------------------------------------------
#- Samples related to IOIPSL
#---------------------------------------------------------------------
SHELL = /bin/ksh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
SCPDIR = ../../../scripts
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
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
#- $Id: AA_make 107 2007-07-26 10:10:41Z bellier $
#-
USER_DIR = $(shell /bin/pwd)
MAKE_DIR = /Users/clargeron/ORCHIDEE_MICT3_Cconso_CH4/modeles/IOIPSL/example
ifneq ($(USER_DIR),$(MAKE_DIR))
 @$(error You must invoke this Makefile from its directory)
endif
#-
IOIPSL_LIB = -L$(LIBDIR) -lioipsl
ifeq ($(L_X),1)
IOIPSL_LIB = -L$(LIBDIR) -lsxioipsl
endif
#-
all: testhist1 testhist2 testrest testopp testcalendar testflio
	echo don t forget to do : $(M_K) clean !!!

testhist1:
	$(M_K) libioipsl
	$(M_K) testhist1.x
	./testhist1.x \

	@echo ' '
	@echo 'file testhist1.nc should be created'
	@echo 'to see it, type : ncdump -h testhist1.nc'
	@echo ' '

testhist1.x: testhist1.o
	$(F_L) -o testhist1.x testhist1.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testhist2:
	$(M_K) libioipsl
	$(M_K) testhist2.x
	./testhist2.x \

	@echo ' '
	@echo 'file testhist2.nc should be created'
	@echo 'to see it, type : ncdump -h testhist2.nc'
	@echo ' '
	@echo 'file testhist2_bis.nc should be created'
	@echo 'to see it, type : ncdump -h testhist2_bis.nc'
	@echo ' '

testhist2.x: testhist2.o
	$(F_L) -o testhist2.x testhist2.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testrest:
	$(M_K) libioipsl
	$(M_K) testrest.x
	$(RM) restfile.nc restfilebis.nc
	./testrest.x \

	@echo ' '
	@echo 'file restfile.nc should be created'
	@echo 'to see it, type : ncdump -h restfile.nc'
	@echo ' '
	@echo 'file restfilebis.nc should be created'
	@echo 'to see it, type : ncdump -h restfilebis.nc'
	@echo ' '

testrest.x: testrest.o
	$(F_L) -o testrest.x testrest.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testopp:
	$(M_K) libioipsl
	$(M_K) testopp.x
	./testopp.x \

testopp.x: testopp.o
	$(F_L) -o testopp.x testopp.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testcalendar:
	$(M_K) libioipsl
	$(M_K) testcalendar.x
	echo; \
        ./testcalendar.x; \
        echo; \
        ./testcalendar.x noleap; \
        echo; \
        ./testcalendar.x 360d; \

testcalendar.x: testcalendar.o
	$(F_L) -o testcalendar.x testcalendar.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testflio:
	$(M_K) libioipsl
	$(M_K) testflio.x

	./testflio.x \

	@echo ' '
	@echo 'file testflio_0001.nc should be created'
	@echo 'to see it, type : ncdump -h testflio_0001.nc'
	@echo ' '

testflio.x: testflio.o
	$(F_L) -o testflio.x testflio.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testconfig:
	$(M_K) libioipsl
	$(M_K) testconfig.x
	./testconfig.x \

testconfig.x: testconfig.o
	$(F_L) -o testconfig.x testconfig.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

testconfig2:
	$(M_K) libioipsl
	$(M_K) testconfig2.x
	./testconfig2.x \

testconfig2.x: testconfig2.o
	$(F_L) -o testconfig2.x testconfig2.o \
               $(L_O) $(IOIPSL_LIB) $(NCDF_LIB)

config: libioipsl
	(cd ../src; $(M_K) -f Makefile config)
	$(BINDIR)/Fparser -main -name testconfig testconfig.f90
	$(SCPDIR)/build_tkmenus
	./config.tk
	$(RM) Config.*

libioipsl:
	(cd ../src; $(M_K) -f Makefile)

clean:
	$(RM) *.*~ Makefile~ core *.o *.nc Out_t* Job_t* *.x *_run.def

testhist1.o: testhist1.f90
	$(F_C) $(F_O) testhist1.f90

testhist2.o: testhist2.f90
	$(F_C) $(F_O) testhist2.f90

testrest.o: testrest.f90
	$(F_C) $(F_O) testrest.f90

testopp.o: testopp.f90
	$(F_C) $(F_O) testopp.f90

testcalendar.o: testcalendar.f90
	$(F_C) $(F_O) testcalendar.f90

testflio.o: testflio.f90
	$(F_C) $(F_O) testflio.f90

testconfig.o: testconfig.f90
	$(F_C) $(F_O) testconfig.f90

testconfig2.o: testconfig2.f90
	$(F_C) $(F_O) testconfig2.f90

