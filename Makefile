##   Makefile
##
##   default variable setting
##   
MAKE=$(shell which make)
PYTHON=$(shell which python)

WHOLE_ARCHIVE = --whole_archive
NO_WHOLE_ARCHIVE = ,--no_whole_archive

SWIG=$(shell which swig)
SWIGFLAG = -Wall -c++ -python

MFEM=/usr/local/mfem-3.1
MFEMLIB = mfem
MFEMINCDIR = $(MFEM)
MFEMLNKDIR = $(MFEM)

MFEMSER=/usr/local/mfem-3.1ser
MFEMSERLIB = mfem
MFEMSERINCDIR = $(MFEMSER)
MFEMSERLNKDIR = $(MFEMSER)

# HYPRE
HYPRE=/usr/local/hypre-2.11.0
HYPRELIB = HYPRE
HYPREINCDIR = $(HYPRE)/include
HYPRELNKDIR = $(HYPRE)/lib

METIS4=/usr/local/metis
METIS4LIB = metis
METIS4LIBA   = $(METIS4)/lib/libmetis.a	

#MPI
MPIINCDIR= /opt/local/include/mpich-mp         #mpi.h
MPICHINCDIR    = /opt/local/include/mpich-mp
MPICHLNKDIR    = /opt/local/lib/mpich-mp
MPILIB = mpi
MPICC = mpicc
MPICXX = mpicxx
MPIFC = mpifort
MPIFL = mpifort
MPI4PYINCDIR = $(shell $(PYTHON) -c "import mpi4py;print mpi4py.get_include()")

#numpy
NUMPYINCDIR = $(shell $(PYTHON) -c "import numpy;print numpy.get_include()")

OUTC    = -o 
OPTF    = -O  -DALLOW_NON_INIT
OPTL    = -O 
OPTC    = -O
NOCOMPACTUNWIND = 
include ./Makefile.local

MFEMINC  = -I$(MFEMINCDIR)
MFEMSERINC  = -I$(MFEMSERINCDIR)
HYPREINC = -I$(HYPREINCDIR)
HYPRELNK = -L$(HYPRELNKDIR) -l$(HYPRELIB)
MPIINC  = -I$(MPIINCDIR)
MPI4PYINC  = -I$(MPI4PYINCDIR)

# export everything so that it is avaialbe in setup.py
export

SUBDIRS = mfem/par mfem/ser

.PHONEY:clean par ser  subdirs subdirs_cxx parcxx sercxx

default: setup_local.py 
#default: setup_local.py
all: par ser
cxx: parcxx sercxx
par: setup_local.py
	$(MAKE) -C mfem/par
ser: setup_local.py
	$(MAKE) -C mfem/ser
parcxx: setup_local.py
	$(MAKE) -C mfem/par cxx
sercxx: setup_local.py
	$(MAKE) -C mfem/ser cxx
setup_local.py: Makefile.local
	$(PYTHON) write_setup_local.py
cleancxx:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs cleancxx;\
	done
clean:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs clean;\
	done
	rm -f setup_local.py

