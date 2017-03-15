##   Makefile
##
##   default variable setting
##   
MAKE=$(shell which make)
PYTHON=$(shell which python)

# serial compiler
CXX_SER = g++
CC_SER = g++

# parallel compiler(CC_PAR)/linker(CXX_PAR)
CXX_PAR = mpicxx
CC_PAR = g++

WHOLE_ARCHIVE = --whole_archive
NO_WHOLE_ARCHIVE = --no-whole-archive

SWIG=$(shell which swig)
SWIGFLAG = -Wall -c++ -python

MFEM=/usr/local/mfem-3.2
MFEMLIB = mfem
MFEMINCDIR = $(MFEM)
MFEMLNKDIR = $(MFEM)

MFEMSER=/usr/local/mfem-3.2ser
MFEMSERLIB = mfem
MFEMSERINCDIR = $(MFEMSER)
MFEMSERLNKDIR = $(MFEMSER)

# HYPRE
HYPRE=/usr/local/hypre-2.11.0
HYPRELIB = HYPRE
HYPREINCDIR = $(HYPRE)/include
HYPRELNKDIR = $(HYPRE)/lib

#metis
# METISLIB will become -lmetis
# METISLNKDIR will become -L<dir>
# overwrite METISLIBA to black in Makefile.local if metis is provide as .so
METIS=/usr/local/
METISLIB = metis
METISLNKDIR = $(METIS)/lib/
METISLIBA   = $(METIS)/libmetis.a 

#MPI
MPIINCDIR= /opt/local/include/mpich-mp         #mpi.h
MPICHINCDIR    = /opt/local/include/mpich-mp
MPICHLNKDIR    = /opt/local/lib/mpich-mp
MPILIB = mpi
MPI4PYINCDIR = $(shell $(PYTHON) -c "import mpi4py;print mpi4py.get_include()")

#numpy
NUMPYINCDIR = $(shell $(PYTHON) -c "import numpy;print numpy.get_include()")

#Boost
BOOSTINCDIR = /opt/local/include
BOOSTLIBDIR = /opt/local/lib
BOOSTLIB = boost_iostream-mt

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
	cp setup_local.py mfem/.       
ser: setup_local.py
	$(MAKE) -C mfem/ser
	cp setup_local.py mfem/.       
parcxx: setup_local.py
	$(MAKE) -C mfem/par cxx
sercxx: setup_local.py
	$(MAKE) -C mfem/ser cxx
setup_local.py: Makefile.local
	$(PYTHON) write_setup_local.py
	cp setup_local.py mfem/.       
cleancxx:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs cleancxx;\
	done
clean:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs clean;\
	done
	rm -f setup_local.py

