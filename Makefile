##   Makefile
##
##   default variable setting
##   
MAKE=$(shell which make)
PYTHON=$(shell which python)
INSTALL_PREFIX=/usr/local/lib

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

#
# MFEM path:
#
#   MFEMBUILDDIR : directory of MFEM build. Need to find config/config.hppx
#   MFEMINCDIR : include files
#   MFEMLNKDIR : path to mfem.so

MFEM=/usr/local
MFEMLIB = mfem
MFEMBUILDDIR = $(HOME)/src/mfem
MFEMINCDIR = $(MFEM)/include
MFEMLNKDIR = $(MFEM)/lib

MFEMSER=/usr/local/mfem_ser
MFEMSERLIB = mfem
MFEMSERBUILDDIR = $(HOME)/src/mfem_ser
MFEMSERINCDIR = $(MFEMSER)/include
MFEMSERLNKDIR = $(MFEMSER)/lib

# HYPRE
HYPRE=/usr/local
HYPRELIB = HYPRE
HYPREINCDIR = $(HYPRE)/include
HYPRELNKDIR = $(HYPRE)/lib

#metis
# METISLIB will become -lmetis
# METISLNKDIR will become -L<dir>
# overwrite METISLIBA to black in Makefile.local if metis is provide as .so
METIS=/usr/local
METISLIB = metis
METISINCDIR = $(METIS)/include
METISLNKDIR = $(METIS)/lib
#METISLIBA   = $(METIS)/libmetis.a 

#MPI
MPIINCDIR= /usr/local/include/mpich-mp         #mpi.h
MPICHINCDIR    = /usr/local/include/mpich-mp
MPICHLNKDIR    = /usr/local/lib/mpich-mp
MPILIB = mpi
MPI4PYINCDIR = $(shell $(PYTHON) -c "import mpi4py;print mpi4py.get_include()")

#numpy
NUMPYINCDIR = $(shell $(PYTHON) -c "import numpy;print numpy.get_include()")

#Boost
BOOSTINCDIR = /usr/local/include
BOOSTLIBDIR = /usr/local/lib
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

.PHONEY:clean par ser  subdirs subdirs_cxx parcxx sercxx pyinstall

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
pyinstall:
	$(PYTHON) setup.py build
	$(PYTHON) setup.py install --prefix=$(INSTALL_PREFIX)
cleancxx:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs cleancxx;\
	done
clean:
	for dirs in $(SUBDIRS); do\
		$(MAKE) -C $$dirs clean;\
	done
	rm -f setup_local.py

