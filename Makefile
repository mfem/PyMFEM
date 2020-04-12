##   Makefile
##
##   default variable setting
##   
MAKE=$(shell which make)
PYTHON ?= $(shell which python)
PREFIX ?= /usr/local

# serial compiler
CXX_SER ?= $(CXX)
CC_SER ?= $(CC)
# parallel compiler(CC_PAR)/linker(CXX_PAR)
CXX_PAR ?= $(MPICXX)
CC_PAR ?=  $(MPICC)

WHOLE_ARCHIVE = --whole_archive
NO_WHOLE_ARCHIVE = --no-whole-archive

SWIG=$(shell which swig)
SWIGFLAG = -Wall -c++ -python -fastproxy -olddefs -keyword

#
# MFEM path:
#
#   MFEMBUILDDIR : directory of MFEM build. Need to find config/config.hpp
#   MFEMINCDIR : include files
#   MFEMLNKDIR : path to mfem.so


MFEM ?=/usr/local
MFEMLIB = mfem
MFEMBUILDDIR ?= $(HOME)/src/mfem/cmbuild_par
MFEMINCDIR = $(MFEM)/include/mfem
MFEMLNKDIR = $(MFEM)/lib

MFEMSER ?=/usr/local/mfem_ser
MFEMSERLIB = mfem
MFEMSERBUILDDIR ?= $(HOME)/src/mfem_ser/cmbuild_ser
MFEMSERINCDIR = $(MFEMSER)/include/mfem
MFEMSERLNKDIR = $(MFEMSER)/lib

# HYPRE
HYPREINC ?= /usr/local/include
HYPRELIB ?= /usr/local/lib

#metis
METIS5INC ?= /usr/local/include
METIS5LIB ?= /usr/local/lib

#MPI
#MPICHINC  ?= /usr/local/include/mpich-mp
#MPICHLIB  ?= /usr/local/lib/mpich-mp
MPI4PYINC = $(shell $(PYTHON) -c "import mpi4py;print(mpi4py.get_include())")

#numpy
#NUMPYINC = $(shell $(PYTHON) -c "import numpy;print(numpy.get_include())")

NOCOMPACTUNWIND = 
include ./Makefile.local

MFEMINCFLAG  = -I$(MFEMINCDIR)
MFEMSERINCFLAG  = -I$(MFEMSERINCDIR)
HYPREINCFLAG = -I$(HYPREINC)
HYPRELNKFLAG = -L$(HYPRELIB) -lHYPRE
#MPIINCFLAG  = -I$(MPIINC)
MPI4PYINCFLAG  = -I$(MPI4PYINC)

ADD_STRUMPACK ?= $(ENABLE_STRUMPACK)
STRUMPACK_INCLUDE ?=
ADD_PUMI ?= $(ENABLE_PUMI)

# export everything so that it is avaialbe in setup.py
export

SUBDIRS = mfem/_par mfem/_ser

.PHONEY:clean par ser  subdirs subdirs_cxx parcxx sercxx pyinstall

default: setup_local.py 
#default: setup_local.py
all: par ser
cxx: parcxx sercxx
par: 
	$(PYTHON) write_setup_local.py
	$(MAKE) -C mfem/_par
	cp setup_local.py mfem/.       
ser: 
	$(PYTHON) write_setup_local.py
	$(MAKE) -C mfem/_ser
	cp setup_local.py mfem/.       
parcxx: setup_local.py
	$(MAKE) -C mfem/_par cxx
sercxx: setup_local.py
	$(MAKE) -C mfem/_ser cxx
setup_local.py: Makefile.local
	$(PYTHON) write_setup_local.py
	cp setup_local.py mfem/.
pyinstall:
	#$(PYTHON) clean_import.py -x
	$(PYTHON) setup.py build
	$(PYTHON) setup.py install --prefix=$(PREFIX)
cleancxx:
	for dirs in $(SUBDIRS); do\
	        if [ -d $$dirs ]; then \
	           $(MAKE) -C $$dirs cleancxx;\
		fi; \
	done
clean:
	for dirs in $(SUBDIRS); do\
	        if [ -d $$dirs ]; then \
	           $(MAKE) -C $$dirs clean;\
		fi; \
	done
	rm -f setup_local.py
	rm -rf build

