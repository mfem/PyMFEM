IFILE = $(wildcard *.i)
IFILE := $(filter-out hypre_int.i, $(IFILE))
ifeq ($(ADD_STRUMPACK), )
   IFILE := $(filter-out strumpack.i, $(IFILE))
endif

ALLCXX = $(IFILE:.i=_wrap.cxx)

.PHONY: all clean cxx 
ifeq ($(MFEM),)
$(info Parallel version is skipped)

all:
default:
cxx:
allso:
%_wrap.cxx: %.i

else

defaut: allso

all: allso

cxx: $(ALLCXX)

allso: 
	$(PYTHON) setup.py build_ext --inplace
%_wrap.cxx: %.i
	$(SWIG) -I$(MFEMBUILDDIR) $(SWIGFLAG) $(MPI4PYINCFLAG) $(MFEMINCFLAG) -I$(STRUMPACK_INCLUDE) $<
endif
cleancxx:
	rm -f *.cxx
clean:
	rm -r -f build
	rm -f *.so
