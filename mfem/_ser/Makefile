IFILE = $(wildcard *.i)
ALLCXX = $(IFILE:.i=_wrap.cxx)
# uncomment if you want to do always make clean + make all
#.PHONY: all clean
.PHONEY:clean cxx allso

ifeq ($(MFEMSER),)
$(info Serial version is skipped)
all:
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
	$(SWIG) -I$(MFEMSERBUILDDIR) $(SWIGFLAG) $(MFEMSERINCFLAG) $<
endif
cleancxx:
	rm -f *.cxx
clean:
	rm -r -f build
	rm -f *.so
