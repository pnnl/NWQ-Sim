### User Configurable Options

## [1] Choose your compiler
# Use Cray’s C++ wrapper (NVIDIA toolchain) with C++17 and PIC
CCCOM = g++ -std=c++17 -fPIC

# If you ever switch off of Cray you could re-enable one of these:
# GNU GCC
# CCCOM = g++ -std=c++17 -fPIC
# Clang
# CCCOM = clang++ -std=c++17 -fPIC -Wno-gcc-compat

#########
## [2] BLAS/LAPACK setup
## We’re using Cray LibSci as a LAPACK provider here.
PLATFORM = lapack

## Link flags for Cray LibSci’s BLAS+LAPACK
## (these libraries come from `module load cray-libsci`)
BLAS_LAPACK_LIBFLAGS   = -lblas -llapack
BLAS_LAPACK_INCLUDEFLAGS =   # nothing special needed for headers

#########
## [3] HDF5 support (keep commented out if you don’t need it)
#HDF5_PREFIX = /path/to/hdf5

#########
## [4] Native ITensor OpenMP (optional)
#ITENSOR_USE_OMP = 1

#########
## [5] Optimization vs Debug flags
OPTIMIZATIONS=-O3 -DNDEBUG -Wall
DEBUGFLAGS    = -g -DDEBUG -Wall -pedantic

## Build dynamic libraries? (usually 0 on Cray)
ITENSOR_MAKE_DYLIB = 0

### Internal make‑vars (no need to edit below here) ###
PREFIX           = $(THIS_DIR)
ITENSOR_LIBDIR   = $(PREFIX)/lib
ITENSOR_INCLUDEDIR = $(PREFIX)
ITENSOR_LIBNAMES = itensor

ITENSOR_LIBFLAGS  = $(patsubst %,-l%,$(ITENSOR_LIBNAMES)) $(BLAS_LAPACK_LIBFLAGS)
ITENSOR_LIBGFLAGS = $(patsubst %,-l%,-g$(ITENSOR_LIBNAMES)) $(BLAS_LAPACK_LIBFLAGS)
ITENSOR_LIBS      = $(patsubst %,$(ITENSOR_LIBDIR)/lib%.a,$(ITENSOR_LIBNAMES))
ITENSOR_GLIBS     = $(patsubst %,$(ITENSOR_LIBDIR)/lib%-g.a,$(ITENSOR_LIBNAMES))

ITENSOR_INCLUDEFLAGS = -I$(ITENSOR_INCLUDEDIR) $(BLAS_LAPACK_INCLUDEFLAGS)

ifdef HDF5_PREFIX
  ITENSOR_USE_HDF5 = 1
  ITENSOR_INCLUDEFLAGS += -I$(HDF5_PREFIX)/include -DITENSOR_USE_HDF5
  ITENSOR_LIBFLAGS     += -L$(HDF5_PREFIX)/lib -lhdf5 -lhdf5_hl
  ITENSOR_LIBGFLAGS    += -L$(HDF5_PREFIX)/lib -lhdf5 -lhdf5_hl
endif

ifndef CCCOM
  $(error CCCOM not defined in options.mk! Please set your compiler wrapper above.)
endif

ifdef ITENSOR_USE_OMP
  ITENSOR_INCLUDEFLAGS += -fopenmp -DITENSOR_USE_OMP
  ITENSOR_LIBFLAGS     += -fopenmp
  ITENSOR_LIBGFLAGS    += -fopenmp
endif

CCFLAGS  = -I. $(ITENSOR_INCLUDEFLAGS)  $(OPTIMIZATIONS)
CCGFLAGS = -I. $(ITENSOR_INCLUDEFLAGS)  $(DEBUGFLAGS)
LIBFLAGS = -L$(ITENSOR_LIBDIR) $(ITENSOR_LIBFLAGS)
LIBGFLAGS= -L$(ITENSOR_LIBDIR) $(ITENSOR_LIBGFLAGS)

# Shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  DYLIB_EXT  = dylib
  DYLIB_FLAGS = -dynamiclib
else
  DYLIB_EXT  = so
  DYLIB_FLAGS = -shared
endif


