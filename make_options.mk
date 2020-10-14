# User Configurable Options
PROJECT_ROOT=/Users/amedhi/Projects/Codes/ExactDiag

#-------------------------------------------------------------
# need mpi version
#MPI=HAVE_BOOST_MPI

#-------------------------------------------------------------
# 1. Set compiler option
PRJ_CXX=clang++ -std=c++11 # Clang compiler 
#PRJ_CXX=g++ -std=c++11 # GNU GCC compiler

ifeq ($(MPI), HAVE_BOOST_MPI)
PRJ_CXX=/opt/openmpi-2.0.2/bin/mpicxx -std=c++11
PRJ_CPPFLAGS=-D$(MPI)
endif

#-------------------------------------------------------------
DEBUG_MODE=YES
ifeq ($(DEBUG_MODE), YES)
PRJ_OPTFLAGS= -Wall -O2
else
PRJ_OPTFLAGS= -Wall -O3
PRJ_CPPFLAGS=-DEIGEN_NO_DEBUG -DEIGEN_USE_MKL_ALL
endif

# other preprocessor directiives
#-------------------------------------------------------------
# 2. Compile flags 
# Flags to give the compiler for "release mode"
#PRJ_OPTFLAGS=-Wall -pg
#PRJ_OPTFLAGS=-DDEBUG_MODE -g -Wall -pedantic
# Flags to give the compiler for "debug mode"
#PRJ_DEBUGFLAGS=-DDEBUG_MODE -g -Wall -pedantic

#-------------------------------------------------------------
# 3. Boost and Eigen library
# Flags to give the compiler for "release mode"
BOOST_INCLUDE=-I/usr/local/include
EIGEN_INCLUDE=-I/usr/local/include

# Boost MPI library
ifeq ($(MPI), HAVE_BOOST_MPI)
BOOST_LIBS=-lboost_filesystem -lboost_system -lboost_mpi -lboost_serialization
BOOST_LDFLAGS=-L/usr/local/lib
else 
BOOST_LIBS=-lboost_filesystem -lboost_system
BOOST_LDFLAGS=-L/usr/local/lib
endif

#EIGEN_USE_MKL=EIGEN_USE_MKL_ALL
MKL_INCLUDE=-I/opt/intel/mkl/include/intel64/lp64
MKL_LDFLAGS=-L/opt/intel/mkl/lib
MKL_LIBS=-lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

#NLOPT_INCLUDE=-I/Users/amedhi/projects/Codes/vmc++/libs/include
#NLOPT_LDFLAGS=-L/Users/amedhi/projects/Codes/vmc++/libs/lib
#NLOPT_LIBS=-lnlopt

INCLUDE = $(BOOST_INCLUDE) #$(MKL_INCLUDE)
ifneq ($(BOOST_INCLUDE), $(EIGEN_INCLUDE))
INCLUDE += $(EIGEN_INCLUDE)
endif

#ifeq ($(EIGEN_USE_MKL), USE_INTEL_MKL_ALL)
#INCLUDE += $(MKL_INCLUDE)
#endif
PRJ_CXXFLAGS=$(PRJ_CPPFLAGS) $(PRJ_OPTFLAGS) $(INCLUDE) #$(NLOPT_INCLUDE)
PRJ_LDFLAGS=$(BOOST_LDFLAGS) #$(NLOPT_LDFLAGS)  #$(MKL_LDFLAGS)
PRJ_LIBS=$(BOOST_LIBS) $(MKL_LIBS) #$(NLOPT_LIBS)  #$(MKL_LIBS)
#ifeq ($(EIGEN_USE_MKL), USE_INTEL_MKL_ALL)
#PRJ_CXXBFLAGS += -D$(EIGEN_USE_MKL)
#PRJ_BLDFLAGS += $(MKL_LDFLAGS)
#PRJ_BLIBS += $(MKL_LIBS)
#endif

#-------------------------------------------------------------
# 4. Where to put the 'cmc' library & the includes
PREFIX=$(PROJECT_ROOT)
BUILD_DIR=$(PREFIX)/build
PRJ_LIBDIR=$(PREFIX)/lib
PRJ_INCLUDE=$(PREFIX)/include
#PRJ_CXXFLAGS= $(PRJ_OPTFLAGS) $(INCLUDE) -I$(PRJ_INCLUDE)
#PRJ_LDFLAGS=$(BOOST_LDFLAGS) -L$(PRJ_LIBDIR)
#PRJ_LIBS=$(BOOST_LIBS) -lvmc++

#-------------------------------------------------------------
