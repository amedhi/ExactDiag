#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./make_options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/mpi_comm.cpp 
SRCS+= scheduler/cmdargs.cpp 
SRCS+= scheduler/inputparams.cpp 
SRCS+= scheduler/taskparams.cpp 
SRCS+= scheduler/worker.cpp 
SRCS+= scheduler/master_scheduler.cpp
SRCS+= scheduler/scheduler.cpp
#SRCS+= xml/pugixml.cpp 
SRCS+= expression/complex_expression.cpp
SRCS+= utils/utils.cpp 
SRCS+= lattice/lattice.cpp
SRCS+= lattice/latticelibrary.cpp
SRCS+= lattice/graph.cpp
SRCS+= model/hilbertspace.cpp
SRCS+= model/strmatrix.cpp
SRCS+= model/hamiltonian_term.cpp
SRCS+= model/model.cpp
SRCS+= model/modellibrary.cpp
SRCS+= diag/datafile.cpp
SRCS+= diag/basis.cpp
SRCS+= diag/diag.cpp
SRCS+= main.cpp
PRJ_SRCS = $(addprefix src/,$(SRCS))
#-------------------------------------------------------------
# Headers
HDRS=    scheduler/mpi_comm.h \
         scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         expression/complex_expression.h \
         utils/utils.h \
         utils/curve_fit.h \
         lattice/constants.h lattice/lattice.h lattice/graph.h \
         model/hilbertspace.h \
         model/modelparams.h  model/quantum_op.h \
	 model/strmatrix.h \
	 model/hamiltonian_term.h \
	 model/model.h \
	 diag/datafile.h \
	 diag/basis.h \
	 diag/diag.h \
	 vmc/random.h  vmc/basisstate.h \
	 main.h \
PRJ_HDRS = $(addprefix src/,$(HDRS))
MUPARSER_LIB = $(PROJECT_ROOT)/src/expression/muparserx/libmuparserx.a
#-------------------------------------------------------------
TAGT=a.out

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
  $(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(PRJ_SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 
# compiler flags

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS) $(MUPARSER_LIB)
	$(PRJ_CXX) -o $(TAGT) $(OBJS) $(PRJ_LDFLAGS) $(PRJ_LIBS) $(MUPARSER_LIB)


%.o: %.cpp
	$(PRJ_CXX) -c $(PRJ_CXXFLAGS) -o $@ $<

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(PRJ_CXX) -c $(PRJ_CXXFLAGS) -o $(@F) $(<F)"
	@$(PRJ_CXX) -MMD -c $(PRJ_CXXFLAGS) -o $@ $<

$(PRJ_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

$(MUPARSER_LIB):
	@cd ./src/expression/muparserx/ && $(MAKE)

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(PRJ_LIBDIR) and $(PRJ_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
