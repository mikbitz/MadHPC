# Repast HPC Manual Build Makefile

# VARIABLES (Supply values for these; for definitions and examples, see INSTALL)
CXX=mpicxx
CXXLD=mpicxx
BASE_DIR=$(HOME)/repastHPC
BOOST_INCLUDE_DIR=$(BASE_DIR)/Boost/Boost_1.61/include
BOOST_LIB_DIR=$(BASE_DIR)/Boost/Boost_1.61/lib
BOOST_INFIX=-mt
NETCDF_INCLUDE_DIR=$(BASE_DIR)/NetCDF/include
NETCDF_LIB_DIR=$(BASE_DIR)/NetCDF/lib64
NETCDF_CXX_INCLUDE_DIR=$(BASE_DIR)/NetCDF-cxx/include
NETCDF_CXX_LIB_DIR=$(BASE_DIR)/NetCDF-cxx/lib64
CURL_INCLUDE_DIR=$(BASE_DIR)/CURL/include
CURL_LIB_DIR=$(BASE_DIR)/CURL/lib64
INSTALL_DIR=$(BASE_DIR)/repast_hpc-2.2.0

# change to true to build static libraries
STATIC=false

# DIRECTORIES
SRC_DIR=../src
REPAST_HPC_DIR=../src/repast_hpc

# DERIVED VARIABLES (Do not modify)
I_REPAST=-I$(REPAST_HPC_DIR)
I_BOOST=-I$(BOOST_INCLUDE_DIR)
L_BOOST=-L$(BOOST_LIB_DIR)
l_BOOST=-lboost_mpi$(BOOST_INFIX) -lboost_serialization$(BOOST_INFIX) -lboost_system$(BOOST_INFIX) -lboost_filesystem$(BOOST_INFIX)
I_NETCDF=-I$(NETCDF_INCLUDE_DIR) -I$(NETCDF_CXX_INCLUDE_DIR)
L_NETCDF=-L$(NETCDF_LIB_DIR) -L$(NETCDF_CXX_LIB_DIR)
l_NETCDF=-lnetcdf_c++ -lnetcdf
I_CURL=-I$(CURL_INCLUDE_DIR)
L_CURL=-L$(CURL_LIB_DIR)
l_CURL=-lcurl

RPATHS :=
ifneq ($(NETCDF_LIB_DIR),)
  RPATHS += -Wl,-rpath -Wl,$(NETCDF_LIB_DIR)
endif

ifneq ($(NETCDF_CXX_LIB_DIR),)
  RPATHS += -Wl,-rpath -Wl,$(NETCDF_CXX_LIB_DIR)
endif

ifneq ($(BOOST_LIB_DIR),)
  RPATHS += -Wl,-rpath -Wl,$(BOOST_LIB_DIR)
endif

# These are are set in the modules included below
REPAST_HPC_SRC :=

include ../src/repast_hpc/module.mk

# SOURCE FILES
# Repast HPC
REPAST_HPC_SOURCES=$(addprefix ../, $(REPAST_HPC_SRC))

# OBJECT FILES
REPAST_HPC_OBJECTS=$(patsubst %.cpp,%.o,$(REPAST_HPC_SOURCES))

REPAST_VERSION=2.2.0
REPAST_HPC_NAME=repast_hpc-$(REPAST_VERSION)
l_REPAST_HPC=-l$(REPAST_HPC_NAME)

DYN_FLAG = -shared -fPIC
LIB_EXT = so
LIB_CXXLD = $(CXXLD)
LIB_L_FLAGS = $(L_BOOST) $(L_NETCDF) $(L_CURL)
LIB_l_FLAGS = $(l_BOOST) $(l_NETCDF) $(l_CURL)

# need to able to set this to empty for the static compile
OUTPUT_FLAG = -o

EXE_FLAG =
ifeq ($(STATIC), true)
	DYN_FLAG =
	LIB_EXT = a
	LIB_CXXLD = ar rcs
	LIB_L_FLAGS =
	LIB_l_FLAGS =
	
	RELOGO_L_FLAGS =
	RELOGO_l_FLAGS =
	OUTPUT_FLAG =
	RPATHS =
	RELOGO_RPATHS =
	EXE_FLAG = -static
else
	UNAME := $(shell uname)
	ifeq ($(UNAME), Darwin)
		DYN_FLAG = -dynamiclib
		LIB_EXT = dylib
	endif
endif

# TARGETS
REPAST_HPC_LIB=$(INSTALL_DIR)/lib/lib$(REPAST_HPC_NAME).$(LIB_EXT)

# Default rule makes everything
.PHONY : all
all: copy_includes repast_hpc

.PHONY : copy_includes
copy_includes:
	mkdir -p $(INSTALL_DIR)/include/repast_hpc
	cp $(REPAST_HPC_DIR)/*.h $(INSTALL_DIR)/include/repast_hpc/

$(REPAST_HPC_LIB): $(REPAST_HPC_OBJECTS)
	mkdir -p $(INSTALL_DIR)/lib
	rm -f $(REPAST_HPC_LIB)
	$(LIB_CXXLD) $(DYN_FLAG) $(LIB_L_FLAGS) $(OUTPUT_FLAG) $(REPAST_HPC_LIB) $(REPAST_HPC_OBJECTS) $(RPATHS)

# User-friendly names
.PHONY : repast_hpc
repast_hpc: $(REPAST_HPC_LIB)

# Cleanup
.PHONY : repast_hpc_clean
repast_hpc_clean:
	rm -f $(REPAST_HPC_LIB)
	rm -f $(REPAST_HPC_DIR)/*.o

.PHONY : clean
clean: repast_hpc_clean

.PHONY : distclean
distclean: clean
	rm -fr ../include
	rm -fr ../lib
	rm -fr ../bin


# Object files are built from cpp files
%.o:
	$(CXX) -fPIC -std=c++11 -c $(I_BOOST) $(I_NETCDF) $(I_CURL) -I$(SRC_DIR) -g0 -O3 -o $@ $*.cpp

# Cancel implicit rule for name with no suffix, to avoid confusing output:
%: %.o
