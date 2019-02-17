REPAST_HPC_VERSION = 2.2.0

ENABLE_SHARED = @ENABLE_SHARED@
ENABLE_STATIC = @ENABLE_STATIC@


LDFLAGS :=  -L/home/mb425/repastHPC/HDF/lib64 -L/home/mb425/repastHPC/NetCDF/lib64 -L/home/mb425/repastHPC/NetCDF-cxx/lib64 -L/home/mb425/repastHPC/CURL/lib -L/home/mb425/repastHPC/Boost/Boost_1.61/lib
CPPFLAGS := -I/home/mb425/repastHPC/HDF/include -I/home/mb425/repastHPC/NetCDF/include -I/home/mb425/repastHPC/NetCDF-cxx/include -I/home/mb425/repastHPC/CURL/include -I/home/mb425/repastHPC/Boost/Boost_1.61/include


#LDFLAGS :=   -L/home/mb425/repastHPC/NetCDF-cxx/lib -L/home/mb425/repastHPC/NetCDF/lib -L/home/mb425/repastHPC/CURL/lib -L/home/mb425/repastHPC/Boost/Boost_1.61/lib
#CPPFLAGS :=  -I/home/mb425/repastHPC/CURL/include -I/home/mb425/repastHPC/Boost/Boost_1.61/include

#LDFLAGS :=  
#CPPFLAGS :=  


CXX = mpicxx
CXXLD = mpicxx
CXXFLAGS = -Wall -g -O2 -std=c++11 -MMD -MP $(CPPFLAGS) -Wno-reorder -Wno-unused-variable -Wno-sign-compare -I./src

LIB_LD_FLAGS = $(LDFLAGS)
LIBS = -lboost_mpi-mt -lboost_serialization-mt -lboost_filesystem-mt -lboost_system-mt
LIBS +=-lhdf5 -lhdf5_hl -lnetcdf -lnetcdf_c++4

NETCDF_LIB_DIR = /home/mb425/repastHPC/NetCDF/lib
NETCDF_CXX_LIB_DIR = /home/mb425/repastHPC/NetCDF-cxx/lib
BOOST_LIB_DIR = /home/mb425/repastHPC/Boost/Boost_1.61/lib

AR      = ar
ARFLAGS = cr
RANLIB  = ranlib

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

USE_MAC = no
SO_SUFFIX=so
ifeq ($(USE_MAC),no)
	LIB_CXX_FLAGS += -fPIC
 	LIB_LD_FLAGS += -shared
else
	LIB_LD_FLAGS += -dynamiclib
	SO_SUFFIX=dylib
endif

HAVE_CP_U = yes
ifeq ($(HAVE_CP_U),yes)
	CP_ARGS = -uvf
else
	CP_ARGS = -vf
endif

ifeq ($(V),)
  # Prints a short description of the action, does not show command
  Q=@echo
  E=@
else
ifeq ($(V),1)
  # Echoes the entire command
  Q=@echo >/dev/null
  E=
else # V=2
  # Echoes the description and the command
  Q=@echo ; echo
  E=
endif
endif

INSTALL_PREFIX = /home/mb425/repastHPC/repast_hpc-2.2.0
INSTALL_BIN = $(INSTALL_PREFIX)/bin
INSTALL_INCLUDE = $(INSTALL_PREFIX)/include
REPAST_HPC_INSTALL_INCLUDE = $(INSTALL_INCLUDE)/repast_hpc
INSTALL_LIB = $(INSTALL_PREFIX)/lib

REPAST_HPC_LIB := repast_hpc
REPAST_HPC_SO = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).$(SO_SUFFIX)
REPAST_HPC_A = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).a
REPAST_HPC_HEADERS = src/repast_hpc/*.h
REPAST_HPC_DEPS = src/repast_hpc/*.d

REPAST_HPC_MAC_NAME=
ifeq ($(USE_MAC),yes)
	REPAST_HPC_MAC_NAME= -install_name "$(INSTALL_PREFIX)/$(REPAST_HPC_SO)"
endif


REPAST_HPC_SRC :=

include src/repast_hpc/module.mk

# OBJECT FILES
REPAST_HPC_OBJECTS=$(patsubst %.cpp,%.o,$(REPAST_HPC_SRC))

-include $(REPAST_HPC_OBJECTS:.o=.d)

.PHONY: all clean
all: repast_hpc

DIRECTORIES = $(INSTALL_PREFIX) $(INSTALL_LIB) $(INSTALL_INCLUDE) $(INSTALL_BIN) \
       $(REPAST_HPC_INSTALL_INCLUDE)

install: install_rhpc

install_rhpc: repast_hpc
	$(Q) "  INSTALL REPAST HPC $(INSTALL_PREFIX)/lib"
	$(E) mkdir -pv $(DIRECTORIES)
	$(E) cp $(CP_ARGS) $(REPAST_HPC_SO) $(INSTALL_LIB)
	$(E) cp $(CP_ARGS) $(REPAST_HPC_A) $(INSTALL_LIB)
	$(E) cp $(CP_ARGS) $(REPAST_HPC_HEADERS) $(REPAST_HPC_INSTALL_INCLUDE)

repast_hpc: $(REPAST_HPC_SO) $(REPAST_HPC_A)

$(REPAST_HPC_SO): $(REPAST_HPC_OBJECTS)
	@-mkdir -p lib
	$(CXXLD) $(REPAST_HPC_MAC_NAME) $(LIB_LD_FLAGS) $(REPAST_HPC_OBJECTS) -o $(REPAST_HPC_SO) $(LIBS)  $(RPATHS)

$(REPAST_HPC_A): $(REPAST_HPC_OBJECTS)
	@-mkdir -p lib
	$(Q) "  AR		$(@)"
	$(E) $(AR) $(ARFLAGS) $(@) $(REPAST_HPC_OBJECTS)
	$(E) $(RANLIB) $(@)


%.o : %.cpp
	$(Q) "  CC		$(@)"
	$(E) $(CXX) $(CXXFLAGS) -Wno-deprecated-declarations $(LIB_CXX_FLAGS) -c $< -o $@

clean::
	@echo
	@echo "  CLEAN"
	@rm -fv $(REPAST_HPC_OBJECTS)
	@rm -fv $(REPAST_HPC_SO)
	@rm -fv $(REPAST_HPC_DEPS)
