REPAST_HPC_VERSION = 2.2.0

ENABLE_SHARED = @ENABLE_SHARED@
ENABLE_STATIC = @ENABLE_STATIC@
BASE_DIR = /home/mb425
#LDFLAGS :=  -L$(BASE_DIR)/repastHPC/NetCDF/lib64 -L$(BASE_DIR)/repastHPC/NetCDF-cxx/lib64 -L$(BASE_DIR)/repastHPC/CURL/lib -L$(BASE_DIR)/repastHPC/Boost/Boost_1.61/lib -L$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/lib/ -L$(BASE_DIR)/repastHPC/MPICH/lib64
#CPPFLAGS :=  -I$(BASE_DIR)/repastHPC/NetCDF/include -I$(BASE_DIR)/repastHPC/NetCDF-cxx/include -I$(BASE_DIR)/repastHPC/CURL/include -I$(BASE_DIR)/repastHPC/Boost/Boost_1.61/include -I$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/include/ -I$(BASE_DIR)/repastHPC/MPICH/include/
#LDFLAGS :=  -L/usr/local/netcdfcxx-4.2.1/lib64 -L/usr/local/netcdfc-4.3.3.1/lib64 -L/usr/local/hdf-1.8.15/lib64  -L$(BASE_DIR)/repastHPC/Boost/Boost_1.61/lib -L$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/lib/ -L$(BASE_DIR)/repastHPC/MPICH/lib64
#CPPFLAGS :=  -I/usr/local/netcdfc-4.3.3.1/include -I/usr/local/netcdfcxx-4.2.1/include   -I$(BASE_DIR)/repastHPC/Boost/Boost_1.61/include -I$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/include/ -I$(BASE_DIR)/repastHPC/MPICH/include/
LDFLAGS :=  -L/usr/lib64  -L$(BASE_DIR)/repastHPC/Boost/Boost_1.61/lib -L$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/lib/ -L$(BASE_DIR)/repastHPC/MPICH/lib64
CPPFLAGS :=  -I/usr/include -I$(BASE_DIR)/repastHPC/Boost/Boost_1.61/include -I$(BASE_DIR)/repastHPC/repast_hpc-2.2.0/include/ -I$(BASE_DIR)/repastHPC/MPICH/include/
CXX = $(BASE_DIR)/repastHPC/MPICH/bin/mpicxx
CXXLD = $(BASE_DIR)/repastHPC/MPICH/bin/mpicxx
CXXFLAGS = -g0 -Wall -O2 -std=c++11 -MMD -MP $(CPPFLAGS) -Wno-reorder -Wno-unused-variable -Wno-sign-compare -I./src

LIB_LD_FLAGS = $(LDFLAGS)
LIBS = -lboost_mpi-mt -lboost_serialization-mt -lboost_filesystem-mt -lboost_system-mt
#LIBS += -lnetcdf -lnetcdf_c++4 -lhdf5 -lhdf5_hl -lcurl
LIBS += -lnetcdf -lnetcdf_c++4 -lcurl

#NETCDF_LIB_DIR = $(BASE_DIR)/repastHPC/NetCDF/lib64
#NETCDF_CXX_LIB_DIR = $(BASE_DIR)/repastHPC/NetCDF-cxx/lib64
BOOST_LIB_DIR = $(BASE_DIR)/repastHPC/Boost/Boost_1.61/lib

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

INSTALL_PREFIX = $(BASE_DIR)/repastHPC/mad_modified
INSTALL_BIN = $(INSTALL_PREFIX)/bin

REPAST_HPC_LIB := repast_hpc
REPAST_HPC_SO = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).$(SO_SUFFIX)
REPAST_HPC_A = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).a
REPAST_HPC_HEADERS = src/repast_hpc/*.h
REPAST_HPC_DEPS = src/repast_hpc/*.d

RELOGO_LIB := relogo
RELOGO_SO = lib/lib$(RELOGO_LIB)-$(REPAST_HPC_VERSION).$(SO_SUFFIX)
RELOGO_A = lib/lib$(RELOGO_LIB)-$(REPAST_HPC_VERSION).a
RELOGO_HEADERS = src/relogo/*.h
RELOGO_DEPS = src/relogo/*.d

RELOGO_LD_FLAGS = $(LIB_LD_FLAGS) -L./lib
RELOGO_LIBS = $(LIBS) -lrepast_hpc-$(REPAST_HPC_VERSION)
RELOGO_RPATHS = $(RPATHS) -Wl,-rpath -Wl,$(INSTALL_LIB)


MADMODEL_EXE = mad_model
MADMODEL_LD_FLAGS = $(LDFLAGS) -L./lib
MADMODEL_INSTALL_DIR = $(INSTALL_BIN)
MADMODEL_PROPS = model.props
MADMODEL_CONFIG = config.props
MADMODEL_DEPS = *.d
MADMODEL_LIBS = $(RELOGO_LIBS) -lrelogo-$(REPAST_HPC_VERSION)

DIR := .

MADMODEL_SRC += $(DIR)/main.cpp
MADMODEL_SRC += $(DIR)/FunctionalGroupDefinitions.cpp
MADMODEL_SRC += $(DIR)/Groups.cpp
MADMODEL_SRC += $(DIR)/Stock.cpp
MADMODEL_SRC += $(DIR)/Cohort.cpp
MADMODEL_SRC += $(DIR)/CohortMerger.cpp
MADMODEL_SRC += $(DIR)/CohortPair.cpp
MADMODEL_SRC += $(DIR)/model.cpp
MADMODEL_SRC += $(DIR)/CohortSum.cpp
MADMODEL_SRC += $(DIR)/Environment.cpp
MADMODEL_SRC += $(DIR)/AutotrophProcessor.cpp
MADMODEL_SRC += $(DIR)/TerrestrialCarbon.cpp
MADMODEL_SRC += $(DIR)/HANPP.cpp
MADMODEL_SRC += $(DIR)/Parameters.cpp
MADMODEL_SRC += $(DIR)/Convertor.cpp
MADMODEL_SRC += $(DIR)/Processor.cpp
MADMODEL_SRC += $(DIR)/DataCoords.cpp
MADMODEL_SRC += $(DIR)/DataIndices.cpp
MADMODEL_SRC += $(DIR)/Variable.cpp
MADMODEL_SRC += $(DIR)/Maths.cpp
MADMODEL_SRC += $(DIR)/FileReader.cpp
MADMODEL_SRC += $(DIR)/ClimateVariablesCalculator.cpp
MADMODEL_SRC += $(DIR)/DataRecorder.cpp
MADMODEL_SRC += $(DIR)/InputData.cpp
MADMODEL_SRC += $(DIR)/InputDatum.cpp
MADMODEL_SRC += $(DIR)/GridDatum.cpp
MADMODEL_SRC += $(DIR)/BasicDatum.cpp
MADMODEL_SRC += $(DIR)/DataLayerSet.cpp
MADMODEL_SRC += $(DIR)/DataLayerProcessor.cpp
MADMODEL_SRC += $(DIR)/TimeStep.cpp
MADMODEL_SRC += $(DIR)/DataLayer.cpp
MADMODEL_SRC += $(DIR)/DataLayer2D.cpp
MADMODEL_SRC += $(DIR)/DataLayer3D.cpp
MADMODEL_SRC += $(DIR)/DataLayer2DwithTime.cpp
MADMODEL_SRC += $(DIR)/DataLayer3DwithTime.cpp
MADMODEL_SRC += $(DIR)/UtilityFunctions.cpp
MADMODEL_SRC += $(DIR)/NonStaticSimpleRNG.cpp
MADMODEL_SRC += $(DIR)/randomizer.cpp
MADMODEL_SRC += $(DIR)/repastRandom.cpp

# OBJECT FILES
MADMODEL_OBJECTS=$(patsubst %.cpp,%.o,$(MADMODEL_SRC))

-include $(MADMODEL_OBJECTS:.o=.d)

.PHONY: all clean
all: mad_model

DIRECTORIES = $(MADMODEL_INSTALL_DIR)


install_mad: mad_model
	$(Q) "  INSTALL MADMODEL MODEL $(INSTALL_PREFIX)/bin"
	$(E) mkdir -pv $(DIRECTORIES)
	$(E) cp $(CP_ARGS) $(MADMODEL_EXE) $(MADMODEL_INSTALL_DIR)
	$(E) cp $(CP_ARGS) $(MADMODEL_PROPS) $(MADMODEL_INSTALL_DIR)
	$(E) cp $(CP_ARGS) $(MADMODEL_CONFIG) $(MADMODEL_INSTALL_DIR)

mad_model: $(MADMODEL_OBJECTS)
	@-mkdir -p bin
	$(CXXLD) $(MADMODEL_LD_FLAGS) $(MADMODEL_OBJECTS) -o $(MADMODEL_EXE) $(MADMODEL_LIBS) $(RELOGO_RPATHS)


%.o : %.cpp
	$(Q) "  CC		$(@)"
	$(E) $(CXX) $(CXXFLAGS) -g -Wno-deprecated-declarations $(LIB_CXX_FLAGS) -c $< -o $@

clean::
	@echo
	@echo "  CLEAN"
	@rm -fv $(MADMODEL_OBJECTS)
	@rm -fv $(MADMODEL_EXE)
	@rm -fv $(MADMODEL_DEPS)
