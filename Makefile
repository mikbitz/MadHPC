REPAST_HPC_VERSION = 2.2.0

ENABLE_SHARED = @ENABLE_SHARED@
ENABLE_STATIC = @ENABLE_STATIC@
BASE_DIR = /home/mb425/repastHPC
LDFLAGS :=  -L$(BASE_DIR)/HDF/lib64 -L$(BASE_DIR)/NetCDF/lib64 -L$(BASE_DIR)/NetCDF-cxx/lib64 -L$(BASE_DIR)/CURL/lib -L$(BASE_DIR)/Boost/Boost_1.61/lib -L$(BASE_DIR)/repast_hpc-2.2.0/lib/ -L$(BASE_DIR)/MPICH/lib64
CPPFLAGS := -I$(BASE_DIR)/HDF/include -I$(BASE_DIR)/NetCDF/include -I$(BASE_DIR)/NetCDF-cxx/include -I$(BASE_DIR)/CURL/include -I$(BASE_DIR)/Boost/Boost_1.61/include -I$(BASE_DIR)/repast_hpc-2.2.0/include/ -I$(BASE_DIR)/MPICH/include/

CXX = $(BASE_DIR)/MPICH/bin/mpicxx
CXXLD = $(BASE_DIR)/MPICH/bin/mpicxx
CXXFLAGS = -g -O2 -Wall  -std=c++11 -MMD -MP $(CPPFLAGS) -Wno-reorder -Wno-unused-variable -Wno-sign-compare -I./src

LIB_LD_FLAGS = $(LDFLAGS)
LIBS = -lboost_mpi-mt -lboost_serialization-mt -lboost_filesystem-mt -lboost_system-mt
LIBS += -lnetcdf -lnetcdf_c++4 -lhdf5 -lhdf5_hl -lcurl

BOOST_LIB_DIR = $(BASE_DIR)/Boost/Boost_1.61/lib

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

INSTALL_PREFIX = $(BASE_DIR)/mad_modified
INSTALL_BIN = $(INSTALL_PREFIX)/bin

REPAST_HPC_LIB := repast_hpc
REPAST_HPC_SO = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).$(SO_SUFFIX)
REPAST_HPC_A = lib/lib$(REPAST_HPC_LIB)-$(REPAST_HPC_VERSION).a
REPAST_HPC_HEADERS = src/repast_hpc/*.h
REPAST_HPC_DEPS = src/repast_hpc/*.d

MADMODEL_EXE = mad_model
MADMODEL_LD_FLAGS = $(LDFLAGS) -L./lib
MADMODEL_INSTALL_DIR = $(INSTALL_BIN)
MADMODEL_PROPS = model.props
MADMODEL_CONFIG = config.props
MADMODEL_DEPS = *.d
MADMODEL_LIBS = $(LIBS) -lrepast_hpc-$(REPAST_HPC_VERSION)

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
	$(CXXLD) $(MADMODEL_LD_FLAGS) $(MADMODEL_OBJECTS) -o $(MADMODEL_EXE) $(MADMODEL_LIBS) $(RPATHS)


%.o : %.cpp
	$(Q) "  CC		$(@)"
	$(E) $(CXX) $(CXXFLAGS) -g -Wno-deprecated-declarations $(LIB_CXX_FLAGS) -c $< -o $@

clean::
	@echo
	@echo "  CLEAN"
	@rm -fv $(MADMODEL_OBJECTS)
	@rm -fv $(MADMODEL_EXE)
	@rm -fv $(MADMODEL_DEPS)
