SHELL := /bin/sh

APP := $(abspath ..)
MAIN := $(APP)/cpp/Main_program
TOOLS := $(APP)/cpp/OpenMX_tools_work1
BUILD := $(APP)/.build
BIN := $(APP)/bin

UNAME_S := $(shell uname -s)
PKG_CONFIG ?= pkg-config
HDF5_PREFIX ?=
LIBOMP_PREFIX ?=

ifneq ($(strip $(HDF5_PREFIX)),)
HDF5_CFLAGS ?= -I$(HDF5_PREFIX)/include
HDF5_LIBS ?= $(HDF5_PREFIX)/lib/libhdf5_hl_cpp.a \
 $(HDF5_PREFIX)/lib/libhdf5_cpp.a $(HDF5_PREFIX)/lib/libhdf5_hl.a \
 $(HDF5_PREFIX)/lib/libhdf5.a -lz -ldl -lm
else
HDF5_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags hdf5 2>/dev/null)
HDF5_LIBS ?= -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl \
 $(shell $(PKG_CONFIG) --libs hdf5 2>/dev/null)
endif

ifeq ($(UNAME_S),Darwin)
OPENMP_CFLAGS ?= -Xpreprocessor -fopenmp
ifneq ($(strip $(LIBOMP_PREFIX)),)
OPENMP_CPPFLAGS ?= -I$(LIBOMP_PREFIX)/include/libomp
OPENMP_LIBS ?= -L$(LIBOMP_PREFIX)/lib/libomp \
 -Wl,-rpath,$(LIBOMP_PREFIX)/lib/libomp -lomp
else
OPENMP_CPPFLAGS ?=
OPENMP_LIBS ?= -lomp
endif
MATH_LIBS ?= -framework Accelerate
STACK_SIZE_FLAG ?= -Wl,-stack_size,0x10000000
else
OPENMP_CFLAGS ?= -fopenmp
OPENMP_CPPFLAGS ?=
OPENMP_LIBS ?= -fopenmp
MATH_LIBS ?= -lblas -llapack
STACK_SIZE_FLAG ?=
endif

CPPFLAGS := -DOLD_HEADER_FILENAME -DHDF_NO_NAMESPACE -DNO_STATIC_CAST \
 $(HDF5_CFLAGS) $(OPENMP_CPPFLAGS)
CXXFLAGS ?= -O3 -std=gnu++17
CXXFLAGS += $(OPENMP_CFLAGS)
LIBS := $(OPENMP_LIBS) $(HDF5_LIBS) $(MATH_LIBS)

LEGACY_SOURCES := $(filter-out $(MAIN)/FPFS_analysis_main.cpp \
 $(MAIN)/spadexp_standard_axes.cpp,$(wildcard $(MAIN)/*.cpp))
LEGACY_OBJECTS := $(patsubst $(MAIN)/%.cpp,$(BUILD)/main/%.o,$(LEGACY_SOURCES))

STANDARD_SOURCES := $(filter-out $(MAIN)/FPFS_analysis_main.cpp \
 $(MAIN)/calculation_pad.cpp $(MAIN)/spadexp_standard_axes.cpp,$(wildcard $(MAIN)/*.cpp))
STANDARD_OBJECTS := $(patsubst $(MAIN)/%.cpp,$(BUILD)/main_standard/%.o,$(STANDARD_SOURCES)) \
 $(BUILD)/main_standard/spadexp_standard_axes.o

FOURD_SOURCES := $(filter-out $(MAIN)/FPFS_analysis_main.cpp \
 $(MAIN)/calculation_pad.cpp $(MAIN)/spadexp_standard_axes.cpp,$(wildcard $(MAIN)/*.cpp))
FOURD_OBJECTS := $(patsubst $(MAIN)/%.cpp,$(BUILD)/main_4d/%.o,$(FOURD_SOURCES)) \
 $(BUILD)/main_4d/calculation_pad_batch.o

TOOLS_COMMON := $(BUILD)/tools/inputTools.o $(BUILD)/tools/kfold.o
POST_COMMON := $(BUILD)/tools/inputTools.o $(BUILD)/tools/HDF5_tools.o $(BUILD)/tools/kfold.o
POST_4D_OBJECTS := $(BUILD)/tools_4d/postproc_heap.o $(BUILD)/tools_4d/inputTools.o \
 $(BUILD)/tools_4d/HDF5_tools.o $(BUILD)/tools_4d/kfold.o

.PHONY: all check clean

all: $(BIN)/spadexp $(BIN)/spadexp_standard_axes $(BIN)/spadexp_4d_plus \
 $(BIN)/openmx_preproc $(BIN)/openmx_postproc $(BIN)/openmx_postproc_large \
 $(BIN)/openmx_postproc_standard_axes $(BIN)/openmx_postproc_4d_plus \
 $(BIN)/kfold_coefficients_report

$(BIN)/spadexp: $(LEGACY_OBJECTS) | $(BIN)
	$(CXX) $^ $(LIBS) -o $@

$(BIN)/spadexp_standard_axes: $(STANDARD_OBJECTS) | $(BIN)
	$(CXX) $^ $(LIBS) -o $@

$(BIN)/spadexp_4d_plus: $(FOURD_OBJECTS) | $(BIN)
	$(CXX) $^ $(LIBS) -o $@

$(BUILD)/main/%.o: $(MAIN)/%.cpp $(wildcard $(MAIN)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(MAIN) -c $< -o $@

$(BUILD)/main_standard/%.o: $(MAIN)/%.cpp $(wildcard $(MAIN)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(MAIN) -c $< -o $@

$(BUILD)/main_standard/spadexp_standard_axes.o: $(MAIN)/spadexp_standard_axes.cpp \
 $(MAIN)/calculation_pad.cpp $(wildcard $(MAIN)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(MAIN) -c $< -o $@

$(BUILD)/main_4d/%.o: $(MAIN)/%.cpp $(wildcard $(MAIN)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSPADEXP_STANDARD_AXES=1 -I$(MAIN) -c $< -o $@

$(BUILD)/main_4d/calculation_pad_batch.o: calculation_pad_batch.cpp batch_grid.hpp \
 $(wildcard $(MAIN)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSPADEXP_STANDARD_AXES=1 -I$(MAIN) -I. -c $< -o $@

$(BIN)/openmx_preproc: $(BUILD)/tools/preproc_main.o $(TOOLS_COMMON) | $(BIN)
	$(CXX) $^ -o $@

$(BIN)/openmx_postproc_large: $(BUILD)/tools/postproc_main.o $(POST_COMMON) | $(BIN)
	$(CXX) $^ $(HDF5_LIBS) $(STACK_SIZE_FLAG) -o $@

$(BIN)/openmx_postproc: $(BIN)/openmx_postproc_large | $(BIN)
	cp $< $@

$(BIN)/openmx_postproc_standard_axes: $(BUILD)/tools_standard/postproc_standard_axes_main.o \
 $(POST_COMMON) | $(BIN)
	$(CXX) $^ $(HDF5_LIBS) $(STACK_SIZE_FLAG) -o $@

$(BIN)/kfold_coefficients_report: $(BUILD)/tools/kfold_coefficients_main.o \
 $(BUILD)/tools/kfold.o | $(BIN)
	$(CXX) $^ $(HDF5_LIBS) -o $@

$(BUILD)/tools/%.o: $(TOOLS)/%.cpp $(wildcard $(TOOLS)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -O3 -std=gnu++17 -I$(TOOLS) -c $< -o $@

$(BUILD)/tools_standard/postproc_standard_axes_main.o: \
 $(TOOLS)/postproc_standard_axes_main.cpp $(TOOLS)/postproc_main.cpp $(wildcard $(TOOLS)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) -O3 -std=gnu++17 -I$(TOOLS) -c $< -o $@

$(BIN)/openmx_postproc_4d_plus: $(POST_4D_OBJECTS) | $(BIN)
	$(CXX) $^ $(LIBS) -o $@

$(BUILD)/tools_4d/postproc_heap.o: postproc_heap.cpp $(wildcard $(TOOLS)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSPADEXP_STANDARD_AXES=1 -I$(TOOLS) -c $< -o $@

$(BUILD)/tools_4d/%.o: $(TOOLS)/%.cpp $(wildcard $(TOOLS)/*.hpp)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(TOOLS) -c $< -o $@

$(BIN):
	mkdir -p $@

check: all
	@file $(BIN)/spadexp $(BIN)/spadexp_standard_axes $(BIN)/spadexp_4d_plus \
 $(BIN)/openmx_preproc $(BIN)/openmx_postproc $(BIN)/openmx_postproc_large \
 $(BIN)/openmx_postproc_standard_axes $(BIN)/openmx_postproc_4d_plus \
 $(BIN)/kfold_coefficients_report

clean:
	rm -rf $(BUILD)
