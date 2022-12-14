.suffxes:i
.oneshell:

FC = ifort
SHELL := /bin/bash

MKLROOT = $(shell echo $$MKLROOT)
export MKL_VERBOSE = 1 

ifeq ($(FC), gfortran)

FFLAGS = -fbackslash -g -fcheck=all -cpp -I./lib -J./lib

else

FFLAGS = -assume bscc -assume nobuffered_stdout  -g -par-runtime-control2 -cpp -module ./lib -L$(MKLROOT)/include/ -I${MKLROOT}/include/intel64/ilp64 -I"${MKLROOT}/include"

LDFLAGS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a -lpthread -lm 

endif
	

$(OBJS): env_setup

SRC_DIR = src
BUILD_DIR = build

SRCS90 = $(wildcard $(SRC_DIR)/*.f90)
SRCS95 = $(wildcard $(SRC_DIR)/*.f95)

SRCS = $(SRCS90) 
SRCS += $(SRCS95)

OBJS90 = $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(SRCS90))
OBJS95 = $(patsubst $(SRC_DIR)/%.f95, $(BUILD_DIR)/%.o, $(SRCS95))

OBJS = $(OBJS90)
OBJS +=$(OBJS95)

TARGET = $(BUILD_DIR)/reserve_sim.so

all: $(TARGET)
	
# env:
# 	$(shell source /opt/intel/oneapi/setvars.sh;)
# 	$(shell source /opt/intel/oneapi/mkl/latest/env/vars.sh;)
# 	$(eval MKLROOT = $(shell echo $$MKLROOT))


-include $(BUILD_DIR)/.depend

depend .depend:
	$(shell makedepf90 -b $(BUILD_DIR) $(SRCS) > $(BUILD_DIR)/.depend)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -fPIC -c $< -o $@


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f95
	$(FC) $(FFLAGS) -fPIC -c $< -o $@


$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -fPIC -shared $^ -o $@ $(LDFLAGS)

.PHONY: test clean

test: test/test.f90 $(OBJS)
	$(FC) $(FFLAGS) -Wl,-rpath=./build $^ -o test/test -I./lib $(LDFLAGS)

clean:
	rm -rf build/* lib/* build/.depend