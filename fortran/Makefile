.suffixes:
.ONESHELL:

FC = ifort
SHELL := /bin/bash

MKLROOT = $(shell echo $$MKLROOT)

ifeq ($(FC), gfortran)
FFLAGS = -fbackslash -fcheck=all -cpp -I./lib -J./lib
else
FFLAGS = -assume bscc -assume nobuffered_stdout -mcmodel=large -par-runtime-control2 -fpp -module ./lib -L$(MKLROOT)/include/ -I${MKLROOT}/include/intel64/ilp64 -I"${MKLROOT}/include"

LDFLAGS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group $(MKLROOT)/lib/intel64/libmkl_blas95_ilp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_ilp64.a -lpthread -lm 
endif

SRC_DIR = src
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(SRCS))

TARGET = $(BUILD_DIR)/reserve_sim.so

.phony: all
all: $(TARGET)

include .depend

.phony: depend
depend: .depend

.depend: $(SRCS)
	makedepf90 -b $(BUILD_DIR) $(SRCS) > .depend

$(OBJS): depend

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -fPIC -c $< -o $@

$(TARGET): $(filter-out $(BUILD_DIR)/reserve_boot.o, $(OBJS))
	$(FC) $(FFLAGS) -fPIC -shared $^ -o $@ $(LDFLAGS)

reserve_boot: $(SRC_DIR)/reserve_boot.f90 $(BUILD_DIR)/dispmodule.o
	$(FC) $(FFLAGS) -fPIC -shared $^ -o $(BUILD_DIR)/reserve_boot.so $(LDFLAGS)

.PHONY: test clean debug

test: FFLAGS += -DDEBUG -g -O0

test: test/test.f90 $(filter-out $(BUILD_DIR)/reserve_boot.o, $(OBJS)) 
	$(FC) $(FFLAGS) -g -Wl,-rpath=./build $^ -o test/test -I./lib $(LDFLAGS)

clean:
	rm -rf build/* lib/* .depend

debug: FFLAGS += -DDEBUG -g -O0 -warn -check
debug: .depend $(TARGET)