.suffxes:i

FC = ifort

ifeq ($(FC), gfortran)
FFLAGS = -fbackslash -g -fcheck=all -cpp -I./lib -J./lib
else
FFLAGS = -assume bscc -g -par-runtime-control2 -cpp -module ./lib -mkl
endif
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

-include $(BUILD_DIR)/.depend

depend .depend:
	$(shell makedepf90 -b $(BUILD_DIR) $(SRCS) > $(BUILD_DIR)/.depend)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -fPIC -c $< -o $@


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f95
	$(FC) $(FFLAGS) -fPIC -c $< -o $@


$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -fPIC -shared $^ -o $@

.PHONY: test clean

test: test/test.f90 $(OBJS)
	$(FC) $(FFLAGS) -Wl,-rpath=./build $^ -o test/test -I./lib

clean:
	rm -rf build/* lib/* build/.depend