FC = gfortran
FCFLAGS = -g -fcheck=all -cpp -fbounds-check

SRC_DIR = src
BUILD_DIR = build

SRCS = $(wildcard $(SRC_DIR)/*.f90)
OBJS = $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(SRCS))

all: $(OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FCFLAGS) -fPIC -c $^ -o $@ -J./lib

$(BUILD_DIR)/reserve_boot.so: $(SRCS) $(OBJS)
	$(FC) $(FCFLAGS) -fPIC -shared $^ -o $@ -I./lib -J./lib

.PHONY: test clean

test: test/test.f90 $(OBJS)
	$(FC) $(FCFLAGS) -Wl,-rpath=./build $^ -o test/test -I./lib

clean:
	rm -rf build/* lib/*