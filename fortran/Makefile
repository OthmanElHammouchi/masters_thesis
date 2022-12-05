FC=gfortran
FCFLAGS= -g -fcheck=all

build/reserve_boot.so: src/reserve_boot.f90 build/random.o
	$(FC) $(FCFLAGS) -fPIC -shared $^ -o $@ -I./lib 

build/random.o: src/random.f90
	$(FC) $(FCFLAGS) -c $^ -o $@ -J./lib 

build/ranlib.o: src/ranlib77.f
	$(FC) $(FCFLAGS) -fPIC -std=legacy -c $^ -o $@ -J./lib

.PHONY: test
test: test/test.f90 build/reserve_boot.so
	$(FC) $(FCFLAGS) -Wl,-rpath=./build $^ -o test/test

.PHONY: clean
clean:
	rm -rf build/* lib/*