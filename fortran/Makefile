FC=gfortran

build/reserve_boot.so: src/reserve_boot.f90 build/random.o
	$(FC) -fPIC -shared $^ -o $@ -I./lib

build/random.o: src/random.f90
	$(FC) -c $^ -o $@ -J./lib 
	
.PHONY: test
test: test/test.f90 build/reserve_boot.so
	$(FC) -Wl,-rpath=./build $^ -o test/test

.PHONY: clean
clean:
	rm -rf build/* lib/*