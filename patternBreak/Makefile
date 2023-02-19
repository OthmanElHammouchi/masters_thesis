R ?= R

.PHONY: all
all:
	$(MAKE) clean
	$(MAKE) build
	$(MAKE) install

.PHONY: build
build:
	$(R) CMD build .

.PHONY: install
install: 
	$(R) CMD INSTALL patternBreak_1.0.tar.gz
	
.PHONY: clean
clean:
	$(RM) patternBreak_1.0.tar.gz
	$(RM) src/*.o
	$(RM) src/*.so