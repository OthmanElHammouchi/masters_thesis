R ?= R

.PHONY: all
all:
	$(MAKE) clean
	$(MAKE) attrs
	$(MAKE) docs
	$(MAKE) build
	$(MAKE) install

.PHONY: attrs
attrs :
	Rscript -e "library(Rcpp) ; compileAttributes(verbose=TRUE)"

.PHONY: docs
docs:
	Rscript -e "library(roxygen2) ; roxygenise()"

.PHONY: build
build:
	$(R) CMD build .

.PHONY: install
install: 
	$(R) CMD INSTALL claimsBoot_1.0.tar.gz
	
.PHONY: clean
clean:
	$(RM) claimsBoot_1.0.tar.gz
	$(RM) src/*.o
	$(RM) src/*.so
	$(RM) src/*.
	$(RM) src/Makevars

.PHONY: check
check:
	_R_CHECK_CRAN_INCOMING_REMOTE_=false $(R) CMD check claimsBoot_1.0.tar.gz --as-cran --ignore-vignettes --no-stop-on-test-error
