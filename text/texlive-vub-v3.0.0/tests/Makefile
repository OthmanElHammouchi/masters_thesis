departments=etro brubotics soft ai indi
outputname=test.pdf \
		   test_bruface.pdf \
		   $(departments:=-beamer.pdf)

BEAMERSTY=../beamerthemevub.sty \
          ../beamercolorthemevub.sty \
          ../beamerfontthemevub.sty \
          ../beamerinnerthemevub.sty \
          ../beamerouterthemevub.sty \
		  ../vubprivate.sty \

TITLESTY=../vub.sty ../bruface.sty ../bruface.png ../vub_logo_cmyk.pdf

.PHONY: clean all

all: $(outputname)

%.pdf: %.tex $(BEAMERSTY) $(TITLESTY)
	export TEXINPUTS="$(PWD)/..":$$TEXINPUTS; \
	latexmk -pdf -xelatex -latexoption='--shell-escape' -use-make $*

%-beamer.tex: beamer-department.tex $(BEAMERSTY)
	cp beamer-department.tex $*-beamer.tex
	sed -ie s/DEPARTMENT/$*/ $*-beamer.tex


clean:
	latexmk -c -C
	rm _minted-main -rf
