.PHONY: all clean simulations plots shiny test example example-plots example-test

TEXTDIR := text
DOC := $(TEXTDIR)/thesis.tex

RERUN := "(undefined references|Rerun to get (cross-references|the bars|point totals) right|Table widths have changed. Rerun LaTeX.|Linenumber reference failed)"
RERUNBIB := "No file.*\.bbl|Citation.*undefined"

PDFTEXOPTS := -shell-escape -output-directory $(TEXTDIR)
BIBEROPTS := --output-directory $(TEXTDIR)
MAKEINDEXOPTS := $(TEXTDIR)/thesis.nlo -s $(TEXTDIR)/nomencl.ist -o $(TEXTDIR)/thesis.nls

all: clean example-plots doc

doc: clean $(DOC:.tex=.pdf)

$(TEXTDIR)/%.pdf: $(TEXTDIR)/%.tex
	pdflatex $(PDFTEXOPTS) $<
	@egrep -q $(RERUNBIB) $(TEXTDIR)/$*.log && biber $(BIBEROPTS) $* && pdflatex $(PDFTEXOPTS) $<; true
	makeindex $(MAKEINDEXOPTS) && pdflatex $(PDFTEXOPTS) $<; true
	@egrep -q $(RERUN) $(TEXTDIR)/$*.log && pdflatex $(PDFTEXOPTS) $<; true
	@egrep -q $(RERUN) $(TEXTDIR)/$*.log && pdflatex $(PDFTEXOPTS) $<; true

clean:
	-rm -f $(DOC:.tex=.pdf)
	-rm -f $(TEXTDIR)/*.{aux,dvi,log,bbl,blg,brf,fls,toc,thm,out,fdb_latexmk}

EXAMPLEDIR = $(wildcard scripts/example/*.r)

example:
	for file in $(EXAMPLEDIR); do\
			$${file} -r --nboot 1000 --nsim 100;\
	done

example-plots:
	-rm -f plots/example/*
	for file in $(EXAMPLEDIR); do\
			$${file};\
	done

example-test:
	for file in $(EXAMPLEDIR); do\
			$${file} -r --nboot 10 --nsim 10;\
	done

simulations:
	Rscript scripts/simulations.r

plots:
	Rscript scripts/plots.r

shiny:
	R -e 'shiny::runApp("scripts/shiny.r")'

test:
	R -e 'source("patternBreak/tests/run_tests.r")'