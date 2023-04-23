simulations:
	Rscript scripts/simulations.r

plots:
	Rscript scripts/plots.r

shiny:
	R -e 'shiny::runApp("scripts/shiny.r")'

test:
	R -e 'source("patternBreak/tests/run_tests.r")'

example:
	Rscript scripts/example/diag_plots.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/mack_pairs.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/mack_semiparam.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/mack_param.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/mack_bench.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/boot_est.r -r --nboot 1000 --nsim 1000
	Rscript scripts/example/boot_reg.r -r --nboot 1000 --nsim 1000

example-plots:
	Rscript scripts/example/diag_plots.r
	Rscript scripts/example/mack_pairs.r
	Rscript scripts/example/mack_semiparam.r
	Rscript scripts/example/mack_param.r
	Rscript scripts/example/mack_bench.r
	Rscript scripts/example/boot_est.r
	Rscript scripts/example/boot_reg.r

example-test:
	Rscript scripts/example/diag_plots.r -r --nboot 10 --nsim 10
	Rscript scripts/example/mack_pairs.r -r --nboot 10 --nsim 10
	Rscript scripts/example/mack_semiparam.r -r --nboot 10 --nsim 10
	Rscript scripts/example/mack_param.r -r --nboot 10 --nsim 10
	Rscript scripts/example/mack_bench.r -r --nboot 10 --nsim 10
	Rscript scripts/example/boot_est.r -r --nboot 10 --nsim 10
	Rscript scripts/example/boot_reg.r -r --nboot 10 --nsim 10

TEXFILE=thesis
RDIR= .
FIGDIR= ./figs

$(RDIR)/%.Rout: $(RDIR)/%.R $(RDIR)/functions.R
	R CMD BATCH $< 

$(TEXFILE).pdf: $(TEXFILE).tex $(OUT_FILES) $(CROP_FILES)
	latexmk -pdf -quiet $(TEXFILE)


