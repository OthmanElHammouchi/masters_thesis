simulations:
	Rscript scripts/simulations.r

plots:
	Rscript scripts/plots.r

shiny:
	R -e 'shiny::runApp("scripts/shiny.r")'

test:
	R -e 'source("patternBreak/tests/run_tests.r")'