simulations:
	Rscript scripts/simulations.r

plots:
	Rscript scripts/plots.r

shiny:
	R -e 'shiny::runApp("scripts/shiny.r")'