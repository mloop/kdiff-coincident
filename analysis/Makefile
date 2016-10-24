# Analyze the results

bias-plot.png: bias-plot.R bias.txt
	R CMD BATCH --no-save --no-restore bias-plot.R

type1-error-plot.png: type1-error-plot.R type1-error.txt
	R CMD BATCH --no-save --no-restore type1-error-plot.R
	
power-plot.png: power-plot.R power.txt
	R CMD BATCH --no-save --no-restore power-plot.R

bias.txt: bias.R results.txt
	R CMD BATCH --no-save --no-restore bias.R

type1-error.txt: type1-error.R results.txt
	R CMD BATCH --no-save --no-restore type1-error.R

power.txt: power.R results.txt
	R CMD BATCH --no-save --no-restore power.R

# Perform simulation to generate raw results
results.txt: simulate.R simulate.job stitch.job ../data/datasets/
	date
	qsub simulate.job
	date
	qsub -hold_jid simulate stitch.job
	date
