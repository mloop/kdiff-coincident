# Create the datasets

datasets/: genesis.R genesis.job check-datasets.job conditions.txt
	date
	qsub genesis.job
	qsub -hold_jid genesis check-datasets.job
	date
