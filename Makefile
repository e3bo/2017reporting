all: mean-and-sfm.eps numerical-verification.eps example-sims.eps quantiles.eps

%.eps: %.pdf
	/usr/bin/pdftops -eps $(<F) $%

example-sims.pdf: example-sims-plot.R individual-sims.rds plotting.R
	./$< >$(<F).so 2>$(<F).se

mean-and-sfm.pdf: mean-and-sfm-plot.R moment-equations.R plotting.R
	./$< >$(<F).so 2>$(<F).se

numerical-verification.pdf: numerical-verification-plot.R parsp.rds plotting.R
	./$< >$(<F).so 2>$(<F).se

parsp.rds: numerical-verification-moment-calcs.R flag-ts-data
	./$< >$(<F).so 2>$(<F).se

flag-ts-data: simulate-time-series.R
	./$< >$(<F).so 2>$(<F).se && touch flag-ts-data
