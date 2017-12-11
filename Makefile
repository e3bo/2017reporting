all: mean-and-sfm.eps numerical-verification.eps example-sims.eps quantiles.eps

%.eps: %.pdf
	/usr/bin/pdftops -eps $(<F) $%

quantiles.pdf: quantile-plot.R flag-mwe-data plotting.R
	./$< >$(<F).so 2>$(<F).se

example-sims.pdf: example-sims-plot.R flag-mwe-data plotting.R
	./$< >$(<F).so 2>$(<F).se

mean-and-sfm.pdf: mean-and-sfm-plot.R moment-equations.R plotting.R
	./$< >$(<F).so 2>$(<F).se

numerical-verification.pdf: numerical-verification-plot.R parsp.rds plotting.R
	./$< >$(<F).so 2>$(<F).se

flag-mwe-data: sim-study-moment-calcs.R flag-sim-data
	./$< >$(<F).so 2>$(<F).se && touch $@

parsp.rds: numerical-verification-moment-calcs.R flag-sim-data
	./$< >$(<F).so 2>$(<F).se

flag-sim-data: simulate-time-series.R
	./$< >$(<F).so 2>$(<F).se && touch $@
