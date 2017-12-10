all: mean-and-sfm.eps numerical-verification.eps

%.eps: %.pdf
	/usr/bin/pdftops -eps $(<F) $%

mean-and-sfm.pdf: mean-and-sfm-plot.R moment-equations.R plotting.R
	./$< >$(<F).so 2>$(<F).se

numerical-verification.pdf: numerical-verification-plot.R parsp.rds plotting.R
	./$< >$(<F).so 2>$(<F).se

parsp.rds: numerical-verification-moment-calcs.R flag-ts-data
	./$< >$(<F).so 2>$(<F).se

flag-ts-data: simulate-time-series.R
	./$< >$(<F).so 2>$(<F).se && touch flag-ts-data
