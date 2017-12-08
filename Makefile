all: numerical-verification.pdf

numerical-verification.pdf: numerical-verification.R flag-ts-data
	$< >$(<F).so 2>$(<F).se

flag-ts-data: simulate-time-series.R
	$< >$(<F).so 2>$(<F).se && touch flag-ts-data
