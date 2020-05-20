.PHONY:	main check build doc scp


main:	check

check:	doc
	R -e 'devtools::test()'

build:	check
	R -e 'devtools::build()'

doc:
	R -e 'devtools::document()'
	R -e 'devtools::build_vignettes()'

scp:	build
	scp ../rucdr_0.1.tar.gz ulb-ws:/tmp/

