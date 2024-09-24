install:
	R CMD INSTALL .

reinstall:
	R CMD REMOVE tractable
	R CMD INSTALL .

check:
	R CMD check .

test:
	Rscript -e "devtools::test()"

docs:
	Rscript -e "devtools::document()"
	Rscript -e "devtools::build_readme()"

examples:
	Rscript -e "devtools::build_rmd('vignettes/tractable-single-bundle.Rmd')"
	Rscript -e "devtools::build_rmd('vignettes/changing-k.Rmd')"

clean:
	rm vignettes/*.html
	rm vignettes/*.png
