default: ODEs_in_R.pdf ODEs_in_R.R

ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

%.pdf: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"pdf_document\")"

%.html: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.md\",output_format=\"html_document\")"

%.R: %.Rmd
	Rscript --vanilla -e "knitr::purl(\"$*.Rmd\",output=\"$*.R\",documentation=2)"

%.R: %.Rnw
	Rscript --vanilla -e "knitr::purl(\"$*.Rnw\",output=\"$*.R\",documentation=2)"

%.pdf: %.Rnw
	Rscript --vanilla -e "knitr::knit(\"$*.Rnw\",output=\"$*.pdf\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure


