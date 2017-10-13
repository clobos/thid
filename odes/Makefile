PDFLATEX = pdflatex
BIBTEX = bibtex
RM = rm -f
CP = cp

default: ODEs_in_R.pdf ODEs_in_R.R Parameter_Estimation.pdf Parameter_Estimation.R

ifneq ($(INSTALLDIR),)
	mkdir -p $(INSTALLDIR)
	install -m0600 $^ $(INSTALLDIR)
endif

%.pdf: %.Rmd
	PATH=/usr/lib/rstudio/bin/pandoc:$$PATH \
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"pdf_document\")"

%.pdf: %.tex
	$(PDFLATEX) $*
	-$(BIBTEX) $*
	$(PDFLATEX) $*
	$(PDFLATEX) $*

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

%.tex: %.Rnw
	Rscript --vanilla -e "knitr::knit(\"$*.Rnw\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak *.bbl *.blg
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure


