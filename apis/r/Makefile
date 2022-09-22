rscript := Rscript --no-save --no-restore

rmds     := $(wildcard vignettes/*.Rmd)
mds      := $(rmds:%.Rmd=%.md)
ipynbs   := $(rmds:%.Rmd=%.ipynb)

jupyter: $(ipynbs)

%.ipynb: %.md
	@echo "Converting $< to $@"
	pandoc --from markdown --to ipynb --output $@ $<

%.md: %.Rmd
	@echo "Converting $< to $@"
	@$(rscript) -e "rmarkdown::render(input = '$<', output_file = '$(@F)', output_dir = '$(@D)', output_format = rmarkdown::md_document(variant = 'commonmark', preserve_yaml = TRUE))"
	@sed -i '' 's/``` r/``` code/g' $@

clean:
	rm -f vignettes/*.{md,ipynb}
