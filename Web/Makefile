# Makefile for Sphinx documentation
#
#
############################################################
# Generate and copy the include files for the documetation #
############################################################
#
INCOMP_INCLUDES_SOURCE:=$(realpath $(CURDIR)/../dev/incompressible_liquids )
#INCOMP_INCLUDES_TARGET:=$(realpath $(CURDIR)/fluid_properties )
INCOMP_INCLUDES_TARGET:=fluid_properties
INCOMP_INCLUDES_LOGFIL:=$(INCOMP_INCLUDES_TARGET)/generator.log.txt
# JSON, PDF and JPG all depend on Python, regenerate all of them if Python sources change
INCOMP_INCLUDES_PYTHON:=$(shell find $(INCOMP_INCLUDES_SOURCE)/CPIncomp -iname '*.py') # recursively
INCOMP_INCLUDES_PYTHON:=$(INCOMP_INCLUDES_PYTHON) $(INCOMP_INCLUDES_SOURCE)/all_incompressibles.py
INCOMP_INCLUDES_PYTHON:=$(INCOMP_INCLUDES_PYTHON) $(INCOMP_INCLUDES_SOURCE)/Incompressibles_include_generator.py
#INCOMP_INCLUDES_JSONLI:=$(shell find $(INCOMP_INCLUDES_SOURCE)/json -iname '*.json') # recursively
#INCOMP_INCLUDES_FNAMES:=$(basename $(notdir $(INCOMP_INCLUDES_JSONLI) ) ) # all fluid IDs
## All the required files
#INCOMP_INCLUDES_PDFREP:=$(addprefix $(INCOMP_INCLUDES_TARGET)/report/, $(addsuffix _fitreport.pdf, $(INCOMP_INCLUDES_FNAMES) ) )
INCOMP_INCLUDES_PDFEXT:= all_incompressibles.pdf
INCOMP_INCLUDES_PDFEXT:=$(addprefix $(INCOMP_INCLUDES_TARGET)/, $(INCOMP_INCLUDES_PDFEXT) )
INCOMP_INCLUDES_TXTTAB:= pure-fluids.txt mass-based-fluids.txt mole-based-fluids.txt volume-based-fluids.txt
INCOMP_INCLUDES_TXTTAB:=$(addprefix $(INCOMP_INCLUDES_TARGET)/, $(INCOMP_INCLUDES_TXTTAB) )
# Special figures
INCOMP_INCLUDES_COMFIG:=$(INCOMP_INCLUDES_TARGET)/report2up
INCOMP_INCLUDES_COMFIG_SRC:=$(INCOMP_INCLUDES_TARGET)/report/NBS_fitreport.pdf $(INCOMP_INCLUDES_TARGET)/report/AEG_fitreport.pdf


# the script generates the overview PDF (and more)
$(INCOMP_INCLUDES_PDFEXT): run_incomp_script

# Make sure this expensive step only gets carried out once
INCOMP_INCLUDES_DUMMY:=$(shell date; cd $(INCOMP_INCLUDES_SOURCE); python all_incompressibles.py -ns -fr )
.PHONY : run_incomp_script
run_incomp_script: $(INCOMP_INCLUDES_PYTHON)
	@echo "\nOutput from fitting the incompressibles:\n\
	  $(INCOMP_INCLUDES_DUMMY)\n" > $(INCOMP_INCLUDES_LOGFIL)
	@rsync -am --include='*.pdf' --include='*/' --exclude='*' --delete $(INCOMP_INCLUDES_SOURCE)/ $(INCOMP_INCLUDES_TARGET)/

# Note that we have to use a static pattern to avoid multiple invocations of the generator!
$(INCOMP_INCLUDES_TXTTAB): %.txt: $(INCOMP_INCLUDES_SOURCE)/Incompressibles_include_generator.py $(INCOMP_INCLUDES_PYTHON)
	@cd $(<D); python Incompressibles_include_generator.py

# Create the image that is used to describe the fitting reports
$(INCOMP_INCLUDES_COMFIG).jpg: $(INCOMP_INCLUDES_COMFIG).pdf
	@convert -density 300 $< $@ # Convert the PDF to JPG
	@convert -crop 100%x47%+0+30 -resize '850x' -quality 75 $@ $@ # Resize it 

$(INCOMP_INCLUDES_COMFIG).pdf: $(INCOMP_INCLUDES_COMFIG_SRC)
	@pdfnup --quiet --nup 2x1 --delta "1cm 0cm" $^ -o $(INCOMP_INCLUDES_COMFIG).pdf

# Define one single target to be referenced 
.PHONY: incompressible_includes
incompressible_includes: $(INCOMP_INCLUDES_TXTTAB) $(INCOMP_INCLUDES_PDFEXT) $(INCOMP_INCLUDES_COMFIG).jpg







# You can set these variables from the command line.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
PAPER         =

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d _build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml pickle json htmlhelp qthelp latex changes linkcheck doctest

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  dirhtml   to make HTML files named index.html in directories"
	@echo "  pickle    to make pickle files"
	@echo "  json      to make JSON files"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  qthelp    to make HTML files and a qthelp project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview of all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  doctest   to run all doctests embedded in the documentation (if enabled)"

clean:
	-rm -rf _build/*

html: incompressible_includes
	$(SPHINXBUILD) -j3 -b html $(ALLSPHINXOPTS) _build/html
	@echo
	@echo "Build finished. The HTML pages are in _build/html."

dirhtml:
	$(SPHINXBUILD) -b dirhtml $(ALLSPHINXOPTS) _build/dirhtml
	@echo
	@echo "Build finished. The HTML pages are in _build/dirhtml."

pickle:
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) _build/pickle
	@echo
	@echo "Build finished; now you can process the pickle files."

json:
	$(SPHINXBUILD) -b json $(ALLSPHINXOPTS) _build/json
	@echo
	@echo "Build finished; now you can process the JSON files."

htmlhelp:
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) _build/htmlhelp
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in _build/htmlhelp."

qthelp:
	$(SPHINXBUILD) -b qthelp $(ALLSPHINXOPTS) _build/qthelp
	@echo
	@echo "Build finished; now you can run "qcollectiongenerator" with the" \
	      ".qhcp project file in _build/qthelp, like this:"
	@echo "# qcollectiongenerator _build/qthelp/sampledoc.qhcp"
	@echo "To view the help file:"
	@echo "# assistant -collectionFile _build/qthelp/sampledoc.qhc"

latex:
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) _build/latex
	@echo
	@echo "Build finished; the LaTeX files are in _build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

changes:
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) _build/changes
	@echo
	@echo "The overview file is in _build/changes."

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) _build/linkcheck
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in _build/linkcheck/output.txt."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) _build/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in _build/doctest/output.txt."

sf:
	cd _build/html; \
	rsync -avz . jdh2358,matplotlib@web.sf.net:/home/groups/m/ma/matplotlib/htdocs/sampledoc -essh --cvs-exclude


sfpdf:
	cd _build/latex; \
	scp sampledoc.pdf jdh2358,matplotlib@web.sf.net:/home/groups/m/ma/matplotlib/htdocs/sampledoc/


rst: FluidInfoGenerator.py
	python FluidInfoGenerator.py
	mv Fluids/FluidInformation.rst Fluids/FluidInformation.rst.2
	rm Fluids/*.rst
	mv Fluids/FluidInformation.rst.2 Fluids/FluidInformation.rst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) _build/html




