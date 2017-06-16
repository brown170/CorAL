PAPER = paper
FIGS =  

$(PAPER).pdf: $(PAPER).ps
	ps2pdf $(PAPER).ps $(PAPER).pdf

$(PAPER).ps:  $(PAPER).dvi
	dvips -t letter -o $(PAPER).ps $(PAPER).dvi

$(PAPER).dvi: $(PAPER).tex $(FIGS)
	latex $(PAPER); latex $(PAPER)

clean:
	rm -f *.dvi *.aux *.bbl *.log *.toc *.fff *.lof *.fgx *.tbx *.mtc* *.maf *.blg $(PAPER).ps $(PAPER).pdf

