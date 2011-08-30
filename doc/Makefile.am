
doc_DATA = MOAB-UG.doc \
           MOAB-performance.doc \
           seq.pdf 

seq_dia = seq-uml-class.dia \
          seq-data.dia \
          MeshSet.dia
           
EXTRA_DIST  = $(doc_DATA) \
	      seq.tex $(seq_dia)

SUFFIXES = .dia .eps .dvi .pdf .ps .tex

.dia.eps: 
	dia -e $@ -t eps $<

#.dvi.pdf:
#	dvipdf $< $@	

.dvi.ps:
	dvips $< $@	

.tex.dvi:
	latex $<
	latex $<
	latex $<

seq.dvi: seq.tex $(seq_dia:.dia=.eps)

MOSTLYCLEANFILES = *.aux *.bbl *.blg *.lof *.log *.log *.toc *.dvi 
