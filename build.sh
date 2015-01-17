DOC=simdjs
LATEX_BATCH=""
LATEX_QUITE=""

pdflatex $LATEX_BATCH -shell-escape $DOC.tex $LATEX_QUIET
bibtex $DOC
pdflatex $LATEX_BATCH -shell-escape $DOC.tex $LATEX_QUIET
pdflatex $LATEX_BATCH -shell-escape $DOC.tex $LATEX_QUIET
