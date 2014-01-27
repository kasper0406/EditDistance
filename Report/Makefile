thesis:
	latex thesis.tex
	bibtex A
	bibtex B
	latex thesis.tex
	dvips -t letter -o thesis.tmp thesis.dvi
	mv -f thesis.tmp thesis.ps
	ps2pdf14 -dPDFSETTINGS=/prepress -dEmbedAllFonts=true thesis.ps thesis.tmp
	mv -f thesis.tmp thesis.pdf

latex:
	latex thesis
	bibtex A
	bibtex B
	latex thesis && (pkill -USR1 xdvi\* || echo "No xdvi!")

xdvi:
	xdvi -background white thesis.dvi -s 5 -sidemargin 0.0cm -topmargin 1.0cm -hush -expertmode 0 -geometry 750x987+0+0 >& /dev/null &

