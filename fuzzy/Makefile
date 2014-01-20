view: fuzzy.pdf
	acroread fuzzy.pdf &

fuzzy: fuzzy.C
	#g++ -ggdb fuzzy.C -o fuzzy
	g++ -O3 fuzzy.C -o fuzzy

fuzzy.text: fuzzy
	nice fuzzy

fuzzy.eps: fuzzy.text fuzzy.py
	python fuzzy.py

fuzzy.dvi: fuzzy.eps fuzzy.tex
	latex fuzzy.tex

fuzzy.pdf: fuzzy.dvi
	dvipdf fuzzy.dvi

profile: fuzzy.C
	g++ -pg -ftest-coverage -fprofile-arcs fuzzy.C -o fuzzy
	./fuzzy 1e-3
	gcov fuzzy.C
	gprof -l -b fuzzy

clean:
	rm -f fuzzy *.text *.eps *.aux *.dvi *.pdf *.gcov *.bb *.bbg *.out *.da *.log

