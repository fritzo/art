
rubiks.pdf: rubiks.eps
	epstopdf rubiks.eps -o=rubiks.pdf

rubiks.eps: rubiks.py
	python rubiks.py

view: rubiks.pdf
	acroread rubiks.pdf &

clean:
	rm -f *.eps *.pdf *.pyc
