
all: test.pdf

test.pdf: skrag.py
	python skrag.py

view: test.pdf
	kpdf test.pdf &

clean:
	rm -f *.eps *.pdf

