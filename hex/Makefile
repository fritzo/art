
all: hex-14x14.pdf hex-19x19.pdf

#P = ps2pdf -sPAPERSIZE="11x17"
P = ps2pdf -sPAPERSIZE="letter"

hex-14x14.eps: hex.py
	python hex.py 14
hex-14x14.pdf: hex-14x14.eps
	$P hex-14x14.eps

hex-19x19.eps: hex.py
	python hex.py 19
hex-19x19.pdf: hex-19x19.eps
	$P hex-19x19.eps

view: hex-19x19.pdf
	kghostview hex-19x19.pdf &

clean:
	rm -f *.eps *.pdf

