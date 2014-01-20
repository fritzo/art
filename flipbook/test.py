from pyx import *

#draw original image
c1 = canvas.canvas()
c1.text(0, 0, "Hello, world!")
c1.stroke(path.line(0, 0, 2, 0))

#calculate rescaling
paperSize = [10.75, 11.75]
b = c1.bbox()
canvasSize = [unit.toinch(b.height()), unit.toinch(b.width())]
s = min([paperSize[i]/canvasSize[i] for i in range(2)])
t = trafo.scale(s)

#create rescaled canvas
c2 = canvas.canvas()
c2.insert(c1, [t])

c2.writetofile("hello")
c2.writeEPSfile("hello")

