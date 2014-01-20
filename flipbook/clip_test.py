
from pyx import *

diamond = path.normpath(path.path(
  path.moveto(-3,0),
  path.lineto(0,3),
  path.lineto(3,0),
  path.lineto(0,-3),
  path.closepath()
))

square = path.normpath(path.path(
  path.moveto(-2,-2),
  path.lineto(2,-2),
  path.lineto(2,2),
  path.lineto(-2,2),
  path.closepath()
))

#original canvas
c1 = canvas.canvas()
c1.stroke(diamond)

#clipped canvas
clipper = canvas.clip(square)
c2 = canvas.canvas()
c2.insert(c1, [clipper])

c2.writeEPSfile("clip_test.eps")
