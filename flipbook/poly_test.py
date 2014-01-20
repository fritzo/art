
from pyx import *

p1 = path.path(path.moveto(0, 0),
               path.lineto(2, 0),
               path.lineto(2, 2),
               path.lineto(0, 2),
               path.closepath())

p2 = path.path(path.moveto(1, 1),
               path.lineto(3, 1),
               path.lineto(3, 3),
               path.lineto(1, 3),
               path.closepath())

c = canvas.canvas()
c.fill(p1, [color.rgb.blue])
c.fill(p2, [color.rgb.black])
c.writeEPSfile("poly_test")

