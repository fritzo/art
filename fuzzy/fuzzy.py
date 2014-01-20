
from math import sqrt, pow
from pyx import *

def drawPoint(canv, x, y, mass):
  circ = path.circle(x, y, sqrt(mass))
  canv.fill(circ, [color.rgb.black])
  
def run(thresh=1e-6):
  R = 10.0
  canv = canvas.canvas()
  points = open("fuzzy.text")
  point = points.readline().split()
  while point:
    mass = float(point[2])
    if mass > thresh:
      x = float(point[0])
      y = float(point[1])
      radius = sqrt(x*x + y*y)
      if radius > 0:
        factor = R*4/sqrt(radius)
        x *= factor
        y *= factor
      mass = sqrt(mass)
      drawPoint(canv, x, y, mass)
    point = points.readline().split()
  points.close()
  canv.writetofile("fuzzy.eps")

  #draw border
  bord = path.path(path.moveto(R,R),
                   path.lineto(R,-R),
                   path.lineto(-R,-R),
                   path.lineto(-R,R),
                   path.closepath())
  canv.stroke(bord)

if __name__ == "__main__":
  run()
