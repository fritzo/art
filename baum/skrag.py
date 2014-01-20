

import random
from math import *
from pyx import *

def draw_path (c,pos):
  c.stroke()

def rand (): return random.normalvariate(0,1)

def walk (N=100):
  "nearly constant position"
  x = [[0,0] for i in range(N)]
  for i in range(N-1):
    x[1+i][0] += x[i][0] + rand()
    x[1+i][1] += x[i][1] + rand()
  return x

def run (N=100,t=5):
  "nearly constant velocity"
  x = [[0,0] for i in range(N)]
  v = [[0,0] for i in range(N)]
  e = exp(-1.0/t)
  for i in range(N-1):
    v[1+i][0] += e * v[i][0] + rand() / t
    v[1+i][1] += e * v[i][1] + rand() / t
    x[1+i][0] += x[i][0] + v[i][0]
    x[1+i][1] += x[i][1] + v[i][1]
  return x

def curve (xs):
  p = path.path(path.moveto(*(xs[0])),
                *[path.lineto(*x) for x in xs[1:]])
  p = deformer.smoothed(100).deform(p)
  return p

if __name__ == "__main__":
  c = canvas.canvas()
  c.stroke(curve(run(200)))
  c.writePDFfile("test.pdf")
  
