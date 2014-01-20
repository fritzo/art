

from exceptions import Exception
from math import *
from pyx import *

LATER = Exception('unfinished code in mandala.py')

density = 20
barwidth = 0.333

_log_2 = log(2)
log2 = lambda x: log(x)/_log_2
heapHash = lambda x: log2(2*x+1)

c = canvas.canvas()

left = lambda n: 2*n
right = lambda n: 2*n+1

def polar2cart(r,a):
  rho = 16.0*(1.0-2.0**(-r))
  return rho*cos(2*pi*a),rho*sin(2*pi*a)

def pointPos(n):
  t = log2(2*n+1)
  return polar2cart(t,t)

pointSize = lambda n: 0.5*(1.0-barwidth)*(1.0-exp(-4.0*sqrt(2.0)*heapHash(n)/n))


#==========[ drawing the Nodes ]==========
def drawNodes(depth):
  N = 2**depth
  for n in range(1,N+1):
    x,y = pointPos(n)
    r = pointSize(n)
    c.stroke(path.circle(x,y,r),canvas.linewidth(r/6.0))


#==========[ drawing the archimedian spiral ]==========
def drawSpiral(depth):
 #define times
  secants = 6*density
  times = [(1.0*n)/secants for n in range(secants*(depth+1)+1)]
 #draw points
  p = path.path(path.moveto(0,0))
  for t in times:
    x,y = polar2cart(t,t)
    p.append(path.lineto(x,y))
  c.stroke(p,canvas.linewidth(barwidth),color.rgb(0.9,0.9,0.9))


#==========[ drawing the tree ]==========
sinkLeft = lambda t0,t: log2(2*t0*2**t+1)
sinkRight = lambda t0,t: log2(2*((t0+1)*2**t-1)+1)

def drawLeftBranch(n,depth):
 #define times
  times = [0.0]
  multiplier = 2.0**(0.4/density)
  while times[-1] < depth:
    times.append(min(depth,(n+times[-1])*multiplier-n))
  for t in range(1,depth):
    times.append(float(t))
  times.sort()
 #draw points
  t = times[0]
  s = sinkLeft(n,t)
  x,y = polar2cart(s,s-t)
  p = path.path(path.moveto(x,y))
  for t in times:
    s = sinkLeft(n,t)
    x,y = polar2cart(s,s-t)
    p.append(path.lineto(x,y))
  c.stroke(p,[color.rgb.blue])

def drawLeftBranch_alt(n,depth):
  x,y = pointPos(n)
  p = path.path(path.moveto(x,y))
  for t in range(depth):
    n = left(n)
    x,y = pointPos(n)
    p.append(path.lineto(x,y))
  c.stroke(p,[color.rgb.blue])

def drawRightBranch(n,depth):
 #define times
  times = [0.0]
  multiplier = 2.0**(0.4/density)
  while times[-1] < depth:
    times.append(min(depth,(n+times[-1])*multiplier-n))
  for t in range(1,depth):
    times.append(float(t))
  times.sort()
 #draw points
  t = times[0]
  s = sinkRight(n,t)
  x,y = polar2cart(s,s-t)
  p = path.path(path.moveto(x,y))
  for t in times:
    s = sinkRight(n,t)
    x,y = polar2cart(s,s-t)
    p.append(path.lineto(x,y))
  c.stroke(p,[color.rgb.red])

def drawRightBranch_alt(n,depth):
  x,y = pointPos(n)
  p = path.path(path.moveto(x,y))
  for t in range(depth):
    n = right(n)
    x,y = pointPos(n)
    p.append(path.lineto(x,y))
  c.stroke(p,[color.rgb.red])

def fillTree(n,depth):
  if depth > 1:
    drawRightBranch(left(n),depth-1)
    fillTree(left(n),depth-1)
    drawLeftBranch(right(n),depth-1)
    fillTree(right(n),depth-1)

def drawTree(n,depth):
  drawLeftBranch(n,depth)
  drawRightBranch(n,depth)
  fillTree(n,depth)



#==========[ drawing the hook ]==========
def drawHook():
 #define times
  N = 10
  secants = 6*density
  times = [-float(n)/secants for n in range(N*secants)]
 #draw points
  t = times[0]
  s = sinkLeft(1,t)
  x,y = polar2cart(s,s-t)
  p = path.path(path.moveto(x,y))
  for t in times:
    s = sinkLeft(1,t)
    x,y = polar2cart(s,s-t)
    p.append(path.lineto(x,y))
  c.stroke(p,[color.rgb.blue])


#==========[ drawing the mandala ]==========
def drawMandala(depth):
  #drawSpiral(depth)
  drawTree(1,depth)
  #drawHook()
  #drawNodes(depth)
  #c.stroke(path.circle(0,0,depth+3))
  c.writetofile('mandala.pdf')
              

if __name__ == '__main__':
  drawMandala(10)
