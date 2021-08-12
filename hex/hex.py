
from math import *
from pyx import *
import sys

def draw_polygon (ps,c, canv):
  outline = (
      [path.moveto(*(ps[0]))]
    + [path.lineto(*p) for p in ps[1:]]
    + [path.closepath()]
  )
  canv.fill(path.path(*outline),[c])

def draw_hexagon (x,y,r,c, canv):
  p = [(x+r*cos(pi*t/3.0), y+r*sin(pi*t/3.0)) for t in range(6)]
  draw_polygon(p,c, canv)

bit = 0.1
more = 1.0+bit
less = 1.0-bit

def draw_space(x,y, canv):
  r = 1.0/sqrt(3.0)
  ro = more * r
  ri = less * r
  draw_hexagon(x,y,ro,color.cmyk.Black, canv)
  draw_hexagon(x,y,ri,color.cmyk.White, canv)

def draw_board(size, canv):
  #axes
  u = (0.0,1.0)
  v = (sqrt(3.0)/2.0, 0.5)

  #draw orienters
  s = size + 1
  U = (s * u[0], s * u[1])
  V = (s * v[0], s * v[1])
  board = [
      (v[0],v[1]),
      (u[0],u[1]),
      (U[0],U[1]),
      (U[0]+V[0]-v[0],U[1]+V[1]-v[1]),
      (U[0]+V[0]-u[0],U[1]+V[1]-u[1]),
      (V[0],V[1]),
  ]
  draw_polygon(board, color.cmyk.Black, canv)
  #shrink board
  for i in range(6):
    x,y = board[i]
    t = pi*i/3.0
    x += bit*sin(t)
    y += bit*cos(t)
    board[i] = x,y
  #draw_polygon(board, color.cmyk.White, canv)
  def mid(ab, cd):
      a, b = ab
      c, d = cd
      return ((a+c)/2.0,(b+d)/2.0)
  draw_polygon(
      [
        mid(board[0],board[1]),
        board[1],
        board[2],
        mid(board[2],board[5])
      ],
      color.cmyk.White, canv
  )
  draw_polygon(
      [
        mid(board[3],board[4]),
        board[4],
        board[5],
        mid(board[2],board[5])
      ],
      color.cmyk.White, canv
  )

  #draw spaces
  x0 = u[0]+v[0]
  y0 = u[1]+v[1]
  for i in range(size):
    for j in range(size):
      x = x0 + i*u[0] + j*v[0]
      y = y0 + i*u[1] + j*v[1]
      draw_space(x,y, canv)

def rotate (old_c, degs):
  t = trafo.rotate(degs)
  new_c = canvas.canvas()
  new_c.insert(old_c, [t])
  return new_c

def resize (old_c, new_width=10.5, new_height=16.5):
  b = old_c.bbox()
  old_height = unit.toinch(b.height())
  old_width  = unit.toinch(b.width())
  s = min(new_height/old_height, new_width/old_width)
  t = trafo.scale(s)
  new_c = canvas.canvas()
  new_c.insert(old_c, [t])
  return new_c

def recenter (old_c, border=0.25):
  b = old_c.bbox()
  dx = border-unit.toinch(b.left())
  dy = border-unit.toinch(b.bottom())
  t = trafo.translate(dx,dy)
  new_c = canvas.canvas()
  new_c.insert(old_c, [t])
  return new_c

def save_board(size):
  assert 0 < size <= 100, "size out of range: %s" % size
  c = canvas.canvas()
  draw_board(size,c)
  c = recenter(resize(rotate(c,-10.0),8.0,10.5)) #letter
  #c = recenter(resize(c,10.5,16.5)) #11x17
  c.writeEPSfile("hex-%ix%i.eps" % (size,size))

if __name__ == "__main__":
  args = sys.argv
  if len(args) < 2:
    print("usage: python hex.py SIZE")
  else:
    size = int(args[1])
    save_board(size)

