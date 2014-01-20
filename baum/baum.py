#discrete tree models

import math, random, pyx

#-----[ forced trees ]-----

class Tree:
  def __init__ (self):
    self.keys = [0]
    self.verts = {}
  def branch (self):
    key = random.choice(self.keys)
    id = len(self.keys)
    self.keys.append(id)
    self.verts[id] = key
  def to_branching (N):
    T = math.log(len(self.keys), 2)
    L = N/T
    ts = [math.ceil(random.expovariate(L)) for _ in self.keys]


def tree (V):
  t = Tree()
  for v in range(V):
    t.branch()
  return t

#-----[ poisson trees ]-----

class Branching:
  def __init__ (self):
    self.id = 0
    self.levels = [[(0,0,0)]]
    #self.update()
  def update (self):
    self.prev = dict([((n,j),i) for n,i,j in self.levels])
    self.parent = dict([(j,i) for _,i,j in self.levels if i!=j])
    self.keys = range(len(self.levels)+1)

    #girth
    self.girth = dict([(k,0) for k in self.keys])
    keys = keys[:]
    keys.reverse()
    for k in keys:
      if self.girth[k] == 0: self.girth[k] = 1
      if k != 0: self.girth[self.parent[k]] += self.girth[k]

  def new_id (self):
    self.ids += 1
    return self.ids

  def branch (self, Pmore, Pless):
    assert 0<=Pmore and 0<=Pless and Pmore+Pless<=1, "invalid Pmore,Pless"
    prev = self.levels[-1]
    next = []
    for n,i,j in prev:
      n += 1
      w = random.unif(0.0,1.0)
      if w > Pless:
        continue
      if 1-w <= Pmore:
        next.append((n,j,self.new_id()))
        next.append((n,j,self.new_id()))
      else:
        next.append((n,j,j))
      self.levels.append(next)
    #self.update()

def branching (N,Pmore,Pless):
  "Pmore,Pless are functions : [0,1]->[0,1]"
  b = Branching()
  for n in range(1,N+1):
    n = n/float(N)
    Pm,Pl = Pmore(n),Pless(n)
    b.branch(Pm,Pl)
  b.update()
  return b

#-----[ data on trees ]-----

class Vect:
  def __init__ (self, branch, data=None):
    self.branch = branch
    if data is not None:
      self.data = data
    else:
      self.data = dict([((n,j),0.0) for l in branch.levels for n,j,_ in l])
  
  def __getitem__ (self,key):
    return self.data[key]
  def __setitem__ (self,key,val):
    return self.data[key] = val
  
  def rand (self,sigma=1.0):
    for key in self.data.keys():
      self.data[key] = random.normalvariate(0.0,sigma)
    return self
  
  def __imul__ (self, s):
    for key in self.data.keys():
      self.data[key] *= s
  def __mul__ (self,s):
    v = Vect(self)
    self *= s
    return v

  def __iadd__ (self, v):
    for key in self.data.keys():
      self.data[key] += v[key]
  def __add__ (u,v):
    uv = Vect(u)
    u += v
    return v
  
  def conv (self):
    data = {(0,0):0.0}
    for n in self.branch.keys[1:]:
      raise "LATER"
    return self

#-----[ biology ]-----

def death (N): return lambda t: pow((1.0 - t)/4.0, 1.0/N)
def birth (N): return lambda t: pow(1.0 - t/4.0, 1.0/N)

class Tree:
  def __init__ (N=128):
    self.b = branching(N,death(N),birth(N))
    sigma = math.sqrt(float(N))
    self.x = Vect(b).rand(sigma).conv()
    self.y = Vect(b).rand(sigma).conv()
    self.z = Vect(b).rand(sigma); z[0,0] = 1.0; z.conv()

  def plot (self, canv):
    return "LATER"

#-----[ main ]-----

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

if __name__ == "__main__":
  t = Tree()
  c = canvas.canvas()
  t.plot(c)
  c = recenter(resize(c))
  c.writeEPSfile("tree")

