
import pdb

Sum = sum #before Numeric redefines it

from Numeric import *
from math import *
import random
from pyx import *
from exceptions import *
from copy import copy, deepcopy
from time import time

random.seed(int(time()))

class EmptyClass: pass

Marray = (lambda x: 1.0*array(x))
norm = (lambda x: sqrt(innerproduct(x,x)))
normalize = (lambda x: x / sqrt(innerproduct(x,x)))
crossprod = (lambda x,y: Marray([x[1]*y[2]-x[2]*y[1],
                                 x[2]*y[0]-x[0]*y[2],
                                 x[0]*y[1]-x[1]*y[0]]))


#==========[ geometric polygons ]==========
class Poly3D:
  def __init__(self, points, normal, center, color=0):
    self.color = color
    self.points = points
    self.normal = normal
    self.center = center

  def translate(self, shift):
    for point in self.points:
      point += shift
    self.center += shift

  def rotate(self, rotator):
    self.points = [rotatePos(rotator, point) for point in self.points]
    self.normal = rotatePos(rotator, self.normal)
    self.center = rotatePos(rotator, self.center)
    return self

  def transformed(self, shift, rot=None):
    points = [matrixmultiply(rot,p)+shift for p in self.points]
    normal = matrixmultiply(rot,self.normal)
    center = matrixmultiply(rot,self.center)
    return Poly3D(points, normal, center, self.color)

  def rotated(self, rotator):
    points = [matrixmultiply(rotator,p) for p in self.points]
    normal = matrixmultiply(rotator,self.normal)
    center = matrixmultiply(rotator,self.center)
    return Poly3D(points, normal, center, self.color)

  def projectToPlane(self, proj):
    if innerproduct(proj.Z, self.normal) < 0:
      points = [proj.projectPoint(point) for point in self.points]
      depth = innerproduct(proj.Z, self.center-proj.pos)
      return Poly2D(points, depth, self.color)
    else:
      return None
 
class Poly2D:
  def __init__(self, points, depth, color):
    self.points = deepcopy(points)
    self.depth = depth
    self.color = color

  def draw(self, c):
    outline = (
      [path.moveto(*self.points[0].tolist())]
      + [path.lineto(*pt.tolist()) for pt in self.points[1:]]
      + [path.closepath()]
    )
    p = path.path(*outline)
    c.fill(p, [colorTable[self.color]])

  def transform(self, shift, scale=1.0):
    for p in self.points:
      p += shift
      p *= scale

class Path2D:
  def __init__(self, points, color=0):
    self.points = deepcopy(points)
    self.color = color

  def transform(self, shift, scale=1.0):
    for p in self.points:
      p += shift
      p *= scale

  def draw(self, c):
    N = len(self.points)
    if N > 0:
      p = [path.moveto(*self.points[0])]
      for n in range(1,N):
        p.append(path.lineto(*self.points[n]))
      c.stroke(path.path(*p), [colorTable[self.color], style.linewidth.thin])


#==========[ abstract cube ]==========
count = lambda i,j,k: abs(i) + abs(j) + abs(k)
dist = lambda (i1,j1,k1),(i2,j2,k2): abs(i1-i2) + abs(j1-j2) + abs(k1-k2)
Cube = EmptyClass()
PM = [+1,-1]
PZM = [+1,0,-1]

Cube.verts = [(i,j,k) for i in PM for j in PM for k in PM]
Cube.edges = [(i,j,k) for i in PZM for j in PZM for k in PZM if count(i,j,k) == 2]
Cube.faces = [(i,j,k) for i in PZM for j in PZM for k in PZM if count(i,j,k) == 1]

Cube.v2e = dict([(v,[e for e in Cube.edges if dist(v,e) == 1]) for v in Cube.verts])
Cube.v2f = dict([(v,[f for f in Cube.faces if dist(v,f) == 2]) for v in Cube.verts])
Cube.e2v = dict([(e,[v for v in Cube.verts if dist(e,v) == 1]) for e in Cube.edges])
Cube.e2f = dict([(e,[f for f in Cube.faces if dist(e,f) == 1]) for e in Cube.edges])
Cube.f2v = dict([(f,[v for v in Cube.verts if dist(f,v) == 2]) for f in Cube.faces])
Cube.f2e = dict([(f,[e for e in Cube.edges if dist(f,e) == 1]) for f in Cube.faces])
def face2rotator((i,j,k)):
  t = i+j+k
  if t%2:  (c,s) = (    0,   (-1)**(t-1))
  else:    (c,s) = ((-1)**t,      0     )
  if i != 0:
    rotator = array([[1 , 0 ,  0],
                     [0 , c , -s],
                     [0 , s ,  c],])
  elif j != 0:
    rotator = array([[ c , 0 , s],
                     [ 0 , 1 , 0],
                     [-s , 0 , c],])
  else: #k!=0
    rotator = array([[c , -s , 0],
                     [s ,  c , 0],
                     [0 ,  0 , 1],])
  return rotator
def rotatorPower(rotator, power):
  #assumes rotator is orthogonal matrix
  assert isinstance(power, int), "only integer powers are handled"
  if power > 0:
    return reduce(matrixmultiply, [rotator]*power)
  elif power < 0:
    return reduce(matrixmultiply, [transpose(rotator)]*(-power))
  else:
    return identity(rotator.shape[0])
def matrixExponential(A, steps=20):
  result = 1.0*identity(3)
  a = copy(result)
  for n in range(1,steps):
    a = (1.0/n) * matrixmultiply(A, a)
    result += a
  return result
epsilon3 = Marray([
  [[0,0,0],[0,0,-1],[0,1,0]],
  [[0,0,1],[0,0,0],[-1,0,0]],
  [[0,-1,0],[1,0,0],[0,0,0]]
])
def pureQuatExponential(q):
  return matrixExponential(sum([q[i]*epsilon3[i] for i in range(3)]))
def pointRotator(axis, theta):
  q = theta * Marray(list(axis))
  return pureQuatExponential(q)
def rotatePos(rotator, pos):
  return matrixmultiply(rotator, pos)
def rotateIndex(rotator, vect):
  return tuple(matrixmultiply(rotator,array(list(vect))))
def rotatePolyDict(rotator, polys):
  return dict([(rotateIndex(rotator,i), polys[i].rotate(rotator)) for i in polys])


#==========[ defining the little cube geometry ]==========
fillet = 0.2
F = 1.0-fillet
F2 = 1.0-fillet/2.0
F3 = 1.0-2.0*fillet/3.0
def vert2poly((i,j,k)):
  points = [Marray([ i  , F*j , F*k]),
            Marray([F*i ,  j  , F*k]),
            Marray([F*i , F*j ,  k ])]
  normal = Marray([i,j,k])
  center = F3*normal
  return Poly3D(points, normal, center)
def edge2poly((i,j,k)):
  if i == 0:
    points = [Marray([ F ,  j  , F*k]),
              Marray([-F ,  j  , F*k]),
              Marray([-F , F*j ,  k ]),
              Marray([ F , F*j ,  k ])]
  elif j == 0:
    points = [Marray([ i  ,  F , F*k]),
              Marray([ i  , -F , F*k]),
              Marray([F*i , -F ,  k ]),
              Marray([F*i ,  F ,  k ])]
  else: #k == 0
    points = [Marray([  i  , F*j ,  F]),
              Marray([  i  , F*j , -F]),
              Marray([ F*i ,  j  , -F]),
              Marray([ F*i ,  j  ,  F])]
  normal = Marray([i,j,k])
  center = F2*normal
  return Poly3D(points, normal, center)
def face2poly((i,j,k)):
  if i != 0:
    points = [Marray([i ,  F ,  F]),
              Marray([i ,  F , -F]),
              Marray([i , -F , -F]),
              Marray([i , -F ,  F])]
  elif j != 0:
    points = [Marray([ F , j ,  F]),
              Marray([ F , j , -F]),
              Marray([-F , j , -F]),
              Marray([-F , j ,  F])]
  else: #k != 0
    points = [Marray([ F ,  F , k]),
              Marray([ F , -F , k]),
              Marray([-F , -F , k]),
              Marray([-F ,  F , k])]
  normal = Marray([i,j,k])
  center = 1.0*normal
  return Poly3D(points, normal, center)

GCube = EmptyClass()
GCube.verts = dict([(v,vert2poly(v)) for v in Cube.verts])
GCube.edges = dict([(e,edge2poly(e)) for e in Cube.edges])
GCube.faces = dict([(f,face2poly(f)) for f in Cube.faces])

allBlack = dict([(f,0) for f in Cube.faces])


#==========[ geometric little cubes ]==========
class LittleCube:
  def __init__(self, pos, colors):
    self.pos = pos
    verts = deepcopy(GCube.verts)
    edges = deepcopy(GCube.edges)
    faces = deepcopy(GCube.faces)
    for face in Cube.faces:
      faces[face].color = colors[face]
    self.polys = verts.values() + edges.values() + faces.values()
    for poly in self.polys:
      poly.translate(pos)

  def projectToPlane(self, projector):
    polyList = []
    for poly in self.polys:
      poly2 = projector.projectPoly3D(poly)
      if poly2 is not None:  polyList.append(poly2)
    return polyList

  def rotate(self, rotator):
    self.pos = rotatePos(rotator, self.pos)
    for poly in self.polys:
      poly.rotate(rotator)

  def rotated(self, rotator):
    return [p.rotated(rotator) for p in self.polys]


#==========[ rotation splines ]==========
class Spline:
  def __init__(self, coeffs):
    self.coeffs = coeffs
    self.N = range(len(coeffs))

  def __call__(self, t):
    return sum([self.coeffs[n]*t**n for n in self.N])
smoothSpline = Spline([0.0, 0.0, 3.0,  -2.0]) #stationary frame


#==========[ combinatorial rubiks cube ]==========
dirList = [-2,-1,1,2]
class RCube:
  def __init__(self):
    self.littleCubes = {}
    for i in PZM:
      for j in PZM:
        for k in PZM:
          r = count(i,j,k)
          pos = 2.0*Marray([i,j,k])
          if r == 0: #center of cube
            continue
          if r == 1: #face centers
            face = (i,j,k)
            colorIndex = 1 + Cube.faces.index(face)
            colors = copy(allBlack)
            colors[face] = colorIndex
            self.littleCubes[face] = LittleCube(pos,colors)
          elif r == 2: #edges
            edge = (i,j,k)
            colors = copy(allBlack)
            for face in Cube.e2f[edge]:
              colorIndex = 1 + Cube.faces.index(face)
              colors[face] = colorIndex
            self.littleCubes[edge] = LittleCube(pos,colors)
          elif r == 3: #corners
            vert = (i,j,k)
            colors = copy(allBlack)
            for face in Cube.v2f[vert]:
              colorIndex = 1 + Cube.faces.index(face)
              colors[face] = colorIndex
            self.littleCubes[vert] = LittleCube(pos,colors)

  def stratifyIndices(self, axis):
    index2level = (lambda p: 1 + p[0]*axis[0] + p[1]*axis[1] + p[2]*axis[2])
    stratified = [[],[],[]]
    for index in self.littleCubes:
      stratified[index2level(index)].append(index)
    return stratified
    
  def stratifyCubes(self, axis):
    "returns a poly3D list for each rotation level"
    index2level = (lambda p: 1 + p[0]*axis[0] + p[1]*axis[1] + p[2]*axis[2])
    stratified = [[],[],[]]
    for (index,cube) in self.littleCubes.iteritems():
      stratified[index2level(index)].append(cube)
    return stratified

  def alignAxis(self,axis):
    return normalize(self.littleCubes[axis].pos)

  def rotate(self, axis, dirs):
    newCubes = {}
    levels = self.stratifyIndices(axis)
    averageDir = sum(dirs) / 3.0
    alignedAxis = self.alignAxis(axis)
    for level in range(3):
      idxRotator = rotatorPower(face2rotator(axis), dirs[level])
      theta = (pi/2.0) * (dirs[level]-averageDir)
      ptRotator = pointRotator(alignedAxis, theta)
      for index in levels[level]:
        cube = self.littleCubes[index]
        cube.rotate(ptRotator)
        newCubes[rotateIndex(idxRotator, index)] = cube
    self.littleCubes = newCubes

  def flattened(self):
    return Sum([cube.polys for cube in self.littleCubes.values()],[])

  #def projectToPlane(self, projector):
  #  return Sum([cube.projectToPlane(projector) for cube in self.littleCubes.values()],[])


#==========[ projecting camera ]==========
class Projector:
  def __init__(self,pos,X,Y,Z):
    self.pos = pos
    self.X = X
    self.Y = Y
    self.Z = Z

  def projectPoint(self, pos):
    rel = pos - self.pos
    x = innerproduct(self.X,rel)
    y = innerproduct(self.Y,rel)
    z = innerproduct(self.Z,rel)
    return Marray([x/z, y/z])

  def measureDepth(self, pos):
    rel = pos - self.pos
    z = innerproduct(self.Z,rel)
    #return z
    return norm(rel)

  def isVisible(self, pos, dir):
    rel = pos - self.pos
    youFaceMe = (innerproduct(rel, dir)      < 0)
    iFaceYou  = (innerproduct(rel, self.Z) > 0)
    return iFaceYou and youFaceMe

  def projectPoly3D(self, poly):
    if not self.isVisible(poly.center, poly.normal):
      return None
    else:
      points = [self.projectPoint(p) for p in poly.points]
      depth = self.measureDepth(poly.center)
      color = poly.color
      return Poly2D(points, depth, color)


#==========[ scene setup ]==========
colorTable = [
  color.cmyk.Black,
  color.cmyk.Bittersweet,
  color.cmyk.OliveGreen,
  color.cmyk.Blue,
  color.cmyk.Yellow,
  color.cmyk.Orange,
  color.cmyk.White,
  color.cmyk.Cyan
]

class Mirror:
  def __init__(self, pos, radius=1.0, numRings=24, numSecants=64, normal=None):
    self.pos = pos
    self.radius = radius
    self.numRings = numRings
    self.numSecants = numSecants
    if normal is None:  self.normal = normalize(-pos)
    else:               self.normal = normal
    inv = outerproduct(self.pos,self.pos)/innerproduct(self.pos,self.pos)
    self.reflect = 1.0*identity(3) - 2.0*inv

  def mirrorPoint(self, pos):
    return self.pos + matrixmultiply(self.reflect, (pos-self.pos))

  def mirrorProjector(self, proj):
    pos = self.mirrorPoint(proj.pos)
    X = matrixmultiply(self.reflect, proj.X)
    Y = matrixmultiply(self.reflect, proj.Y)
    Z = matrixmultiply(self.reflect, proj.Z)
    return Projector(pos,X,Y,Z)

  def projectToPlane(self, proj):
    dr = (1.0*self.radius)/self.numRings
    dt = 1.0/self.numSecants
    U = normalize(crossprod(self.pos,[0,0,1])) #fails for horizontal mirrors
    V = normalize(crossprod(self.pos,U))
    pos3D = lambda t: self.pos + dr*t*(cos(2*pi*t)*U + sin(2*pi*t)*V)
    pos2D = lambda t: proj.projectPoint(pos3D(t))
    N = 1 + self.numRings * self.numSecants
    lineColor = 6 #6=white, 7=cyan
    return Path2D([pos2D(n*dt) for n in range(N)], lineColor)

depthCmp = (lambda p,q: cmp(q.depth,p.depth))
class Scene:
  def __init__(self, polys, mirror, proj):
    nonNull = lambda p: p is not None
    mirroredProj = mirror.mirrorProjector(proj)
    bgPolys = filter(nonNull, [p.projectToPlane(mirroredProj) for p in polys])
    bgPolys.sort(depthCmp)
    spiral = mirror.projectToPlane(proj)
    fgPolys = filter(nonNull, [p.projectToPlane(proj) for p in polys])
    fgPolys.sort(depthCmp)
    self.elements = bgPolys + [spiral] + fgPolys

  def transform(self, shift, scale=1.0):
    for e in self.elements:
      e.transform(shift, scale)

  def draw(self, c):
    for e in self.elements:
      e.draw(c)


#==========[ rotation sequences ]==========
_axes = [(1,0,0),(0,1,0),(0,0,1)]
_levels = range(3)
_directions = [-2,-1,-1,1,1,2] #ones are twice as likely, since -2=2
def dirMod(dir):
  while True:
    if dir<-2:    dir += 4
    elif dir>=2:  dir -= 4
    else:         break
  return dir
def dirsMod(dirs):
  shifted = ( [[dirMod(d+i)   for d in dirs] for i in range(4)]
            + [[dirMod(d+i)+1 for d in dirs] for i in range(4)] )
  weigh = (lambda D: (D,norm(D)))
  weighed = [weigh(s) for s in shifted]
  weighed.sort(lambda x,y: cmp(x[1],y[1]))
  return weighed[0][0]

class Move:
  def __init__(self, axis, dirs):
    assert axis in _axes, "invalid axis: %s" % axis
    assert len(dirs) == 3, "wrong number of directions"
    for dir in dirs:
      assert dir is 0 or dir in _directions, "invalid direction: %s" % dir
    self.axis = axis
    self.dirs = dirs

  def __str__(self):
    return "%s: %s" % (_axes.index(self.axis), self.dirs)

  def __add__(lhs,rhs):
    assert lhs.axis == rhs.axis, "cannot add moves on different axes"
    dirs = dirsMod([lhs.dirs[level]+rhs.dirs[level] for level in _levels])
    return Move(lhs.axis, dirs)

  def duration(self):
    averageDir = sum(self.dirs) / 3.0
    return norm([dir-averageDir for dir in self.dirs])

  def actOnCube(self, cube):
    cube.rotate(self.axis, self.dirs)

  def reverse(self):
    self.dirs = [-dir for dir in self.dirs]

class MoveSequence:
  def __init__(self):
    self.moves = []
    self.durations = []
    self.times = []

  def __calcTimes(self):
    self.durations = [move.duration() for move in self.moves]
    self.times = cumsum(self.durations)

  def __str__(self):
    return "\n".join([str(move) for move in self.moves])

  def duration(self):
    if self.moves:  return self.times[-1]
    else:           return 0.0

  def randomlyAdvance(self, steps=20):
    for step in range(steps):
      axis = random.choice(_axes)
      level = random.choice(_levels)
      dirs = [0,0,0]
      dirs[level] = random.choice(_directions)
      self.moves.append(Move(axis, dirs))
    self.__calcTimes()

  def simplify(self):
    "combines consecutive moves on a single axis, reduces null rotations"
    moves = []
    prevAxis = None
    for move in self.moves:
      if move.axis == prevAxis:  moves[-1] = moves[-1] + move
      else:                      moves.append(deepcopy(move))
      prevAxis = move.axis
    self.moves = [m for m in moves if m.duration() > 0]
    self.__calcTimes()

  def reverse(self):
    self.moves.reverse()
    for move in self.moves:
      move.reverse()
    self.__calcTimes()

  def actOnCube(self, cube):
    for move in self.moves:
      move.actOnCube(cube)

  def animateOnCube(self, cube, mirror, proj, interval=0.25):
    "WARNING: this modifies the cube"
    assert interval < 1, "can't act that fast"
    Nframes = ceil(self.duration()/interval)
    t = 0.0
    T = Nframes*interval
    n = 0
    N = len(self.moves)
    print "animating %s moves over %s frames:" % (N,Nframes)
    scenes = []
    while t <= T:
      print "  drawing cube at (t,n) = (%s/%s, %s/%s)" % (t,T,n,N)
      while (n < N) and (t >= self.times[n]):
        self.moves[n].actOnCube(cube)
        n += 1
      if n < N:
        move = self.moves[n]
        if n == 0:  t0 = 0.0
        else:       t0 = self.times[n-1]
        tau = smoothSpline((t-t0) / self.durations[n])
        cubeLevels = cube.stratifyCubes(move.axis)
        complex3D = []
        averageDir = sum(move.dirs)/3.0
        q0 = cube.alignAxis(move.axis)
        for level in range(3):
          theta = (pi/2.0) * tau * (move.dirs[level] - averageDir)
          q = theta*q0
          rotator = pureQuatExponential(q)
          for littleCube in cubeLevels[level]:
            complex3D += littleCube.rotated(rotator)
      else:
        print "  drawing cube at (t,n) = (%s/%s, %s/%s)" % (t,T,n,N)
        complex3D = cube.flattened()
      scenes.append(Scene(complex3D, mirror, proj))
      t += interval
    return scenes
        

#==========[ paper resizing ]==========
def resize(old_c, new_width=16.75, new_height=10.75):
  b = old_c.bbox()
  old_height = unit.toinch(b.height())
  old_width  = unit.toinch(b.width())
  s = min(new_height/old_height, new_width/old_width)
  t = trafo.scale(s)
  new_c = canvas.canvas()
  new_c.insert(old_c, [t])
  return new_c


#==========[ default tests ]==========
def runPicture():
  #define cube
  cube = RCube()

  #define projector
  pos = Marray([20,15,10])
  scale = 30
  #Z = normalize(Marray([-1, -1,    -1  ]))
  #X = scale*normalize(Marray([-1, +1,     0  ]))
  #Y = scale*normalize(Marray([-1, -1, sqrt(2)]))
  Z = normalize(-pos)
  X = scale*normalize(crossprod(Z,[0,0,1]))
  Y = scale*normalize(crossprod(X,Z))
  proj = Projector(pos, X,Y,Z)
  
  #define mirror
  pos = Marray([-8,-16,-7])
  rad = 8.0
  mirror = Mirror(pos, rad)

  #define scene
  scene = Scene(cube.flattened(), mirror, proj)

  #draw scene
  c = canvas.canvas()
  scene.draw(c)
  r = 6.0
  c.stroke(path.path(path.moveto(r,r),
                     path.lineto(r,-r),
                     path.lineto(-2*r,-r),
                     path.lineto(-2*r,r),
                     path.lineto(r,r)))
  c.writeEPSfile("rubiks")

def runAnimation():
  #define move sequence
  moves = MoveSequence()
  moves.randomlyAdvance(1)
  moves.simplify()
  print "move sequence (axis, directions):"
  for move in moves.moves:
    print "%s, %s" % (move.axis, move.dirs)
  #define cube
  cube = RCube()
  moves.actOnCube(cube)
  moves.reverse()

  #define projector
  pos = Marray([20,0,0])
  scale = 6
  #Z = normalize(Marray([-1, -1,    -1  ]))
  #X = scale*normalize(Marray([-1, +1,     0  ]))
  #Y = scale*normalize(Marray([-1, -1, sqrt(2)]))
  Z = normalize(-pos)
  X = scale*normalize(crossprod(Z,[0,0,1]))
  Y = scale*normalize(crossprod(X,Z))
  proj = Projector(pos, X,Y,Z)

  #define mirror
  pos = Marray([-4,-4,0])
  rad = 9.0
  mirror = Mirror(pos, rad)

  #animate the scenes
  scenes = moves.animateOnCube(cube, mirror, proj)

  #draw the scenes
  print "assembling final drawing:"
  N = len(scenes)
  aspectRatio = 0.7
  I = int(round(sqrt(2.0*aspectRatio*N)))
  c = canvas.canvas()
  r = 1.2
  for n in range(N):
    (x,y) = (4.0*r*(n/I), 2.0*r*(I-(n%I)))
    print "  drawing box at (%s, %s)" % (x,y)
    c.text(x-2.9*r, y, str(n))
    s = scenes[n]
    s.transform(Marray([x,y]))
    s.draw(c)
    c.stroke(path.path(path.moveto(x+r,y+r),
                       path.lineto(x+r,y-r),
                       path.lineto(x-3*r,y-r),
                       path.lineto(x-3*r,y+r),
                       path.lineto(x+r,y+r)))
  c = resize(c) #changes paper size
  c.writeEPSfile("rubiks")

if __name__ == "__main__":
  runAnimation()

