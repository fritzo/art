
from exceptions import *
from pyx import *
import random

from sets import Set
Set.__add__ = Set.__or__
Set.__iadd__ = Set.__ior__

#==========[ go board ]==========
def areSequential(i1, i2):
  return abs(i1-i2) <= 1

def areConnected(pos1, pos2):
  return areSequential(pos1[0], pos2[0]) or areSequential(pos1[1], pos2[1])

dotPositions = {
  19:sum([[(3+6*i,3+6*j) for i in range(3)] for j in range(3)], []),
  13:sum([[(2+4*i,2+4*j) for i in range(3)] for j in range(3)], []),
   9:sum([[(2+4*i,2+4*j) for i in range(2)] for j in range(2)], [])
}

colorMap = {
  1 : color.rgb.black,
  2 : color.rgb.white
}

class GoBoard:
  def __init__(self, size=19):
    self.size = size
    self.states = dict(sum([[((i,j),0) for i in range(size)]
                                   for j in range(size)], []))
    self.pos2group = {}
    self.group2pos = {}
    self.groupCounter = 0 #just for naming

  #playing
  def play(self, pos, state):
    assert pos in self.states, "invalid position: %s" % pos
    assert self.states[pos] == 0, "cannot play on top of stone"
    assert state in [1,2], "invalid play state: %s" % state
   #change state
    self.states[pos] = state
   #merge with surrounding groups
    good, evil = Set(), Set()
    for nbh in self.getNeighbors(pos):
      nbhState = self.states[nbh]
      if nbhState:
        if nbhState == state:  good.add(self.pos2group[nbh])
        else:                  evil.add(self.pos2group[nbh])
   #kill surrounded groups
    for group in evil:
      if not self.getLiberties(group):
        self.removeGroup(group)
   #possibly commit suicide
    if good:
      group = self.mergeGroups(good)
      self.pos2group[pos] = group
      self.group2pos[group].add(pos)
    else:
      group = self.makeNewGroup(pos)
    if not self.getLiberties(group):
      #print "suicide doesn't pay"
      self.removeGroup(group)
    print self

  def randomlyPlay(self, numSteps):
    positions = list(self.states)
    for step in range(numSteps):
      state = 1 + step%2
      pos = random.choice(positions)
      while self.states[pos]:
        pos = random.choice(positions)
      self.play(pos, state)

  def getNeighbors(self,pos):
    (x0,y0) = pos
    nbhs = []
    if x0 > 0:            nbhs.append((x0-1,y0))
    if x0 < self.size-1:  nbhs.append((x0+1,y0))
    if y0 > 0:            nbhs.append((x0,y0-1))
    if y0 < self.size-1:  nbhs.append((x0,y0+1))
    return nbhs

  def getLiberties(self, group):
    postions = self.group2pos[group]
    liberties = sum([[nbh for nbh in self.getNeighbors(pos)
                          if self.states[nbh] == 0]
                          for pos in postions], [])
    return Set(liberties)

  def makeNewGroup(self, pos):
    group = self.groupCounter
    self.groupCounter += 1
    self.group2pos[group] = Set([pos])
    self.pos2group[pos] = group
    return group

  def mergeGroups(self, groups):
    merged = min(groups)
    positions = sum([self.group2pos[group] for group in groups], Set())
    self.group2pos[merged] = positions
    for pos in positions:
      self.pos2group[pos] = merged
    return merged

  def removeGroup(self, group):
    for pos in self.group2pos[group]:
      del self.pos2group[pos]
      self.states[pos] = 0

  #display, both ascii & poscript
  def __str__(self):
    symbols = {0:'+', 1:'X', 2:'O'}
    states = [[self.states[i,j] for j in range(self.size)]
                                for i in range(self.size)]
    state2str = lambda line: '--'.join(map(symbols.__getitem__, line))+'\n'
    lines = [state2str(line) for line in states]
    between = (self.size-1)*'|  ' + '|\n'
    return between.join(lines)

  def draw(self):
    c = canvas.canvas()
    I = self.size
   #draw background lines
    for i in range(I):
      v = path.path(path.moveto(i,0),
                    path.lineto(i,I-1))
      c.stroke(v, [style.linewidth(0.02)])
      h = path.path(path.moveto(0,i),
                    path.lineto(I-1,i))
      c.stroke(h, [style.linewidth(0.02)])
   #draw vertex dots
    if I in dotPositions:
      for pos in dotPositions[I]:
        p = path.circle(pos[0], pos[1], 0.1)
        c.fill(p)
   #draw foreground stones
    for i in range(self.size):
      for j in range(self.size):
        state = self.states[i,j]
        if state:
          color = colorMap[state]
          p = path.circle(i, j, 0.40)
          c.fill(p, [color])
          c.stroke(p, [style.linewidth(0.05)])
   #draw border
    c.stroke(path.path(path.moveto(          -0.5,           -0.5),
                       path.lineto(-0.5+self.size,           -0.5),
                       path.lineto(-0.5+self.size, -0.5+self.size),
                       path.lineto(          -0.5, -0.5+self.size),
                       path.lineto(          -0.5,           -0.5)))
    return c


#==========[ runtime tests ]==========
def drawBoard():
  b = GoBoard()
  b.randomlyPlay(19*19*4)
  c = b.draw()
  c.writeEPSfile("go_board")

if __name__ == "__main__":
  drawBoard()

