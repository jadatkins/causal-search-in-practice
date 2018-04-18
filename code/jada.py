#! /usr/bin/env python -i
# -*- coding: utf_8 -*-

from copy import deepcopy
from random import *
from itertools import *
from enum import Enum
import numpy as np
import networkx as nx
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr
r = robjects.r
ggm = importr("ggm")          # R: library(ggm)
stats = importr("stats")      # R: library(stats)


if __name__ == '__main__':
  print(quit)

np.set_printoptions(precision=5, threshold=200, linewidth=110)


def random_dag(n, trace=False):
  if n < 1:
    raise(ValueError('n = ' + str(n)))

  G = nx.empty_graph(n, create_using=nx.DiGraph())

  rota = range(n)
  shuffle(rota)
  if trace:
    print 'rota =', rota

  if n == 1:
    return G

  for k in range(n-1):
    (i, j) = (randint(0, n-2), randint(1, n-1))
    if j > i:
      G.add_edge(rota[i], rota[j])

  for (u, v) in combinations(rota, 2):
    if not G.has_edge(u, v) and random() < 2.0 / ((n-1)*(n-1)):
      G.add_edge(u, v)

  if trace:
    print 'G.edges() =', G.edges(), '+ [',

  while not nx.is_weakly_connected(G):
    (i, j) = (randint(0, n-2), randint(1, n-1))
    if j > i and not G.has_edge(rota[i], rota[j]):
      if trace:
        print str((rota[i], rota[j])) + ',',
      G.add_edge(rota[i], rota[j])

  if trace:
    print ']'

  assert nx.is_directed_acyclic_graph(G)
  return G


def random_dist(G, nonnegative=False):
  if not nx.is_directed_acyclic_graph(G):
    raise(ValueError('G has cycles: ' + str(nx.simple_cycles(G))))

  for i in G.nodes_iter():
    G.node[i]['mu'] = choice((1, -1)) * lognormvariate(1, 0.2)
    G.node[i]['sigma'] = lognormvariate(0, 0.2)

  for (u, v) in G.edges_iter():
    w = lognormvariate(0, 0.2)
    if not nonnegative and random() < 0.5:
      w = -w
    G[u][v]['weight'] = w

  return


def sample(G, rows, trace=False):
  if not nx.is_directed_acyclic_graph(G):
    raise(ValueError('G has cycles: ' + str(nx.simple_cycles(G))))

  # rota specifies when the data are generated, in order from causes to effects
  rota = list()
  for (v, d) in G.in_degree_iter():
    if d == 0:
      rota.append(v)

  waiting = list()  # waiting to be added to rota; children of recent
  current = list()  # currently being added to rota; snapshot of waiting
  recent = rota     # recently added to rota; try to add children now

  while len(rota) < len(G):
    for parent in recent:
      for child in G.successors_iter(parent):
        if child not in rota and child not in waiting:
          waiting.append(child)
    recent = []
    current = deepcopy(waiting)  # we mustn't modify waiting
    for node in current:              # while stepping through it
      if all(parent in rota for parent in G.predecessors_iter(node)):
        rota.append(node)
        recent.append(node)
        waiting.remove(node)

  if trace:
    print 'rota =', rota

  # preprocessing done -- time to generate some data
  table = np.zeros((rows, len(G)))

  for row in range(rows):
    for v in rota:
      (mu, sigma) = (G.node[v]['mu'], G.node[v]['sigma'])
      for p in G.predecessors_iter(v):
        mu = mu + G.edge[p][v]['weight'] * table[row][p]
      table[row][v] = gauss(mu, sigma)

  return table


EdgeT = Enum('none', 'undirected', 'forward', 'back', 'both')

def edge_type(G, a, b):
  if not ( G.has_edge(a, b) or G.has_edge(b, a) ):
    return EdgeT.none

  if G.has_edge(a, b) and G.has_edge(b, a):
    if 'arrow' not in G[a][b] and 'arrow' not in G[b][a]:
      return EdgeT.undirected
    if 'arrow' in G[a][b] and 'arrow' in G[b][a]:
      return EdgeT.both

  if G.has_edge(a, b) and 'arrow' in G[a][b] and (
      not G.has_edge(b, a) or 'arrow' not in G[b][a] ):
    return EdgeT.forward
  if G.has_edge(b, a) and 'arrow' in G[b][a] and (
      not G.has_edge(a, b) or 'arrow' not in G[a][b] ):
    return EdgeT.back

  assert False


def orient(G, u, v):
  direction = edge_type(G, u, v)
  if direction == EdgeT.forward:
    return
  if direction not in (EdgeT.undirected, EdgeT.both):
    raise(RuntimeError(
      'Cannot orient '+ str(u) +' - '+ str(v) +' as '+ str(u) +' -> '+ str(v)
    ))
  G[u][v]['arrow'] = True
  G.remove_edge(v, u)


def phase_1(Gu, table, S, separating, threshold, trace, final):
  (rows, n) = table.shape

  for j in range(n-1):
    for (a, b) in combinations(Gu.nodes(), 2):  # for edges (a, b) in Gu
      if not Gu.has_edge(a, b):
        continue

      a_nb = Gu.neighbors(a)
      b_nb = Gu.neighbors(b)
      neighbours = ( set(a_nb) | set(b_nb) ) - set((a, b))

      if j > len(neighbours):
        continue

      for others in combinations(neighbours, j):
        if len(others) == 0:
          # R: pvalue = cor.test(table[,a], table[,b])$p.value
          pvalue = stats.cor_test(table[:,a], table[:,b])[2][0]
        else:
          u = (a+1, b+1) + tuple((o+1) for o in others)
          # R: parcor = pcor(u, S)
          parcor = ggm.pcor(IntVector(u), S)
          # R: pvalue = pcor.test(parcor, length(others), rows)$pvalue
          pvalue = ggm.pcor_test(parcor, len(others), rows)[2][0]

        if trace and (pvalue > threshold or (final and pvalue > 1e-9)):
          print (a, b), "|", others, "\t", "p =", pvalue, "\t",

        if pvalue > threshold: # if a and b independent conditional on others
          if trace:
            print "independent"
          Gu.remove_edge(a, b)
          separating.add_edge(a, b, s=others)  # others is a minimal set, but
          break                 # there could be more of the same cardinality
        elif trace and final and pvalue > 1e-9:
          print


def pc_search(table, threshold=0.05, firstpass=1.0/3.0, trace=False):

  (rows, n) = table.shape
  Gu = nx.complete_graph(n)       # Gu is undirected
  separating = nx.empty_graph(n)  # data structure to hold separating sets


  # Phase 1: find conditional independencies and delete edges

  # compute the estimated covariance matrix and store it in S
  S = r.var(table)            # R: S = var(table)

  if firstpass is not None:
    phase_1(Gu, table, S, separating, firstpass, trace, final=False)
  phase_1(Gu, table, S, separating, threshold, trace, final=True)

  Gd = nx.DiGraph(Gu)  # Gd is Gu, directed, with edges pointing both ways


  ## Phase 2: identify V-structures

  if n < 3:
    return Gd

  for (a, b, c) in permutations(Gu.nodes(), 3):
    if ( Gu.has_edge(a, b) and Gu.has_edge(c, b) and not Gu.has_edge(a, c) and
         b not in separating[a][c]['s'] and
         ('arrow' not in Gd[a][b] or 'arrow' not in Gd[c][b])
       ):
      if trace:
        print a,'->',b,'<-',c,'because',a,'_||_',c,'|',separating[a][c]['s']
      Gd[a][b]['arrow'] = True
      Gd[c][b]['arrow'] = True

  for (u, v) in Gu.edges_iter():
    direction = edge_type(Gd, u, v)
    if direction == EdgeT.forward:
      Gd.remove_edge(v, u)
    elif direction == EdgeT.back:
      Gd.remove_edge(u, v)


  ## Phase 3: orient remaining edges if possible

  change = True

  while change:
    change = False

    for (a, b, c) in permutations(Gd.nodes(), 3):
      if (edge_type(Gd, a, b) == EdgeT.forward    and
          edge_type(Gd, b, c) == EdgeT.undirected and
          edge_type(Gd, a, c) == EdgeT.none):
        if trace:
          print b,'->',c,'because',a,'->',b,'and no V-structure',a,'->',b,'<-',c
        orient(Gd, b, c)
        change = True

    for (a, b, c) in permutations(Gd.nodes(), 3):
      if (edge_type(Gd, a, b) == EdgeT.forward    and
          edge_type(Gd, b, c) == EdgeT.forward    and
          edge_type(Gd, a, c) in (EdgeT.undirected, EdgeT.both)):
        if trace:
          print a,'->',c,'because',a,'->',b,'->',c,'and cycles not allowed'
        orient(Gd, a, c)
        change = True

    if n < 4:
      continue

    for (a, b, c, d) in permutations(Gd.nodes(), 4):
      if (edge_type(Gd, a, b) in (EdgeT.undirected, EdgeT.both) and
          edge_type(Gd, b, c) == EdgeT.undirected and
          edge_type(Gd, c, d) == EdgeT.forward    and
          edge_type(Gd, d, a) == EdgeT.back       and
          edge_type(Gd, a, c) == EdgeT.none       and
          edge_type(Gd, b, d) in (EdgeT.undirected, EdgeT.both)):
        if trace:
          print b,'->',d,'by rule 3 with a =',str(a)+', c =',c
        orient(Gd, b, d)
        change = True

    for (a, b, c, d) in permutations(Gd.nodes(), 4):
      if (edge_type(Gd, a, b) == EdgeT.undirected and
          edge_type(Gd, b, c) in (EdgeT.undirected, EdgeT.both) and
          edge_type(Gd, c, d) == EdgeT.undirected and
          edge_type(Gd, d, a) == EdgeT.forward    and
          edge_type(Gd, a, c) != EdgeT.none       and
          edge_type(Gd, b, d) == EdgeT.none):
        if trace:
          print a,'->',b,'<-',c,'by rule 4 with d =',d
        orient(Gd, a, b)
        orient(Gd, c, b)
        change = True

  return Gd
