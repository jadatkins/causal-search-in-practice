from itertools import *
import sys
from multiprocessing import Pool, Lock
import numpy as np
import networkx as nx
from edgetypes import *
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr
r = robjects.r
ggm = importr("ggm")          # R: library(ggm)
stats = importr("stats")      # R: library(stats)
lock = Lock()


def phase_1(Gu, table, meta, args, final):
  (rows, n) = table.shape
  pool = Pool()
  try:
    for j in range(n-1):

      # multiple CPU stuff
      result = pool.imap_unordered(test_edge, (
        ((a, b), (Gu, table, meta[a][b], final, rows, n, j), args)
        for (a, b) in combinations(Gu.nodes(), 2)
      ), 4)

      for (a, b, trigger, pvalue, sepset) in result:
        if pvalue > meta[a][b]['pvalue'] and (final or trigger):
          meta[a][b]['pvalue'] = pvalue
          meta[a][b]['sepset'] = sepset
        if trigger and Gu.has_edge(a, b):
          Gu.remove_edge(a, b)

      # the above is somewhat equivalent to:
      # for (a, b) in combinations(Gu.nodes(), 2):
      #   test_edge(( (a, b), (Gu, table, meta[a][b], final, rows, n, j),
      #               (S, threshold, confidence, verbose) ))
      # except replace the line 'trigger = True' in test_edge()
      # with if Gu.has_edge(a, b): 'Gu.remove_edge(a, b)

  finally:
    try:
      lock.release()
    except Exception:
      pass
    pool.terminate()


def test_edge((
  (a, b),
  (Gu, table, meta_ab, final, rows, n, j),
  (S, threshold, confidence, verbose)
)):
  try:
    if (not Gu.has_edge(a, b)) and (j > len(meta_ab['sepset'])):
      return (a, b, False, 0, tuple())

    a_nb = set(Gu.neighbors(a))
    b_nb = set(Gu.neighbors(b))
    neighbours = (a_nb | b_nb) - set((a, b))

    if j > len(neighbours):
      return (a, b, False, 0, tuple())

    trigger = False
    meta_ab['sepset'] = tuple()

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

      if pvalue > meta_ab['pvalue']:
        if final or pvalue > threshold:
          meta_ab['pvalue'] = pvalue
          meta_ab['sepset'] = others

        if verbose and (pvalue > threshold or (final and pvalue > 1e-9)):
          lock.acquire(True)
          sys.stdout.flush()
          print "test", a, "_||_", b, "|", str(others)+":\t", "p =", pvalue,

        if pvalue > threshold: # if a and b independent conditional on others
          if verbose:
            print "\t", "=> independent"
            sys.stdout.flush()
            lock.release()
          trigger = True
          if not confidence:
            break
        elif verbose and final and pvalue > 1e-9:
          print
          sys.stdout.flush()
          lock.release()

    return (a, b, trigger, meta_ab['pvalue'], meta_ab['sepset'])

  except KeyboardInterrupt:
    pass


def causal_search(table, threshold=0.05, firstpass=None,
                  confidence=False, verbose=False):

  (rows, n) = table.shape
  Gu   = nx.complete_graph(n)   # Gu is undirected
  meta = nx.complete_graph(n)   # data structure to hold separating sets
  for (a, b) in meta.edges_iter():
    meta[a][b]['pvalue'] = 0    # highest p-value found so far


  ## Phase 1: find conditional independencies and delete edges

  # compute the estimated covariance matrix and store it in S
  S = r.var(table)            # R: S = var(table)

  if firstpass is not None:
    if verbose:
      print 'first pass:'
    phase_1(Gu, table, meta, (S, firstpass, confidence, verbose), final=False)
    if verbose:
      print 'second pass:'
  phase_1(Gu, table, meta, (S, threshold, confidence, verbose), final=True)

  Gd = nx.DiGraph(Gu)  # Gd is Gu, directed, with edges pointing both ways


  ## Phase 2: identify V-structures

  if n < 3:
    if confidence:
      return (Gd, meta)
    else:
      return Gd

  for (a, b, c) in permutations(Gu.nodes(), 3):
    if ( Gu.has_edge(a, b) and Gu.has_edge(c, b) and not Gu.has_edge(a, c) and
         b not in meta[a][c]['sepset'] and
         ('arrow' not in Gd[a][b] or 'arrow' not in Gd[c][b])
       ):
      if verbose:
        print a,'->',b,'<-',c,'because',a,'_||_',c,'|',meta[a][c]['sepset']
      Gd[a][b] = {'arrow': True, 'pvalue': meta[a][c]['pvalue']}
      Gd[c][b] = {'arrow': True, 'pvalue': meta[a][c]['pvalue']}

  for (u, v) in Gu.edges_iter():
    direction = edge_type(Gd, u, v)
    if direction == EdgeT.forward:
      Gd.remove_edge(v, u)
    elif direction == EdgeT.back:
      Gd.remove_edge(u, v)


  ## Phase 2a: fix orientation contradictions
  if confidence:
    for (u, v) in Gu.edges_iter():
      if edge_type(Gd, u, v) == EdgeT.both:
        if Gd[u][v]['pvalue'] > Gd[v][u]['pvalue']:
          if verbose:
            print u, '->', v, 'more likely than', u, '<-', v
          orient(Gd, u, v)
        elif Gd[v][u]['pvalue'] > Gd[u][v]['pvalue']:
          if verbose:
            print v, '->', u, 'more likely than', v, '<-', u
          orient(Gd, v, u)


  ## Phase 3: orient remaining edges if possible

  change = True

  while change:
    change = False

    for (a, b, c) in permutations(Gd.nodes(), 3):
      if (edge_type(Gd, a, b) == EdgeT.forward    and
          edge_type(Gd, b, c) == EdgeT.undirected and
          edge_type(Gd, a, c) == EdgeT.none):
        if verbose:
          print b, '->', c, 'because', a, '->', b, \
                'and no V-structure', a, '->', b, '<-', c
        orient(Gd, b, c)
        change = True

    for (a, b, c) in permutations(Gd.nodes(), 3):
      if (edge_type(Gd, a, b) == EdgeT.forward    and
          edge_type(Gd, b, c) == EdgeT.forward    and
          edge_type(Gd, a, c) in (EdgeT.undirected, EdgeT.both)):
        if verbose:
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
        if verbose:
          print b, '->', d, 'by rule 3 with [a, c] =', [a, c]
        orient(Gd, b, d)
        change = True

    for (a, b, c, d) in permutations(Gd.nodes(), 4):
      if (edge_type(Gd, a, b) not in (EdgeT.none, EdgeT.back) and
          edge_type(Gd, c, b) not in (EdgeT.none, EdgeT.back) and
          edge_type(Gd, c, d) in (EdgeT.undirected, EdgeT.both) and
          edge_type(Gd, d, a) == EdgeT.forward    and
          edge_type(Gd, a, c) != EdgeT.none       and
          edge_type(Gd, b, d) == EdgeT.none):
        if (edge_type(Gd, b, c) == EdgeT.both and
            edge_type(Gd, c, d) == EdgeT.both):
          continue
        if (edge_type(Gd, a, b) == EdgeT.forward and
            edge_type(Gd, c, b) == EdgeT.forward):
          continue
        if verbose:
          print a, '->', b, '<-', c, 'by rule 4 with d =', d
        orient(Gd, a, b)
        orient(Gd, c, b)
        change = True

  if confidence:
    return (Gd, meta)
  else:
    return Gd
