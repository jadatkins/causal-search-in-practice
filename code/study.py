#! /usr/bin/env python

from copy import deepcopy
from random import *
from itertools import *
from time import clock
import numpy as np
import networkx as nx
from edgetypes import *
from causes import causal_search


def blank():
  return {
    'pscl': {'existence': 0.0, 'orientation': 0.0, 'time': 0.0},
    'twop': {'existence': 0.0, 'orientation': 0.0, 'time': 0.0},
    'conf': {'existence': 0.0, 'orientation': 0.0, 'time': 0.0},
    'both': {'existence': 0.0, 'orientation': 0.0, 'time': 0.0}
  }


def few_tests(graphs_per_test=1, distns_per_g=3, verbosity=2):
  accu_one = blank()
  print '-' * 62
  print "\nreally simple graph:"
  Gtrue = nx.empty_graph(3, create_using=nx.DiGraph())
  nodes = range(3)
  shuffle(nodes)
  [a, b, c] = nodes
  Gtrue.add_edges_from([(a, b), (c, b)])
  print 'Gtrue.edges() =', Gtrue.edges()
  test_one_graph(Gtrue, accu_one, distns_per_g, 250, verbosity)
  print_accuracy(accu_one, 1)

  accu_two = blank()
  for k in range(graphs_per_test):
    print '-' * 62
    print "\nrule 3 test, graph", str(k+1) + ':'
    Gtrue = nx.empty_graph(4, create_using=nx.DiGraph())
    nodes = range(4)
    shuffle(nodes)
    [a, b, c, d] = nodes
    Gtrue.add_edges_from([(a, d), (b, d), (c, d)])
    Gtrue.add_edge(*choice([(a, b), (b, a)]))
    Gtrue.add_edge(*choice([(b, c), (c, b)]))
    print '(a, b, c, d) =', (a, b, c, d)
    print 'Gtrue.edges() =', Gtrue.edges()
    test_one_graph(Gtrue, accu_two, distns_per_g, 10000, verbosity)
  print_accuracy(accu_two, graphs_per_test)

  accu_three = blank()
  for k in range(graphs_per_test):
    print '-' * 62
    print "\nrule 4 test, graph", str(k+1) + ':'
    Gtrue = nx.empty_graph(5, create_using=nx.DiGraph())
    nodes = range(5)
    shuffle(nodes)
    [a, b, c, d, e] = nodes
    Gtrue.add_edges_from([(a, b), (c, b), (d, a)])
    Gtrue.add_edges_from([(c, a), (c, d)])
    Gtrue.add_edge(e, a)  # required in order to detect (d, a)
    # prevents V-structure e -> b <- c:
    Gtrue.add_edge(*choice([(c, e), (e, c)]))
    # a -> b will still be detected by rule 1 due to d -> a
    print '(a, b, c, d, e) =', (a, b, c, d, e)
    print 'Gtrue.edges() =', Gtrue.edges()
    assert nx.is_directed_acyclic_graph(Gtrue)
    test_one_graph(Gtrue, accu_three, distns_per_g, 40000, verbosity)
  print_accuracy(accu_three, graphs_per_test)

  accu_four = blank()
  for k in range(graphs_per_test):
    print '-' * 62
    print "\ncombined test, graph", str(k+1) + ':'
    Gtrue = nx.empty_graph(8, create_using=nx.DiGraph())
    nodes = range(8)
    shuffle(nodes)
    [a3, b3, c3, d3, a4, b4, c4, d4] = nodes
    # start with the rule 3 graph
    Gtrue.add_edges_from([(a3, d3), (b3, d3), (c3, d3)])
    Gtrue.add_edge(*choice([(a3, b3), (b3, a3)]))
    Gtrue.add_edge(*choice([(b3, c3), (c3, b3)]))
    # add the rule 4 graph
    Gtrue.add_edges_from([(a4, b4), (c4, b4), (d4, a4)])
    Gtrue.add_edges_from([(c4, a4), (c4, d4)])
    # adding these two edges between the two subgraphs cannot produce a cycle:
    u = choice([a3, b3, c3, d3])
    Gtrue.add_edges_from([(u, a4), choice([(c4, u), (u, c4)])])
    print '(a3, b3, c3, d3, a4, b4, c4, d4) =', \
          (a3, b3, c3, d3, a4, b4, c4, d4)
    print 'Gtrue.edges() =', Gtrue.edges()
    assert nx.is_directed_acyclic_graph(Gtrue)
    test_one_graph(Gtrue, accu_four, distns_per_g, 90000, verbosity)
  print_accuracy(accu_four, graphs_per_test)


def many_tests(min_n=6, max_n=15, graphs_per_n=3, distns_per_g=3, verbosity=1):
  number_of_tests = graphs_per_n * len(range(min_n, max_n+1))
  accuracy = blank()
  for n in range(min_n, max_n+1):
    for j in range(graphs_per_n):
      print '-' * 62
      if verbosity >= 1:
        print
      print "n =", n, "\tgraph", str(j+1) + ':'
      Gtrue = random_dag(n)
      print 'Gtrue.edges() =', Gtrue.edges()
      test_one_graph(Gtrue, accuracy, distns_per_g, None, verbosity)
  print_accuracy(accuracy, number_of_tests)


def print_accuracy(accuracy, number_of_tests):
  print '-' * 62
  print '#' * 14, 'Overall Accuracy:', '#' * 14
  for (key, message) in [
    ('pscl', 'basic PC algorithm:'),
    ('twop', 'two-pass PC algorithm:'),
    ('conf', 'using confidence levels to make best guesses:'),
    ('both', 'using two passes and confidence levels:')
  ]:
    print
    print message
    print 'mean edge existence accuracy   =',
    print accuracy[key][ 'existence' ] / number_of_tests
    print 'mean edge orientation accuracy =',
    print accuracy[key]['orientation'] / number_of_tests
    print 'total time taken (seconds)     =',
    print accuracy[key]['time']


def test_one_graph(Gtrue, accuracy, batches=1, rows=None, verbosity=1):
  n = Gtrue.number_of_nodes()
  num_edges = Gtrue.number_of_edges()
  if rows is None:
    rows = int(n**(5.0/3.0) * 100)
  if verbosity >= 2:
    verbose = True
  else:
    verbose = False

  for batch in range(1,batches+1):
    if verbosity >= 1:
      print
    print "n =", n, "\trows =", rows, "\tbatch", str(batch)+':'
    random_dist(Gtrue)
    if verbosity >= 1:
      print 'Gtrue.nodes(data=True) =', Gtrue.nodes(data=True)
      print 'Gtrue.edges(data=True) =', Gtrue.edges(data=True)
      print
    data = sample(Gtrue, rows)

    for (test_fn, key) in (
      (test_pc  , 'pscl'),
      (test_twop, 'twop'),
      (test_conf, 'conf'),
      (test_both, 'both')
    ):  # Now test_fn is one of the four functions test_pc, etc
      # and accuracy[key] is one of accuracy['pscl'], etc

      accuracy[key]['time'] = accuracy[key]['time'] - clock()
      Guess = test_fn(data, verbosity, verbose)
      accuracy[key]['time'] = accuracy[key]['time'] + clock()
      if verbosity >= 1:
        print 'Guess.edges() =', Guess.edges()

      # count edge existence errors
      gtrue_edges = set(Gtrue.to_undirected().edges())
      guess_edges = set(Guess.to_undirected().edges())
      exist_errs = len(gtrue_edges ^ guess_edges)
      if verbosity >= 1:
        print 'edge existence accuracy   =', \
              str(num_edges - exist_errs) + '/' + str(num_edges),
      if num_edges > 0:
        exist_accu = float(num_edges - exist_errs) / float(num_edges)
        if verbosity >= 1:
          print "\t=", exist_accu,
        accuracy[key]['existence'] = \
          accuracy[key]['existence'] + exist_accu / batches
      else:
        accuracy[key]['existence'] = \
          accuracy[key]['existence'] + 1.0 / batches
      if verbosity >= 1:
        print

      # determine edge orientation accuracy
      direction_good = 0
      direct_max = num_edges
      for (u, v) in Gtrue.edges_iter():
        direction = edge_type(Guess, u, v)
        if direction == EdgeT.forward:
          direction_good = direction_good + 1
        elif direction == EdgeT.back:
          direction_good = direction_good - 1
        elif direction == EdgeT.none:
          direct_max = direct_max - 1
      if verbosity >= 1:
        print 'edge orientation accuracy =', \
              str(direction_good) + '/' + str(direct_max),
      if direct_max > 0:
        direct_accu = float(direction_good) / float(direct_max)
        if verbosity >= 1:
          print "\t=", direct_accu,
        accuracy[key]['orientation'] = \
          accuracy[key]['orientation'] + direct_accu / batches
      else:
        accuracy[key]['orientation'] = \
          accuracy[key]['orientation'] + 1.0 / batches
      if verbosity >= 1:
        print "\n"

  return accuracy


def test_pc(data, verbosity=2, verbose=True):
  if verbosity >= 1:
    print 'basic PC algorithm:'
  return causal_search(data, 0.05, None, False, verbose)

def test_twop(data, verbosity=2, verbose=True):
  if verbosity >= 1:
    print 'two-pass PC algorithm:'
  return causal_search(data, 0.05, 1.0/3.0, False, verbose)

def test_conf(data, verbosity=2, verbose=True):
  if verbosity >= 1:
    print 'using confidence levels to make best guesses:'
  (guess, meta) = causal_search(data, 0.05, None, True, verbose)
  return guess

def test_both(data, verbosity=2, verbose=True):
  if verbosity >= 1:
    print 'using two passes and confidence levels:'
  (guess, meta) = causal_search(data, 0.05, 1.0/3.0, True, verbose)
  return guess


def random_dag(n, verbose=False):
  if n < 1:
    raise(ValueError('n = ' + str(n)))

  G = nx.empty_graph(n, create_using=nx.DiGraph())

  rota = range(n)
  shuffle(rota)
  if verbose:
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

  if verbose:
    print 'G.edges() =', G.edges(), '+ [',

  while not nx.is_weakly_connected(G):
    (i, j) = (randint(0, n-2), randint(1, n-1))
    if j > i and not G.has_edge(rota[i], rota[j]):
      if verbose:
        print str((rota[i], rota[j])) + ',',
      G.add_edge(rota[i], rota[j])

  if verbose:
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


def sample(G, rows, verbose=False):
  if not nx.is_directed_acyclic_graph(G):
    raise(ValueError('G has cycles: ' + str(nx.simple_cycles(G))))

  for v in G:  # let v be an arbitrary node in G
    if ('mu' not in G.node[v] and 'sigma' not in G.node[v]):
      raise(TypeError('G has no distribution. Call random_dist(G) first.'))
    break

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
    current = deepcopy(waiting)  # we can't modify waiting
    for node in current:         # while stepping through it
      if all(parent in rota for parent in G.predecessors_iter(node)):
        rota.append(node)
        recent.append(node)
        waiting.remove(node)

  if verbose:
    print 'rota =', rota
    if rows > 10000:
      print 'generating data...',

  # preprocessing done -- time to generate some data
  table = np.zeros((rows, len(G)))

  for row in range(rows):
    for v in rota:
      (mu, sigma) = (G.node[v]['mu'], G.node[v]['sigma'])
      for p in G.predecessors_iter(v):
        mu = mu + G.edge[p][v]['weight'] * table[row][p]
      table[row][v] = gauss(mu, sigma)

  if verbose and rows > 10000:
    print 'done.'

  return table


if __name__ == '__main__':
  print 'few_tests(graphs_per_test=1, distns_per_g=1, verbosity=2):'
  few_tests(1, 1, 2)
  print "\n" + ('=' * 62) + "\n"
  print 'many_tests(min_n=4, max_n=23, graphs_per_n=3,', \
        'distns_per_g=3, verbosity=0):'
  many_tests(4, 23, 3, 3, 0)
