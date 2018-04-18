from enum import Enum
import networkx as nx


EdgeT = Enum('none', 'undirected', 'forward', 'back', 'both')

def edge_type(G, a, b):
  if not isinstance(G, (nx.DiGraph, nx.MultiDiGraph)):
    raise(TypeError('G is not a directed graph.'))

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
