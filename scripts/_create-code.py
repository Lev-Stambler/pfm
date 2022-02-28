from venv import create
from networkx.generators import margulis_gabber_galil_graph
import math

from sympy import Ge

"""
Create an expander graph: note it is not a bipartite graph

The graph creation follows from https://networkx.org/documentation/networkx-1.10/reference/generated/networkx.generators.expanders.margulis_gabber_galil_graph.html
The graph has degree 8 and a second largest eigen value of 5 sqrt(2)

"""
def createExpanderGraph(n):
	G = margulis_gabber_galil_graph(int(n))
	return G


"""
Take in an expander graph and get the set of edges
"""
def getExpanderEdges(G):
	edges = {}
	for i in G._adj.items():
		v1 = i[0]
		for edge in i[1].items():
			# if list(edge[1].keys())[1] != 0:
			v2 = edge[0]
			edges[(v1, v2)] = 1
	return edges

"""
Takes in edges of the form (v1, v2), where v1, v2 = (x, y), x, y \in [n]

And returns a set of edges for G' = (L U R, E) following
https://cstheory.stackexchange.com/questions/31751/d-regular-bipartite-expander-graph,
to create a bipartite expander graph

The vertices are now indexed by n^2, where both L, R have n^2 vertices.

with f_{vertex transform} being ((x, y) => x * n + y. n^n will also be added to the right edges
"""
def getBipartiteGraphEdges(edges, n):	
	edgesBipartite = {}
	for x in range(n):
		for y in range(n):
			for xP in range(n):
				for yP in range(n):
					edge = ((x, y), (xP, yP))
					if edge in edges:
						newEdge = (x * n + y, xP * n + y + n**2)
						edgesBipartite[newEdge] = 1
	return edgesBipartite

def getBipartiteGraphRoutine(n):
	GExpander = createExpanderGraph(10)
	expEdges = getExpanderEdges(GExpander)
	return getBipartiteGraphEdges(expEdges, n)

"""
Note that here n is the number of vertices in any one side of the
bipartite partition. For simplicity, the sizes are the same.

First it will take graphs with edges of form a,b \in ([n], [n]) where a \in L and b \in R
and add n to all indices in b. 

returns two lists of edges for Hz, Hx
"""
def hypergraph2BipartiteCode(G1Edges, G2Edges, n):

	stabilizerZQubitPairs = {}

	# build Hz by having 100 <= b < 200, 0 <= \alpha < 100
	for b in range(100, 200):
		for alpha in range(0, 100):
			stabilizer_vertex = (b, alpha)
			for a in range(0, 100):
				# CHECK FOR VERTEX
				if (a, b) in G1Edges:
					stabilizerZQubitPairs[stabilizer_vertex, (a, alpha)] = 1
			for beta in range(100, 200):
				if (alpha, beta) in G2Edges:
					stabilizerZQubitPairs[stabilizer_vertex, (b, beta)] = 1

	# # build Hx by having 100 <= beta < 200, 0 <= a < 100
	# for beta in range(100, 200):
	# 	for a in range(0, 100):
	# 		stabilizer_vertex = (a, beta)
	# 		for alpha in range(0, 100):
	# 			# CHECK FOR VERTEX
	# 			if (a, b) in G1Edges:
	# 				stabilizerZQubitPairs[stabilizer_vertex, (a, alpha)] = 1
	# 		for beta in range(100, 200):
	# 			if (alpha, beta) in G2Edges:
	# 				stabilizerZQubitPairs[stabilizer_vertex, (b, beta)] = 1


	return (stabilizerZQubitPairs, stabilizerZQubitPairs)

"""
Create both parity check matrices (X, Z) from
G1, G2
"""
def createH(G1, G2):
	pass

"""
Save the parity check matrix
"""
def saveH(G, path="./expander-code.json"):
	pass

def main():
	n = 9
	G1 = getBipartiteGraphRoutine(n)
	G2 = getBipartiteGraphRoutine(n)
	(stabilizerZPairs, stabilizerXPairs) = hypergraph2BipartiteCode(G1, G2, n * n)

main()