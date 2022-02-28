"""
While we should ensure expander properties,
for now just take a 2 rigular bipartite graph, b/c, w.h.p.,
its an expander:
Ref: https://theory.epfl.ch/courses/topicstcs/Lecture3.pdf
"""

import networkx as nx
import numpy as np
import random
import math
import torch
import hypernetx as hnx
import itertools
import sys
import json
from scipy import sparse


def createTannerGraph(regularity: tuple[int, int], n: int, m: int):
  (r, c) = regularity
  print("Creating bipartite graph with shape:", n, m)
  B = nx.bipartite.gnmk_random_graph(n, m, c * n)
  return B

# REGULARITY = (5,6) # defines a (r,c) biregular bipartite graph
def gen_hx_hz(k=5, regularity=(5, 6)):
	n = regularity[1] * k
	m = regularity[0] * k

	(g1, g2) = (createTannerGraph(regularity, n, m), createTannerGraph(regularity, n, m))
	hyperGraph = nx.cartesian_product(g1, g2)

	## Get all (i, j) (i, j') edges for H_x
	## This can be done by  iterating over the edges
	## of g1 (resp g2) and populating the matrix.
	## Then for each edge, check if (i, j) in g1 for edge (i, i), (i, j) where j in 0...m - 1
	numb_stabilizers = m * n
	numb_qubits = m * m + n * n
		
	def gen_H_X():
		H_x = torch.zeros(numb_stabilizers, numb_qubits, dtype=torch.int8)
		# for stabilizer_idx in range(m * n):
		# 	# 0 <= a < n
		# 	a = stabilizer_idx % n
		# 	# n <= beta < n + m
		# 	beta = int(stabilizer_idx / n) + n
		# 	stabilizerVertex = (a, beta)

		for x in range(m * m):
			dataVertex = (int(x / m) + n, x % m + n)
			for y in range(n):
				stabilizerVertex = (y, dataVertex[1])
				if hyperGraph.has_edge(dataVertex, stabilizerVertex):
					H_x[m * y + (dataVertex[1] - n)][n * n + x] = 1
		
		for x in range(n * n):
			dataVertex = (int(x / n), x % n)
			for y in range(m):
				stabilizerVertex = (dataVertex[0], y + n) ## hmm not sure
				if hyperGraph.has_edge(dataVertex, stabilizerVertex):
					H_x[m * dataVertex[0] + y][x] = 1

		maxDeg = 0
		for i in range(0, m *m):
			maxDeg = max(maxDeg, torch.count_nonzero(H_x[i]))
		print(maxDeg)

		return H_x
		
	def gen_H_Z():
		H = torch.zeros(numb_stabilizers, numb_qubits, dtype=torch.int8)
		for x in range(m * m):
			dataVertex = (int(x / m) + n, x % m + n)
			for y in range(n):
				stabilizerVertex = (dataVertex[0], y)
				if hyperGraph.has_edge(dataVertex, stabilizerVertex):
					H[m * y + (dataVertex[0] - n)][n * n + x] = 1
		
		for x in range(n * n):
			dataVertex = (int(x / n), x % n)
			for y in range(m):
				stabilizerVertex = (y + n, dataVertex[1])
				if hyperGraph.has_edge(dataVertex, stabilizerVertex):
					H[m * dataVertex[1] + y][x] = 1

		maxDeg = 0
		for i in range(0, m * m):
			maxDeg = max(maxDeg, torch.count_nonzero(H[i]))
		print(maxDeg)
		return H



	H_X = gen_H_X()
	H_Z = gen_H_Z()
	return (H_X, H_Z, numb_qubits, numb_stabilizers)

def save_H(H_X, H_Z, N_QUBITS, N_STABLE, outfile, is_matlab=True):
	triple_savable = {}
	def to_triple_format(H_sparse, is_matlab):
		offset = 1 if is_matlab else 0
		I = [int(x + offset) for x in H_sparse.row]
		J = [int(x + offset) for x in H_sparse.col]
		V = [1 for x in I]
		return I, J, V
	
	H_X_sparse = sparse.coo_matrix(H_X)
	H_Z_sparse = sparse.coo_matrix(H_Z)

	Ix, Jx, Vx = to_triple_format(H_X_sparse, is_matlab)
	triple_savable["Ix"] = Ix
	triple_savable["Jx"] = Jx
	triple_savable["Vx"] = Vx

	Iz, Jz, Vz = to_triple_format(H_Z_sparse, is_matlab)
	triple_savable["Iz"] = Iz
	triple_savable["Jz"] = Jz
	triple_savable["Vz"] = Vz
	triple_savable["N_Qubits"] = N_QUBITS
	triple_savable["N_Stabilizers"] = N_STABLE
	with open(outfile, 'w') as outfile:
			json.dump(triple_savable, outfile)
	
# TODO: test for expander property??
for i in range(10, 15):
	regularity=(2, 3)
	(H_X, H_Z, N_QUBITS, N_STABLE) = gen_hx_hz(k=i, regularity=regularity)
	save_H(H_X, H_Z, N_QUBITS, N_STABLE, f"../matlab/g-{i * regularity[1]}-{i * regularity[0]}.json")