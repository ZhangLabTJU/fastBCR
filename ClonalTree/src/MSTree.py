import numpy as np
import random
import math
from operator import xor
from ete3 import Tree
import numpy as np
import sys 
from BasicSeq import *

INF = float('inf') 

#---------------------------------------------------------------------------
# MSTree - here you functions to deal with trees
#---------------------------------------------------------------------------

#===================================================================================
def primMST(cost, root, labels, abundance, useAb=True): 
	infoTree = ""
	tree = Tree()
	tree.add_child(name=labels[root])
	adjMatrixNP = np.array(cost); np.fill_diagonal(adjMatrixNP, INF)
	visitedNodes = [root]
	it = 0
	while len(visitedNodes) < len(labels):
		if useAb:
			minsI, minsJ =  aminIndex(adjMatrixNP, visitedNodes)
		else:
			minsI, minsJ =  aminIndexFirstFound(adjMatrixNP, visitedNodes)
		nodeA, nodeB = chooseBestNode(minsI, minsJ, visitedNodes, adjMatrixNP, labels, abundance)
		minV = adjMatrixNP[nodeA][nodeB]
		adjMatrixNP[nodeA][nodeB] = INF; adjMatrixNP[nodeB][nodeA] = INF
	
		tree, infoTree = addNodeTree(tree, labels[nodeA], labels[nodeB], minV, infoTree)
		if nodeA not in visitedNodes:
			visitedNodes.append(nodeA)
		if nodeB not in visitedNodes:
			visitedNodes.append(nodeB)

		it +=1
		
		adjMatrixNP = correctMatrix(adjMatrixNP, visitedNodes)
	
	return tree, infoTree

#===================================================================================
def aminIndex(matrix, indices):
	minV = INF
	mins = []; minI = []; minJ = []
	#Find the minimal value in the matrix constrainted by indices
	for i in indices:
		m = np.amin(matrix[i]) #take the min value
		if m < minV:
			minV = m
	
	#Find the j indices with minimal value
	for i in indices:
		for j in range(len(matrix[i])):
			if matrix[i][j] == minV:
				minI.append(i)
				minJ.append(j)
	return minI, minJ
#===================================================================================
def aminIndexFirstFound(matrix, indices):
	minV = INF
	mins = []; minI = []; minJ = []
	#Find the minimal value in the matrix constrainted by indices
	for i in indices:
		m = np.amin(matrix[i]) #take the min value
		j = np.argmin(matrix[i]) #take the index of the min value
		if m < minV:
			minV = m
			minJ = j
			minI = i
	return [minI], [minJ]

#===================================================================================
def chooseBestNode(minsI, minsJ, visitedNodes, adjMatrix, labels, abundance):

	maxAb = -INF; nodeA = 0; nodeB = 0
	
	for i in range(len(minsI)):
		a = minsI[i]; b = minsJ[i];
		
		if (labels[a]!=labels[b]) and (xor(a in visitedNodes, b in visitedNodes)):
			ab = abundance[labels[a]] + abundance[labels[b]]
			
			if  ab > maxAb:
				maxAb = ab
				nodeA = a
				nodeB = b
	if nodeA==0 and nodeB ==0:
		print ("ERROR: Disconnected Tree", len(visitedNodes))
		sys.exit()
	return nodeA, nodeB

#===================================================================================
def addNodeTree(t, a, b, min, infoCost):

	tp = ()
	G = t.search_nodes(name=a)
	
	if (G):
		G[0].add_child(name=b, dist=min) 
		infoCost += a + "," + b + "," + str(min) + "\n"
	else:
		G = t.search_nodes(name=b)
		if (G):
			G[0].add_child(name=a, dist=min)
			infoCost += a + "," + b + "," + str(min) + "\n"
		else:
			print ("Warnning nodes do not exists: ", a, b)
	return t, infoCost






































