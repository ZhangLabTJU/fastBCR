import numpy as np
import random
import random
from ete3 import Tree
import numpy as np
import sys 

#===================================================================================
def pathToRoot(nodeName, tree):
	D = tree&nodeName
	# Get the path from B to the root
	node = D
	path = []
	while node.up:
		nn = node.name
		if node.name == '' or node.name.isdigit():
			nn = 'none'
		path.append(nn)
		node = node.up
	if path and path[len(path)-1] == 'none':
		path[len(path)-1] = 'naive'
	return path
#===================================================================================
def readNKTree(fileName, root=False):
	with open(fileName, 'r') as file:
		data = file.read().replace('\n', '')
	file.close()
	tree = Tree(str(data),  format=1)
	if root:
		rootedTree = Tree()
		A = rootedTree.add_child(name="naive")
		A.add_child(tree)
		return rootedTree
	return tree

#===================================================================================
def costTree(tree, labels, adjMatrix):
	ctotal = 0; count = 0
	for node in tree.traverse("preorder"):

		parent = node.up
		if parent:
			path = pathToRoot(node.name, tree)
			if path:
				i = 1; 
				pn = path[len(path)-1]
				while (i< len(path)):
					if path[i] in labels:
						pn = path[i]
						i = len(path)
					i+=1									

				if pn in labels and node.name in labels:
					cost = adjMatrix[labels.index(node.name)][labels.index(pn)]
					ctotal = ctotal + cost
					count +=1
	if count != len(labels):
		print ('ERROR missing nodes', count,  len(labels))
	return ctotal

#===================================================================================
def getCommonAncestorPath(tree, node):
	path = pathToRoot(node, tree)
	return path[1]

#===================================================================================
def getCommonAncestorPaths(pathA, pathB):
	for i in pathA:
		if i !='none' and i in pathB:
			return i
	return 'naive'

#===================================================================================
def findCommonAncestorLeaves(tree, labels):
	coupleNodes = {}
	for i in range(len(labels)):	  
		for j in range(i+1, len(labels)):
			pathi = pathToRoot(labels[i], tree)
			pathj = pathToRoot(labels[j], tree)
			if 	labels[i] in pathj:
				ancestor = labels[i]
			elif labels[j] in pathi:
				ancestor = labels[j]
			else:
				ancestor = getCommonAncestorPaths(pathi, pathj)
			coupleNodes[labels[i]+ '-'+ labels[j]] = ancestor
	return coupleNodes


#===================================================================================
def checkConsistence(tree, labels):
	seen = {}
	for node in tree.traverse("preorder"):
		if node.name not in seen.keys():
			if node.name in labels:
				seen[node.name] = True
			elif node.name != '':
				print ("ERROR ", node.name, " not in labels")
		else:
			print ("ERROR ", node.name, " several times")
	print (len(seen),len(labels))
	return len(seen)==len(labels)

#===================================================================================
def takeCostAB(adjMatrix, labels, a, b):
	#take cost value of two nodes in the adjMatrix
	ai = labels.index(a)
	bi = labels.index(b)
	return adjMatrix[ai][bi]
	
#===================================================================================
def getDistances(t):
	infoCost = ""
	for node in t.traverse("preorder"):
		children = node.get_children()
		for c in children:
			if (node.name):
				infoCost += node.name + "," + c.name + "," + str(t.get_distance(node, c)) + "\n"
	return infoCost
	
		


#===================================================================================
def trimming(tree, labels, adjMatrix):
	L = []
	#Generating a list of nodes in postorder	
	for node in tree.traverse("postorder"):
		#print (node.name)
		if node.name != '':
			L.append(node.name)

	for node in L:
		path = pathToRoot(node, tree)
		#print (path)
		if len(path) >= 3: #try to up a level
			y = takeCostAB(adjMatrix, labels, path[0], path[1]) #Take cost of parent
			x = takeCostAB(adjMatrix, labels, path[0], path[2]) #Take cost of granParent
			if x <= y: #Move node to the up level
				courrent = tree.search_nodes(name=path[0])[0]
				removed_node = courrent.detach()
				grandParent =  tree.search_nodes(name=path[2])[0]
				grandParent.add_child(removed_node)

	return tree


#===================================================================================
def editTree(t1, adjMatrix, labels):
	
	tr = Tree()
	ardColapsed = []
	compteurNodesInternes=1
	for node in t1.traverse("preorder"):
		if node.is_root():
			children = node.get_children()
			for n in children:
				tr.add_child(name=n.name)
		else:
			G = tr.search_nodes(name=node.name)
			if (G):
				children = node.get_children()
				for n in children:
					
					if n.name not in ardColapsed:
						colapse = colapseNodes(n.name, children, node.name, adjMatrix, labels, ardColapsed); #print ('==', colapse)
						if colapse:
							distances = updateDistances(node.name, n.name, colapse, adjMatrix, labels)
							N = G[0].add_child(name='None'+ str(compteurNodesInternes), dist=distances[0]) # Adds a empty branch or bifurcation
							compteurNodesInternes +=1
							distA = t1.get_distance(node, n); 
							
							compteur=1
							n1 = N.add_child(name=n.name, dist=distances[compteur]) # Adds current node
							ardColapsed.append(n.name)
							compteur += 1
							for c in colapse:
								n2 = N.add_child(name=c, dist=distances[compteur])
								ardColapsed.append(c) # Adds other nodes
								compteur +=1
						else:
							distA = t1.get_distance(node, n); 
							G[0].add_child(name=n.name, dist=distA)
	
	return tr
	

#===================================================================================
def updateDistances(parent, node, sisters, cost, labels):
	distances = []
	idNode = labels.index(node); idPar = labels.index(parent)
	costNone = 1
	for s in sisters:
		idS = labels.index(s)
		costNone = max(cost[idS][idPar] - cost[idNode][idPar], cost[idNode][idPar] - cost[idS][idPar], costNone)
		
	distances.append(costNone)
	distances.append(cost[idNode][idPar] - costNone)
	for s in sisters:
		idS = labels.index(s)
		distances.append(cost[idS][idPar] - costNone)
	return  distances

#===================================================================================
def colapseNodes(node, lnodes, parent, cost, labels, aldColapsed):
	colapse = []
	#print ('---> ', node, lnodes, parent, labels)
	idNode = labels.index(node); idPar = labels.index(parent)
	for i in lnodes:
		idI = labels.index(i.name)
		#print ('nodes ', node, i.name)
		if i.name != node and cost[idNode][idI] <= cost[idNode][idPar] and cost[idNode][idI] <= cost[idI][idPar]:
			#print ('cost ',node, i.name, cost[idNode][idI]);print ('cost P ',node, parent, cost[idNode][idPar])
			if i.name not in aldColapsed:
				colapse.append(i.name)
	#print ("colapse: ", colapse)
	return colapse
