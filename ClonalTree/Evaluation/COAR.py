from Bio import SeqIO
import numpy as np
from optparse import OptionParser
import random
from ete3 import Tree
import sys 
sys.path.insert(1, '../src')

from BasicTree import *
from BasicSeq  import *
from NWS import *

def makeBoolean(var):
	if var == '0':
		return True
	else:
		return False

			
#===================================================================================
def getAllPathLeaves(tree, labels):
	paths = []
	for leaf in labels:
		pathTree = pathToRoot(leaf, tree)
		if "none" in pathTree:
			pathTree = list(filter(("none").__ne__, pathTree))
		paths.append(pathTree)
	return paths

#===================================================================================
def comparePathsTree(paths1, paths2):
	newPaths = []
	for path1 in paths1:
		if path1 not in paths2:
			newPaths.append(path1)
	return newPaths

#===================================================================================
def mapping(paths, hashPath, index):
	mapped = []
	for path in paths:
		newPath = []
		for i in path:
			if i not in hashPath.keys():
				hashPath[i]=chr(index); index +=1;
			newPath.append(hashPath[i])
		mapped.append(newPath)
	return index, hashPath, mapped

#===================================================================================
def sbsMatrix(hashPath, arraySeqs, labels):
	D = len(hashPath.keys())
	#matrix = np.zeros((D, D)); i = 0;
	matrix = {}; minDist = 1000
	for k1, v1 in hashPath.items():
		j = 0
		for k2, v2 in hashPath.items():
			i1 = labels.index(k1)
			i2 = labels.index(k2)
			chaine1 = arraySeqs[i1]
			chaine2 = arraySeqs[i2]
			#matrix[i][j] = (-1)*hamming_distance(chaine1, chaine2)
			dist = (-1)*hamming_distance(chaine1, chaine2) 
			matrix[(v1, v2)] = dist
			if dist < minDist:
				minDist = dist
	return matrix, minDist
#===================================================================================
def computeCOAR(mapped1, mapped2, matrix, minDist):
	scoreT = 0
	#print ("md = ", minDist)
	for seq1 in mapped1:
		minScore = 1000
		for seq2 in mapped2:
			strS1 = ''.join(seq1)
			strS2 = ''.join(seq2)
			#print ("aln ", strS1, strS2)
			seqAln = [("s1", strS1), ("s2", strS2)]
			matrixAln = alignit(strS1, strS2, minDist-1, matrix);
			path = backtrack(matrixAln)
			#print (path)
			sc = compScore(path, strS1, strS2, matrix, minDist)
		   
			score =  sc/float(minDist); 
			#print ('sc = ', sc, 'scNorm ', score)
			if score < minScore:
				minScore = score
			#print ('Minscore', minScore)
		scoreT += minScore
	return scoreT
#===================================================================================
#						Main
#===================================================================================
def main():
	usage = usage = "python MRCA.py -a <nkTree1> -b <nkTree2> -f <fasta> \n"
	parser = OptionParser(usage)
	parser.add_option("-a", "--nkTree1", dest="nkTree1", help="nk file for gcTree")
	parser.add_option("-b", "--nkTree2", dest="nkTree2", help="nk file for clonalTree")
	parser.add_option("-f", "--fasta", dest="fastaFile", help="fasta file")
	
	(options, args) = parser.parse_args()
	if len(sys.argv) < 7:
		parser.error("incorrect number of arguments")
	
	
	nkTree1 = options.nkTree1
	nkTree2 = options.nkTree2
	fastaFile = options.fastaFile
	
	aR = "0"
	bR = "1"

	aRooted = makeBoolean(aR)
	bRooted = makeBoolean(bR)

	
	GCTree = readNKTree(nkTree1, aRooted) #GCTree
	#print (GCTree.get_ascii(show_internal=True))

	clonalTree = readNKTree(nkTree2, bRooted) #clonal Tree
	#print (clonalTree.get_ascii(show_internal=True))

	labels, root, arraySeqs, abundance, dico =  readFastaAbundance(fastaFile)
	
	pathsClonalTree = getAllPathLeaves(clonalTree, labels)
	pathsGCTree = getAllPathLeaves(GCTree, labels)

	totalPaths = len(pathsGCTree)
	#print ('totalPaths ', totalPaths)

	pathsCTClean = comparePathsTree(pathsClonalTree, pathsGCTree)
	#print ("CT= ", pathsCTClean)

	pathsGCClean = comparePathsTree(pathsGCTree, pathsClonalTree)
	#print ("GC= ",pathsGCClean)


	index, hashPath, mappedGC = mapping(pathsGCClean, {}, 97)
	#print(mappedGC)
	
	index, hashPath, mappedCL = mapping(pathsCTClean, hashPath, index)   	

	matrix, minDist = sbsMatrix(hashPath, arraySeqs, labels)
	#print (matrix, minDist)

	print ("COAR= ", computeCOAR(mappedGC, mappedCL, matrix, minDist)/float(totalPaths))


#===================================================================================
if __name__ == "__main__":
	main()

