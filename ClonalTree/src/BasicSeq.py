import numpy as np
import random
import sys 
from Bio import SeqIO

INF = float('inf') 

#---------------------------------------------------------------------------
# Basic sequences - here you find basic functions to deal with sequences
#---------------------------------------------------------------------------

def readFastaAbundance(fastaFile):
	dico = {}
	naive = ""; SEP="@"
	labels = [];  arraySeqs = []; Abundance = {}
	count = 0; root  = 0

	for record in SeqIO.parse(fastaFile, "fasta"):
		
		if SEP in record.id:
			ID = record.id.split(SEP)[0]
			abund = int(record.id.split(SEP)[1])
		else:
			ID = record.id
			abund = 1
		
		if ID not in dico.keys(): 
			dico[ID] = str(record.seq)
			labels.append(ID)
			arraySeqs.append(str(record.seq))
			if ID == 'naive':
				root = count
			count = count + 1
		if ID in Abundance.keys():
			Abundance[ID] += abund
		else:
			Abundance[ID] = abund
		
	return labels, root, arraySeqs, Abundance, dico

#---------------------------------------------------------------------------
def hamming_distance(chaine1, chaine2):
	return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

#---------------------------------------------------------------------------
#create adjacent matrix from colpased sequences by using hamming distance
def createAdjMatrix(arraySeqs):
	adjMatrix = np.zeros((len(arraySeqs), len(arraySeqs)))

	for i in range(len(arraySeqs)):
		#for j in range(i+1, len(arraySeqs)):
		for j in range(0, len(arraySeqs)):
			adjMatrix[i][j] = hamming_distance(arraySeqs[i], arraySeqs[j])
	return adjMatrix

#===================================================================================
def correctMatrix(adjMatrixNP, visitedNodes):
   
	for i in range(len(adjMatrixNP)):
		for j in range(i, len(adjMatrixNP[i])):
			if i in visitedNodes and j in visitedNodes:
				adjMatrixNP[i][j] = INF; adjMatrixNP[j][i] = INF
	return adjMatrixNP


#---------------------------------------------------------------------------
def readFasta(fastaFile):
	dico = {}
	naive = ""
	labels = [];  arraySeqs = [];
	count = 0; root  = 0

	for record in SeqIO.parse(fastaFile, "fasta"):
		if record.id not in dico.keys(): 
			dico[record.id] = str(record.seq)
			labels.append(record.id)
			arraySeqs.append(str(record.seq))
			if record.id == 'naive':
				root = count
			count = count + 1
		
	return labels, root, arraySeqs

#---------------------------------------------------------------------------
def readFasta2(fastaFile):
	dico = {}
	naive = ""
	labels = [];  arraySeqs = [];
	count = 1; root  = 0

	for record in SeqIO.parse(fastaFile, "fasta"):
		if record.id not in dico.keys(): 
			dico[record.id] = str(record.seq)
			labels.append(record.id)
			arraySeqs.append(str(record.seq))
			count = count + 1
		
	return labels, dico

#---------------------------------------------------------------------------
def readFastaRepeat(fastaFile):
	dico = {}; 
	naive = "" ; maps = "";
	labels = [];  arraySeqs = []; Abundance = {}
	count = 1; root  = 0

	for record in SeqIO.parse(fastaFile, "fasta"):
		seq = str(record.seq)
		if seq not in dico.keys(): 
			dico[seq] = str(record.id)
			labels.append(record.id)
		else:
			dico[seq] += "," + str(record.id)
	
	for k,v in dico.items():
		if "naive" in v:
			ID = "naive"
			root = count - 1
		else:
			ID = "seq" + str(count)
			count +=1
		labels.append(ID)
		arraySeqs.append(k)
		Abundance[ID] = v.count(",")
		maps += ID + "\t" + v + "\n"	
	f = open(fastaFile +'.maps', 'w'); f.write(maps); f.close()	
	return labels, root, arraySeqs, Abundance






