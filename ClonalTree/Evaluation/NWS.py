# -*- coding: utf8 -*-
import sys		# for sys.stdout.write()


########################################################################
###Global variables
fastaFile = ""
substMatrixFile = ""
outputFile= ""


#Usage
usage = "python NWS.py -f fastaFile -sM substitution  Matrix File -o outputFile\n"


########################################################################
### Read parameters
def readParameters(args):
	global fastaFile
	global substMatrixFile
	global outputFile
	
	
	for i in range(1,len(args)):
		if (args[i] == "-f"):
			fastaFile = args[i+1]
		elif (args[i] == "-sM"):
			substMatrixFile = args[i+1]
		elif (args[i] == "-o"):
			outputFile = args[i+1]
		elif (args[i] == "-h"):
			print (usage)
########################################################################
### Check parameters
def checkParameters():
	if (fastaFile == ""):
		print ("ERROR::Parameter -f fastaFile is required\n")
		sys.exit(1);
	elif (substMatrixFile == ""):
		print ("ERROR::Parameter -sM substitution  Matrix File is required\n")
		sys.exit(1);
	elif (outputFile == ""):
		print ("ERROR::Parameter -o outputFile is required\n")
		sys.exit(1);
	


###########################################
# 1 Affichage de l'alignement
### sortie d'alignement s1,s2, et string d'alignement
def saveAli(id1, id2, s1, s2, sa):
	f = open(outputFile, "w")
	
	ii = 0
	f.write(id1 + ' ')
	#print sa
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'g':
			f.write(s1[ii])
			ii = ii + 1
		else:
			f.write('-')
	f.write('\n')
	f.write(id2 + " ")
	ii = 0
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'h':
			f.write(s2[ii])
			ii = ii + 1
		else:
			f.write('-')
	f.write('\n')
	f.close()

###########################################
# 1 Affichage de l'alignement
### sortie d'alignement s1,s2, et string d'alignement
def printali(s1,s2,sa):
	ii = 0
	#print sa
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'g':
			sys.stdout.write(s1[ii])
			ii = ii + 1
		else:
			sys.stdout.write('-')
	print( '')
	ii = 0
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'h':
			sys.stdout.write(s2[ii])
			ii = ii + 1
		else:
			sys.stdout.write('-')
	print ('')

###########################################
# 1 Affichage de l'alignement
### sortie d'alignement s1,s2, et string d'alignement
def saveAli2(id1, id2, s1, s2,sa):
	ii = 0
	info = id1 + " ";
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'g':
			info = info + s1[ii]
			ii = ii + 1
		else:
			info = info + '-'
	print (info)
	sys.exit(1);
 
	info = info + "\n" + id2 + " ";
	ii = 0
	for i in range(len(sa)):
		if sa[i] == 'd' or sa[i] == 'h':
			info = info + s2[ii]
			ii = ii + 1
		else:
			 info = info + '-'
	print (info)

#######################################
#2 NWS aligment
#s1 est en haut (donc en colonne)
#s2 est vertical, donc en ligne

### alignement iteratif
def alignit(s1,s2,gap=-1,subs={}):
	gapH = -100; gapG=0; gapI=-100
	# matrice des distances
	m = list(range(len(s2)+1))
	# matrice des chemins
	for i in range(len(s2)+1):
		m[i] = list(range(len(s1)+1))
	# root cell
	m[0][0] = (0, 'o')
	# first line
	for j in range(1,len(s1)+1):
		v=m[0][j-1][0] + gapI
		m[0][j]=(v, 'g')
	# first column
	for i in range(1,len(s2)+1):
		m[i][0]=(m[i-1][0][0] + gapI , 'h')	
	# tab
	for i in range(1,len(s2)+1):
		for j in range(1,len(s1)+1):
			distd = substitue(s2[i-1], s1[j-1], subs) + m[i-1][j-1][0]
			disth = gapH + m[i-1][j][0]
			distg = gapG + m[i][j-1][0]
			if distd >= disth and distd >= distg:
				m[i][j] = (distd,'d')
#				c[i][j] = 'd'  # substitution
			elif disth >= distd and disth >= distg:
				m[i][j] = (disth, 'h')
#				c[i][j] = 'h' # insertion
			else:
				m[i][j] = (distg, 'g')
#				c[i][j] = 'g' # deletion
	#printMatAli(m, s1, s2)
	return (m)
 #


###################
#Ecriture matrice alignement
def printMatAli(matAli, s1="", s2=""):
	ss1="0"+s1
	ss2="0"+s2
	if s1!="" and s2 != "":
		sys.stdout.write("	 ")
		for i in ss1:
			sys.stdout.write("   %c  "%i)
		sys.stdout.write("\n")
	for i in range(0,len(matAli)):
		if s1!="" and s2 != "":
			sys.stdout.write("  %c  "%ss2[i])
		for j in range(0,len(matAli[i])):
			sys.stdout.write("%3d:%c "%(matAli[i][j][0], matAli[i][j][1]))
		sys.stdout.write("\n")



 
############################
def compScore(path, s1, s2, matrix, gap=0):
	score = 0; k=0; j=0
	for i in range(len(path)):
		#print (i, k, j)
		if path[i] == 'd':
			#print ('sb ', s1[k], s2[j], substitue(s1[k], s2[j], matrix))
			score = score + substitue(s1[k], s2[j], matrix)
			k +=1; j+=1
		else:
			score += gap
			if path[i] == 'g':
				k = k + 1

			else:
				j = j + 1
	return score


############################
#Retrouver le chemin
def backtrack(m):
	#print "distance = ",m[len(s1)][len(s2)]
	#print m
	#print c
	# backtrack
	#nb col
	j = len(m[0])-1
	#nb lignes
	i = len(m)-1
	chemin = ''
	while (i != 0 or j != 0):
		if i < 0 or j < 0:
			print ("backtrack:: ERROR i or j <0",i,j)
			exit(1)
		#print i,j,m[i][j]
		if m[i][j][1] == 'd':
			i = i-1
			j = j-1
			chemin = 'd'+chemin
		elif m[i][j][1] == 'g':
			if len (m[i][j]) ==3:
				#print "Saut de ",m[i][j][2]
				for k in xrange(m[i][j][2]):
					j=j-1 
					chemin = 'g'+chemin
			else:
				j = j-1
				chemin = 'g'+chemin
		else:
			if len (m[i][j]) ==3:
				#print "Saut de ",m[i][j][2]
				for k in xrange(m[i][j][2]):
					i = i-1
					chemin = 'h'+chemin
			else:
				i = i-1
				chemin = 'h'+chemin
		#print i,j
	#print chemin
	return chemin


#####################################
#3 Lecture matrice
##lit et retourne une matrice de similarite, ne pas oublier l'affectation devant
### argument: nom de fichier
def readsimmatrix(nomFi):
	f= open(nomFi)
	lignes=f.readlines()
	f.close()
	#print lignes
	#lignes = filter(lambda x:x[0] != '#', lignes)
	#print lignes
	#lignes = map(lambda x:x.strip().split(), lignes)
	#print lignes
	mats = {}
	numL=0
	for l in lignes:
		if l[0] != '#':
			c=l[:-1].split()
			#print "ligne",c
			if numL == 0:
				listeAA2=c
				numL=1
			else:
				aa1=c[0]
				for j in range(1,len(c)):
					mats[(aa1,listeAA2[j-1])]=int(c[j])
	print (mats)
	return mats
 



###########################
def substitue(res1, res2, subs):
	if (type(subs) is dict):
		if subs=={}:				
			# Par default match 0, mismatch -1
			score = -(res2 != res1) 
		else:
			# match/mismatch pris dans la matrice de distance "matrice"
			score= subs[res1,res2]
	else:
		if (type(subs) is int or type(subs) is float):
			if(res1 != res2):
				score = subs 
			else:
				score =  0
		else:
			print ("substitue:: pb de type avec la variable subs", subs, type(subs))
			sys.exit(1)
	return score
	
###########################
#Lecture d'un fichier de plusieurs sequence au format fasta
def readFastaMul(nomFi):
	#Lecture du contenu du fichier
	f=open(nomFi,"r")
	lines=f.readlines()
	f.close()
	#Parsage du du contenu du fichier
	seq=""
	nom=""
	lesSeq=[]
	for l in lines:
		if l[0] == '>':
			if seq != "":
				tmp=(nom,seq)
				lesSeq.append(tmp)
			nom=l[1:-1]
			seq=""
		else:
			#Ne pas oublier d'enlever le retour chariot a la fin des lignes
			seq=seq+l[:-1]
	if seq != "":
		tmp=(nom,seq)
		lesSeq.append(tmp)
	return lesSeq

###########################
########################################################################
#####################		 Main	  ################################
########################################################################
'''
readParameters(sys.argv)
checkParameters()

subs = readsimmatrix(substMatrixFile);
seqs = readFastaMul(fastaFile);
print (seqs)
#print subs
matrix = alignit(seqs[0][1], seqs[1][1], -5, subs);
path = backtrack(matrix)
id1 = seqs[0][0]; id2 = seqs[1][0];
seq1 = seqs[0][1]; seq2 = seqs[1][1];

saveAli(id1, id2, seq1, seq2, path)

'''

