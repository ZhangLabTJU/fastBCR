from networkx import *
from igraph import *
import dendropy
import numpy
import re



#CONSTANTS -> to be a config file in the near future
GLAmst_COL_DELIMITER_WANNA_BE = " "
COL_DELIMITER = "\t"
NCOL_TERMINATION = ".ncol"
SVG_TERMINATION = ".svg"

def fromNewickToNCOL(newickTreeFile : str) -> networkx.Graph:
    ncolTreeFromNewick = networkx.Graph()
    newickTree = dendropy.Tree.get(path=newickTreeFile, schema='newick')

    counter = 0
    for node in newickTree.nodes():
        if (not node is None):
            if (str(node.label) == "None" and str(node.taxon) == "None"):
                node.label = "index" + str(counter)
                counter = counter + 1

    for edge in newickTree.edges():
        if (not edge.head_node is None and not edge.tail_node is None):
            if (not str(edge.tail_node.label) == 'None'):
                tail = edge.tail_node.label
            elif (not str(edge.tail_node.taxon) == 'None'):
                tail = str(edge.tail_node.taxon)
#            elif (not str(edge.tail_node) == 'None'):
#                tail = str(edge.tail_node)
            else:
                continue

            if (not str(edge.head_node.label) == 'None'):
                head = str(edge.head_node.label).replace("'","")
            elif (not str(edge.head_node.taxon) == 'None'):
                head = str(edge.head_node.taxon).replace("'","")
#            elif (not str(edge.head_node) == 'None'):
#                head = edge.head_node
            else:
                continue

            ncolTreeFromNewick.add_edges_from([(head,tail)])

    for node in ncolTreeFromNewick.nodes():
        ncolTreeFromNewick.nodes[node]['label'] = str(node)

    return ncolTreeFromNewick

def fromNCOLtoNewick(ncolTreeFile: str) -> dendropy.Tree:
    parent = -1
    descendant = -1
    newickTree = ""
    descendants = ""
    edges = []
    ncolGraph = numpy.array([0,0])

    with open(ncolTreeFile) as treeFile:
        for row in treeFile:
            edges = row.strip().split("\t")
            parent = str(edges[0])
            descendant = str(edges[1])
            ncolGraph = numpy.vstack([ncolGraph, [parent, descendant]])

    ncolGraph = numpy.delete(ncolGraph, 0, axis = 0)
    ncolGraph.view('i8,i8').sort(order=['f0'], axis=0)

    ncolGraphR = Graph.Read_Ncol(ncolTreeFile)
    ncolGraphR.write(ncolTreeFile + ".svg", "svg")

    parentSet = set(ncolGraph[:,0])
    descendantSet = set(ncolGraph[:,1])
    root = parentSet - descendantSet
    descendantSet = root.copy()
    newickTree = " " + str(root.pop()) + " ;"

    while (len(descendantSet) > 0):
        parent = descendantSet.pop()
        for edge in ncolGraph:
            if (parent == edge[0]):
                descendants = str(edge[1]) + " , " + descendants
                descendantSet.add(edge[1])
                #remover a edge de ncolGRaph
        if (len(descendants) > 0):
            newickTree = newickTree.replace(" " + str(parent) + " ", " ( " + descendants[:-2] + " ) " + str(parent) + " ")
        descendants = ""

    print(newickTree)

    return dendropy.Tree.get(data = newickTree, schema = "newick")



if __name__ == "__main__":
    if (len(sys.argv) < 2):
        raise Exception("Missing argument: ncol tree file.")
    elif (len(sys.argv) > 2):
        raise Exception("Olny ONE argument is required: ncol tree file.")

#    print("#Processing file " + sys.argv[1] + ".")
    fromNCOLtoNewick(sys.argv[1])
