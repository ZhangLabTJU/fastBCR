from networkx import *
import dendropy
import sys
from convertTrees import *

def nodeMatch(a : networkx.Graph, b : networkx.Graph) -> bool:
    result = False
    #print("node a: "+ str(a) + " " + a['label'] + " node b: "+ str(b) + " " + b['label'])

    if (not(str(a).startswith("seq") or str(b).startswith("seq"))):
        if (not(str(a) == "naive" or str(b) == "naive")):
            result = True
    elif (str(a) == str(b)):
        result = True

    return result

def edgeMatch(a, b) -> bool:
    result = False

    for node in a:
        print("node ")
        print(node)

#    if (a[0] == b[0] and a[1] == b[1]):
#        result = True

    return result

def compare(newickTreeFileA : str, newickTreeFileB : str) -> int:
    ncolTreeA = fromNewickToNCOL(newickTreeFileA)
    ncolTreeB = fromNewickToNCOL(newickTreeFileB)

    print("Edges in A: ", ncolTreeA.edges(data=False))
    print("Edges in B: ", ncolTreeB.edges(data=False))

    pathA = networkx.Graph()
    pathB = networkx.Graph()
    score = 0

    for node in ncolTreeA.nodes():
        if (not node.startswith("index")):
            pathA = networkx.shortest_path(ncolTreeA, 'naive', node)
            print("A" + str(pathA))
            pathGraphA = networkx.path_graph(pathA)
            for nodePath in pathGraphA:
                pathGraphA.nodes[nodePath]['label'] = str(nodePath) #pathA[0]

            pathB = networkx.shortest_path(ncolTreeB, 'naive', node)
            print("B" + str(pathB))
            pathGraphB = networkx.path_graph(pathB)
            for nodePath in pathGraphB:
                pathGraphB.nodes[nodePath]['label'] = str(nodePath) #pathB[0]

            score = score + networkx.graph_edit_distance(pathGraphA, pathGraphB, nodeMatch)
            print (score)

    print("GED Path based")
    print(str(score) + " (path comparison score)"); return 0;
    print("GED tree based, it can take many minutes...")
    print(str(networkx.graph_edit_distance(ncolTreeA, ncolTreeB, nodeMatch)) + " (tree comparison score - GED)")

    return 0

if __name__ == "__main__":
    glasMSTtree = sys.argv[0]

    if (len(sys.argv) < 3):
        raise Exception("Missing arguments: newick tree A file and newick tree B file.")
    elif (len(sys.argv) > 3):
        raise Exception("Olny THREE arguments are required: newick tree A file and newick tree B file.")

#    print("#Processing file " + sys.argv[1] + ".")
    compare(sys.argv[1], sys.argv[2])
