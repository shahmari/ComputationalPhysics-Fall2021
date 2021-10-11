import random
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def InitalGraph(dim, OIP):
    Graph_ = nx.grid_2d_graph(dim, dim)
    for i in list(Graph_.nodes):
        if 1-OIP > random.random():
            Graph_.remove_node(i)
    return Graph_


def CheckPath(dim, P):
    G = InitalGraph(dim, P)
    for i in range(dim-1):
        for j in range(dim-1):
            if ((0, i) in G.nodes) and ((dim-1, j) in G.nodes):
                check = nx.has_path(G, (0, i), (dim-1, j))
                if check == True:
                    return 1
    return 0


Plist = np.linspace(0.4, 0.8, num=50)
MeanList = []
STDList = []
num = 0
for p in Plist:
    print("\r", num, end='')
    num += 1
    AllRuns = []
    for i in range(1000):
        AllRuns.append(CheckPath(50, p))
    MeanList.append(np.mean(AllRuns))
    STDList.append(np.std(AllRuns))


plt.figure(dpi=150)
plt.errorbar(Plist, MeanList, yerr=STDList, ecolor='r', color='k')
plt.scatter(Plist, MeanList, s=15, c='k')
plt.title('Avg of results for diffrent P (dimension = 50, 1000 runs)')
plt.xlabel('P')
plt.ylabel('Avg of results')
plt.savefig('../../Figs/Q3/Q3.pdf')
plt.show()
