"""local assortativity statistics"""
from numpy import *
import scipy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
G1 = nx.barabasi_albert_graph(200, 10, seed=2)  #N=200, 10, 7, 2
G2 = nx.barabasi_albert_graph(300, 15, seed=1)  #N=300, 15, 10, 1
G3 = nx.barabasi_albert_graph(500, 20, seed=1)  #N=500, 20, 10 ,1
G4 = nx.barabasi_albert_graph(1000, 20, seed=1) #N=1000，20, 10, 1;


def attribution_node(G, node):
    '''attribution of a node'''
    neighbours = list(G.adj[node])
    neighbours_degree = []
    for i in range(len(neighbours)):
        neighbours_degree.append(G.degree(neighbours[i]))
    neighbours_dict = dict(zip(neighbours, neighbours_degree))
    #print('neighbours dict:'+str(neighbours_dict))
    return neighbours_dict


def nodes_min_degree(G):
    '''find the nodes with minimum degree of the network
    in scale free network there exists a lot of nodes with minimum degrees'''
    a = G.degree()
    b = sorted(a, key=lambda t: t[1])
    degree_sequence = [d for n, d in G.degree()]
    min_degree = min(degree_sequence) #minimum degree of the network
    #print(min_degree)
    min_node1 = []
    for i in range(len(b)):
        if b[i][1] == min_degree:
            min_node1.append(b[i][0])
    #print('nodes with min degree:'+str(min_node1))
    return min_node1


def estimate_eig_sec(G):
    lap_matrix = nx.laplacian_matrix(G).todense()
    eigval, eigvec = scipy.linalg.eig(lap_matrix)
    d = eigval.sort()
    #print(eigval)
    eig_sec = eigval[1]
    #print('laplacian eigen:'+str(eig_sec))
    r = nx.degree_assortativity_coefficient(G)
    #print(r)
    return eig_sec


def link_xmatch(G, node1, node2, node3, node4):
    '''1-2;3-4
        1-3;2-4'''
    G.remove_edge(node1, node2)
    G.add_edge(node1, node3)
    G.remove_edge(node3, node4)
    G.add_edge(node2, node4)
    return G


def key_link_reconection(G, n1, n2):
    neighbour1 = list(G.adj[n1])
    #print(neighbour1)
    neighbour2 = list(G.adj[n2])
    a = 0
    if n1 not in neighbour2:
        i = 0
        while a <1:
            neighbour3 = list(G.adj[neighbour2[i]])
            neighbour3.extend(neighbour2)
            neighbour4 = list(set(neighbour3))
            if neighbour1[i] not in neighbour4:
                link_xmatch(G, n1, neighbour1[i], n2, neighbour2[i])
                a +=1
            else:
                i+=1
                if i>= len(neighbour1):
                    break
    return G


def reconnection1(G, a, b):
    min_node = nodes_min_degree(G)
    #print(len(min_node))
    for i in range(a):
        for j in range(a, b):
            G = key_link_reconection(G, min_node[i], min_node[j])
    eig = float(estimate_eig_sec(G))
    #print(eig)
    r = nx.degree_assortativity_coefficient(G)
    #print(r)
    return G, eig, r
#G1 = reconnection1(G1, 10, 20)[0]


def difference(G, node):
    neighbours = list(G.adj[node])
    dif = []
    for i in range(len(neighbours)):
        dif.append(abs(G.degree(neighbours[i])-G.degree(node)))
    difvalue = sum(dif)/G.degree(node)
    #print(difvalue)
    return difvalue


def assortativity(G, node):
    N = G.number_of_nodes()
    list1 = []
    r = nx.degree_assortativity_coefficient(G)
    f = (r + 1) / N
    for i in range(N):
        list1.append(difference(G, i))
    scaler_dif = difference(G, node)/sum(list1)
    noder = f - scaler_dif
    #print(noder)
    return noder


def average_assor(G):
    min_node = nodes_min_degree(G)
    list1 = []
    for i in range(len(min_node)):
        list1.append(assortativity(G, min_node[i]))
    avediff = np.mean(list1)
    #print(avediff)
    return avediff
#average_assor(G1)


# N=200
list11 = []
list12 = []
for i in range(7):
    eig = reconnection1(G1,i,2*i)[1]
    #r = reconnection1(G1,i,2*i)[2]
    G1 = reconnection1(G1, i, 2*i)[0]
    noder = average_assor(G1)
    list11.append(eig)
    list12.append(noder)
#print(list11)
#print(list12)


# N=300
list21 = []
list22 = []
for i in range(10):
    eig = reconnection1(G2,i,2*i)[1]
    G2 = reconnection1(G2, i, 2*i)[0]
    noder = average_assor(G2)
    list21.append(eig)
    list22.append(noder)
#print(list21)
#print(list22)

# N=500
list31 = []
list32 = []
for i in range(10):
    eig = reconnection1(G3, i, 2*i)[1]
    G3 = reconnection1(G3, i, 2*i)[0]
    noder = average_assor(G3)
    list31.append(eig)
    list32.append(noder)


# N=1000
list41 = []
list42 = []
for i in range(10):
    eig = reconnection1(G4, i,2*i)[1]
    G4 = reconnection1(G4, i, 2*i)[0]
    noder = average_assor(G4)
    list41.append(eig)
    list42.append(noder)

fig = plt.figure(figsize=(8, 6))
#fig = plt.figure()
plt.subplot(2, 2, 1)
plt.plot(list12, list11, label=r'$N=200$')
plt.xlabel('local assortativity')
plt.ylabel(r'$\lambda^{(2)}$')
plt.gca().ticklabel_format(axis='x',style='scientific',scilimits=(-1,1),useMathText=True)
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(list22, list21, label=r'$N=300$')
plt.xlabel('local assortativity')
plt.ylabel(r'$\lambda^{(2)}$')
plt.gca().ticklabel_format(axis='x',style='scientific',scilimits=(-1,1),useMathText=True)
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(list32, list31, label=r'$N=500$')
plt.xlabel('local assortativity')
plt.ylabel(r'$\lambda^{(2)}$')
plt.gca().ticklabel_format(axis='x',style='scientific',scilimits=(-1,1),useMathText=True)
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(list42, list41, label=r'$N=1000$')
plt.xlabel('local assortativity')
plt.ylabel(r'$\lambda^{(2)}$')
plt.gca().ticklabel_format(axis='x',style='scientific',scilimits=(-1,1),useMathText=True)
plt.legend()

plt.savefig('assortativey_statis.pdf')
plt.show()