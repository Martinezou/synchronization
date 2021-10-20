"""Coupled rossler system with different network topology;
The BA model generated here will be used in networksynchronization.py"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import scipy

G = nx.barabasi_albert_graph(100, 10, seed=2)  # BA network

def topology_matrix(GG):
    lap_mat = nx.laplacian_matrix(GG).todense()
    A2 = nx.to_numpy_array(GG)  # adjacent matrix of network
    eigval, eigvec = scipy.linalg.eig(lap_mat)
    eigval.sort()
    #print(eigval)
    eig_sec = round(eigval[-1], 4)
    eig_min = round(eigval[1], 4)
    print('max eigen:' + str(eig_sec))
    print('min eigen:' + str(eig_min))
    return A2
AA = topology_matrix(G)  #original network


def link_xmatch(GG, node1, node2, node3, node4):
    '''1-2;3-4
        1-3;2-4'''
    GG.remove_edge(node1, node2)
    GG.add_edge(node1, node3)
    GG.remove_edge(node3, node4)
    GG.add_edge(node2, node4)
    return GG


def key_link_reconection(GG, n1, n2):
    neighbour1 = list(GG.adj[n1])
    #print(neighbour1)
    neighbour2 = list(GG.adj[n2])
    a = 0
    if n1 not in neighbour2:
        i = 0
        while a <1:
            neighbour3 = list(GG.adj[neighbour2[i]])
            neighbour3.extend(neighbour2)
            neighbour4 = list(set(neighbour3))
            if neighbour1[i] not in neighbour4:
                link_xmatch(GG, n1, neighbour1[i], n2, neighbour2[i])
                a +=1
            else:
                i+=1
                if i>= len(neighbour1):
                    break
    return GG


def nodes_min_degree(GG):
    '''find the nodes with minimum degree of the network
    in scale free network there exists a lot of nodes with minimum degrees'''
    a = GG.degree()
    b = sorted(a, key=lambda t: t[1])
    degree_sequence = [d for n, d in GG.degree()]
    #print(degree_sequence)
    min_degree = min(degree_sequence) #minimum degree of the network
    #print(min_degree)
    min_node1 = []
    for i in range(len(b)):
        if b[i][1] == min_degree:
            min_node1.append(b[i][0])
    #print('nodes with min degree:'+str(min_node1))
    return min_node1


def reconnection1(GG, a, b):
    min_node = nodes_min_degree(GG)
    for i in range(a):
        for j in range(a, b):
            GG = key_link_reconection(GG, min_node[i], min_node[j])
    r = nx.degree_assortativity_coefficient(GG)
    #print(r)
    return GG,  r


for i in range(7): # 7 paper used
    G3 = reconnection1(G, i, 2 * i)[0]

A = nx.to_numpy_array(G3)   #network after connection
topology_matrix(G3)