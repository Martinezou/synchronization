"""the difference between original eigenvalue and eigenvalue after link removal"""
from numpy import *
import scipy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def main(N, K,S):
    G1 = nx.barabasi_albert_graph(N, K, seed=S)
    a = G1.degree()
    b = sorted(a, key=lambda t:t[1])
    lap_matrix = nx.laplacian_matrix(G1).todense()
    eigval, eigvec = scipy.linalg.eig(lap_matrix)
    # print(eigval)
    d = eigval.sort()
    eig_sec1 = float(eigval[1])
    print(eig_sec1)
    degree_sequence = [d for n, d in G1.degree()]
    #print(b)

    def nodes_min_degree():
        '''find the nodes with minimum degree of the network
        in scale free network there exists a lot of nodes with minimum degrees'''
        min_degree = min(degree_sequence) #minimum degree of the network
        #print(min_degree)
        min_node1 = []
        for i in range(len(b)):
            if b[i][1] == min_degree:
                min_node1.append(b[i][0])
        return min_node1
    min_node_list = nodes_min_degree()
    #print(min_node_list)

    def remove_edge(node):
        '''remove edge between min_degree node and its max_degree neighbour'''
        G2 = nx.barabasi_albert_graph(N, K, seed=S)
        a = G2.degree()
        b = sorted(a, key=lambda t: t[1])
        #print(b)
        neighbour = list(G1.adj[node])
        list1 = []
        for i in range(len(neighbour)):
            list1.append(G1.degree(neighbour[i]))
        max_neigh = list1.index(max(list1))
        G2.remove_edge(node, neighbour[max_neigh])
        return G2


    def attribution():
        removal_eigen = []
        for i in range(len(min_node_list)):
            G = remove_edge(min_node_list[i])
            lap_matrix = nx.laplacian_matrix(G).todense()
            eigval, eigvec = scipy.linalg.eig(lap_matrix)
            d = eigval.sort()
            eig_sec1 = float(eigval[1])
            removal_eigen.append(eig_sec1)
        print(removal_eigen)
        #print(len(removal_eigen))
        return removal_eigen
    remove_eig = attribution()
    return eig_sec1, remove_eig
#a = main(500, 20, 1)


#N=1000, K=20
fig1 = plt.figure(figsize=(10,8))
eigen, perturbation_eigen = main(1000, 20, 2)
origional_eigen = [eigen for i in range(len(perturbation_eigen))]
difference_eigen = [-origional_eigen[0]+i for i in perturbation_eigen]
difference_eigen_percent = [i/origional_eigen[0] for i in difference_eigen]
fig=plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax1= fig.add_axes([left, bottom, width, height])
plt.hist(difference_eigen_percent, cumulative=True, bins=50, histtype='step', density=True)
plt.xlabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.ylabel('Probability')

left, bottom, width, height = 0.2, 0.6, 0.25, 0.25
ax2 = fig.add_axes([left, bottom, width, height])
ax2.ticklabel_format(style='sci', scilimits=(-3,0), axis='y')
plt.boxplot([difference_eigen_percent], showfliers=False, labels=[r'$N=1000, k_{\rm min}=20$'])
plt.ylabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.savefig('perturbation probability1.pdf')
plt.show()

#N=1000, K=15
fig1 = plt.figure(figsize=(10,8))
eigen, perturbation_eigen = main(1000, 15, 1)
origional_eigen = [eigen for i in range(len(perturbation_eigen))]
difference_eigen = [-origional_eigen[0]+i for i in perturbation_eigen]
difference_eigen_percent = [i/origional_eigen[0] for i in difference_eigen]
fig=plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax1= fig.add_axes([left, bottom, width, height])
plt.hist(difference_eigen_percent, cumulative=True, bins=50, histtype='step', density=True)
plt.xlabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.ylabel('Probability')

left, bottom, width, height = 0.2, 0.6, 0.25, 0.25
ax2 = fig.add_axes([left, bottom, width, height])
ax2.ticklabel_format(style='sci', scilimits=(-3,0), axis='y')
plt.boxplot([difference_eigen_percent], showfliers=False, labels=[r'$N=1000, k_{\rm min}=20$'])
plt.ylabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.savefig('perturbation probability2.pdf')
plt.show()

#N=500, K=20
fig1 = plt.figure(figsize=(10,8))
eigen, perturbation_eigen = main(500, 20, 1)
origional_eigen = [eigen for i in range(len(perturbation_eigen))]
difference_eigen = [-origional_eigen[0]+i for i in perturbation_eigen]
difference_eigen_percent = [i/origional_eigen[0] for i in difference_eigen]
fig=plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax1= fig.add_axes([left, bottom, width, height])
plt.hist(difference_eigen_percent, cumulative=True, bins=50, histtype='step', density=True)
plt.xlabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.ylabel('Probability')

left, bottom, width, height = 0.25, 0.6, 0.25, 0.25
ax2 = fig.add_axes([left, bottom, width, height])
ax2.ticklabel_format(style='sci', scilimits=(-3,0), axis='y')
plt.boxplot([difference_eigen_percent], showfliers=False, labels=[r'$N=1000, k_{\rm min}=20$'])
plt.ylabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.savefig('perturbation probability3.pdf')
plt.show()

#N=500, K=15
fig1 = plt.figure(figsize=(10,8))
eigen, perturbation_eigen = main(500, 15, 1)
origional_eigen = [eigen for i in range(len(perturbation_eigen))]
difference_eigen = [-origional_eigen[0]+i for i in perturbation_eigen]
difference_eigen_percent = [i/origional_eigen[0] for i in difference_eigen]
fig=plt.figure()
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax1= fig.add_axes([left, bottom, width, height])
plt.hist(difference_eigen_percent, cumulative=True, bins=50, histtype='step', density=True)
plt.xlabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.ylabel('Probability')

left, bottom, width, height = 0.25, 0.6, 0.25, 0.25
ax2 = fig.add_axes([left, bottom, width, height])
ax2.ticklabel_format(style='sci', scilimits=(-3,0), axis='y')
plt.boxplot([difference_eigen_percent], showfliers=False, labels=[r'$N=1000, k_{\rm min}=20$'])
plt.ylabel(r'$\Delta \lambda^{(2)}/\lambda^{(2)}$')
plt.savefig('perturbation probability3.pdf')
plt.show()