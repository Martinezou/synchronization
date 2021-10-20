"""Estimation of lambdaN"""

import matplotlib.pyplot as plt
import networkx as nx
from numpy import *
import scipy


def maximum_eigen():
    list1 = []
    list2 = []
    for i in range(10, 31):
        for j in range(0, 10):
            G1 = nx.barabasi_albert_graph(500, i, seed=j)
            degree_sequence = [d for n, d in G1.degree()]
            maxeigen = max(degree_sequence)
            lap_matrix = nx.laplacian_matrix(G1).todense()
            eigval, eigvec = scipy.linalg.eig(lap_matrix)
            d = eigval.sort()
            eig_max = eigval[-1]
            list1.append(maxeigen)
            list2.append(float(eig_max))
    #print(list1)
    #print(list2)
    return list1, list2
estimation, maxeigen = maximum_eigen()
print(estimation)
print(maxeigen)
plt.scatter(maxeigen, estimation, c='b', s=6., label=r'Estimation of $\lambda^{(N)}$')
plt.plot(maxeigen, maxeigen, c='y', label = r'$\lambda^{(N)}$')
plt.xlabel(r'$\lambda^{(N)}$')
plt.ylabel(r'estimate $\lambda^{(N)}$')
plt.savefig('estimation_accuracy_maxeigen.pdf')
plt.legend()
plt.show()