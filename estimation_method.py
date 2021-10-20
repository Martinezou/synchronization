'''perturbation estimation method
estimate lambda2'''
from numpy import *
import scipy
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time
import pandas as pd

def main(N, K,S):
    G1 = nx.barabasi_albert_graph(N, K, seed=S)
    a = G1.degree()
    b = sorted(a, key=lambda t:t[1])
    degree_sequence = [d for n, d in G1.degree()]
    #print(b)


    def attribution():
        '''attribution of the graph'''
        c = mean(list(G1.degree()))
        min_degree = degree_sequence[0]
        #print(min_degree)
        #print('average degree:' + str(c))
        d = K * (1 - 2 / sqrt(c))
        #print('estimate by mean degree:' + str(d))

    attribution()

    def attribution_node(node):
        '''attribution of a node'''
        neighbours = list(G1.adj[node])
        neighbours_degree = []
        for i in range(len(neighbours)):
            neighbours_degree.append(G1.degree(neighbours[i]))
        #print(neighbours_degree)
        #print(sorted(neighbours_degree))
        #print('neighbours:'+str(neighbours))


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

    nodes_min_degree()


    def average_degree_neighbour():
        '''calculate average degree of nodes with minimum degree'''
        list1 = nodes_min_degree()
        averdegree = []
        vardegree = []
        for i in range(len(list1)):
            neighbours = list(G1.adj[list1[i]])
            list2 = []
            for j in range(len(neighbours)):
                list2.append(G1.degree(neighbours[j])-K)
            averdegree.append(mean(list2))
            vardegree.append(np.var(list2))
        #print('average degree of neighbours'+str(averdegree))
        averdegree_dict1 = dict(zip(list1, averdegree))
        #vardegree_dict1 = dict(zip(list1, vardegree))
        averdegree_dict2 = sorted(averdegree_dict1.items(), key=lambda t:t[1])
        #vardegree_dict2 = sorted(vardegree_dict1.items(), key=lambda t:t[1])
        #print('average degree of neighbours'+str(averdegree_dict2))
        #print('var degree of neighbours'+str(vardegree_dict2))

        return averdegree_dict2
    #averdegree_dict = average_degree_neighbour()

    def estimate_eig_sec(a):
        lap_matrix = nx.laplacian_matrix(a).todense()
        eigval, eigvec = scipy.linalg.eig(lap_matrix)
        d = eigval.sort()
        eig_sec = float(eigval[1])
        #print('laplacian eigen:'+str(eig_sec))
        return eig_sec
    original_eigen = estimate_eig_sec(G1)


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


    def perturbation_estimation():
        list1 = nodes_min_degree()
        #print(list1)
        list3 = []
        list5 = []
        eig_list=[] # second eigenvalue of perturbation
        for i in range(len(list1)):
            G = remove_edge(list1[i])
            lap_matrix = nx.laplacian_matrix(G).todense()
            eigval, eigvec = scipy.linalg.eig(lap_matrix)
            d = eigval.sort()
            eig_sec1 = float(eigval[1])
        #print('eigenvalue of perturbation matrix:' + str(eig_sec1))
            eig_list.append(eig_sec1)
            degree_sequence = [d for n, d in G.degree()]
            min_degree = min(degree_sequence)
            #print(min_degree)
            neighbour = list(G.adj[list1[i]])
            #print(neighbour)
            list2 = []
            for i in range(len(neighbour)):
                a = 1 / (min_degree - G.degree(neighbour[i]))
                list2.append(a)
            d1 = sum(list2)
            list3.append(min_degree)
            list5.append(min_degree+d1)
        dict1 = dict(zip(list1,list3))
        dict3 = dict(zip(list1,list5))
        dict2 = sorted(dict1.items(), key=lambda t:t[1])
        dict4 = sorted(dict3.items(), key=lambda t:t[1])
        return dict2, dict4
    purterb_eig_list = perturbation_estimation()[0]
    first_est = purterb_eig_list[0][1]
    print('real value:'+str(original_eigen))
    print('first estimation:'+str(first_est))
    est_eig_list = perturbation_estimation()[1]
    second_est = est_eig_list[0][1]
    print('second estimation:'+str(second_est))
    return original_eigen, first_est, second_est

#main(500, 12, 2)
list1 = []
list2 = []
list3 = []
start = time.time()
for i in range(10, 31):
    for j in range(0, 10):
        a = main(500, i, j)[0]
        if a > 9:
            list1.append(a)
            list2.append(main(500, i, j)[1])
            list3.append(main(500, i, j)[2])
print(list1)
print(list2)
print(list3)

dict = {'real value': list1, 'first order estimation': list2, 'second order estimation': list3}
df = pd.DataFrame(dict)
df.to_csv('perturbation_estimation.csv')
print(df)
end = time.time()
print('running time:'+str(end-start))
