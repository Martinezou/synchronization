# global synchronization

# About the Project

This project is aim to show the synchronizability of coupling rossler system. We use non-generate perturbation theory to estimate $$\lambda_N$$ and $\lambda_2$, then compare the estimation with the real value. And we have shown the synchronization process of couping rossler system in different network topology.

1.estimation_method.py 
This is mainly used to generate data for estimation of $$\lambda_2$$. The running time of this code is very long. So we have written the data to the csv document perturbation_estimation.csv through this code.

2.lambda2_estimation.py
This code is used to verify the estimation accuracy of $$\lambda_2$$ by our proposed purterbation theory. The data used here is generated from estimation_method.py.

3.lambdaN_estimation.py
This code is used to verify the estimation accuracy of $$\lambda_N$$.

4.remove_accuracy.py
This code is used to show the effect of link removal on $$\lambda_2$$. 

5.lyapunov.py
This code is used to calculate maximum lyapunov exponent of the coupling dynamic system. When the maximum lyapunov exponent is less than 0, nodes tend to synchronize.

6.network_topology.py
Generate different network model for coupling dynamic system. The network model will be used in coupled_rossler_system.py.

7.coupled_rossler_system.py
This code shows the synchronization process of coupling rossler system in different network topology.

8.assortativity_statis.py
Calculte assortativity of the network and show its relationship with $$\lambda_2$$.

# License
Distributed under the MIT License. 

