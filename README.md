# Adiabatic-master-equation-and-quantum-trajectories
Code of adiabatic master equation and quantum trajectories used in paper
These are the codes of adiabatic master equation and quantum trajectories used in paper:
https://journals.aps.org/pra/abstract/10.1103/PhysRevA.97.022116

The Hamiltonian is a 8-qubit ferromagnetic Ising spin chain in a transverse field, with annealing schedule attached as 'DW1_parameters.txt'. This relatively complex problem illustrates a lot of subtle points in the implementation of ame and aqt, and can be easily reduced to smaller and simpler problem. 


![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/8-qubit_chain.png)






The master equation is:

![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame1.png)
![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame2.png)









Both ame.m and aqt.m use sparse matrix for computation.
