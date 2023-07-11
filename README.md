<h1 align="center" style="display: block; font-size: 2.5em; font-weight: bold; margin-block-start: 1em; margin-block-end: 1em;">
<br><br><strong> Time-dependent master equations and quantum trajectories</strong>
</h1>


[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Latest release](https://img.shields.io/github/v/release/kwyip/Adiabatic-master-equation-and-quantum-trajectories)](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/releases/tag/v2.0.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2022a-BLUE.svg)](https://www.mathworks.com/login?uri=%2Fdownloads%2Fprerelease)

---

## Table of contents[![](./docs/img/pin.svg)](#table-of-contents)
1. [Motivation](#motivation)
2. [Code description](#codedescription)
3. [Tutorial](#tutorial)
    - [UCL 4 qubit gadgets](#ucl)
    - [Single qubit](#single)
4. [Reference](#reference)
---

## Motivation[![](./docs/img/pin.svg)](#motivation)
Simulation codes of time-dependent master equation and quantum trajectories.

## Code description <a name="codedescription"></a>
These are the codes of adiabatic master equation (`ame.m`) and quantum trajectories (`aqt.m`) used in this paper [[1]](#1).

The Hamiltonian is a 8-qubit ferromagnetic Ising spin chain in a transverse field, with annealing schedule attached as `DW1_parameters.txt`. This relatively complex problem illustrates a lot of subtle points in the implementation of *ame* and *aqt*, and can be easily reduced to smaller and simpler problems. 


![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/8-qubit_chain.png)



The master equation is:

![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame1.png)
![alt text](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/ame2.png)









Both `ame.m` and `aqt.m` use sparse matrix for computation,

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/sparsem%20(1).png)

and perform basis rotation for every ode step.

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/rotation%20(1).png)

`aqt.m` uses parallel computation, with sorting of jump operators,

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/sortingjump%20(1).png)

and backtracking for error control.

![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/backtracking%20(1).png)


## Tutorial: UCL's 4-qubit gadgets and a single qubit <a name="tutorial"></a>
### UCL's 4-qubit gadgets <a name="ucl"></a>
The folder [tutorial_4_qubit_gadget](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/tree/master/tutorial_4_qubit_gadget) contains a tutorial on how to run the codes of a time-dependent case of the 4-qubit gadgets below. 
![alt text](https://github.com/USCqserver/Adiabatic-master-equation-and-quantum-trajectories/blob/master/images/4-qubit.png)
### Single qubit <a name="single"></a>
A [onequbit_demo](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/tree/master/tutorial_4_qubit_gadget/onequbit_demo) instruction folder was made for how to run codes on various time-dependent single-qubit Hamiltonian simulation. The pdf file [kawacode.pdf](https://github.com/kwyip/Adiabatic-master-equation-and-quantum-trajectories/blob/master/tutorial_4_qubit_gadget/kawacode.pdf) contains a presentation on excecutions.

## Reference <a name="reference"></a>
<a id="1">[1]</a> 
Yip, K.W., Albash, T. and Lidar, D.A., 2018. [Quantum trajectories for time-dependent adiabatic master equations](https://arxiv.org/pdf/1710.03431.pdf). Physical Review A, 97(2), p.022116.
