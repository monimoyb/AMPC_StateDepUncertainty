# A Semi-Definite Programming Approach to Robust Adaptive MPC under State Dependent Uncertainty
Adaptive MPC with convex programming: Handling state dependent additive uncertainty. These set of codes replicate the results of the paper https://ieeexplore.ieee.org/document/9143777. 

A computationally efficient SDP based algorithm is presented for robust MPC design under state dependent additive uncertainties. The additive uncertainty is assumed Lipschitz, with a known Lipschitz constant. The graph of this uncertainty is learned from system's historical closed-loop data, as shown in the following figure. 

<img width="1233" alt="d" src="https://user-images.githubusercontent.com/12418616/114331087-1af89e00-9af8-11eb-8269-f1420977ce8a.png">

The red dots are the collected closed-loop data points, each of which refines the graph envelope estimate constructed using Lipschitz property. As a result of such uncertainty learning and adaptation, the conservatism in control design lowers. This is seen in the figure below.

<img width="1066" alt="xn" src="https://user-images.githubusercontent.com/12418616/114331365-b12cc400-9af8-11eb-8c4f-63ea771b822d.png">

Starting from an empty terminal state with over-approximated uncertainty bounds, a non empty set is obtained, as the uncertainty envelope is refined with each collected data point.   
