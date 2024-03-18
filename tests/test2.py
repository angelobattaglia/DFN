import numpy as np
import sys
sys.path.append('../src')
from DFN import DiscreteFractureNetwork as dfn

# Setting the parameters for out network

N = 4
Xmin = 0
Xmax = 5
Ymin = 0
Ymax = 5
Zmin = 0
Zmax = 5
alpha_pl = 2
radius_l = 1.5
radius_u = 2.5
k = 3
mode_vector = np.array([[0.], [0.], [1.]])
fixed_n_edges = 0

# Let's create a network

## The Following two lines of code will work if used in the same directory of DFN.py
#network = dfn.DiscreteFractureNetwork(N, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                                      #alpha_pl, radius_l, radius_u, k, mode_vector, fixed_n_edges)

network = dfn(N, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                                        alpha_pl, radius_l, radius_u, k, mode_vector, fixed_n_edges)

network.visual3D('1')
network.scrittura1()
network.scrittura2()
print(network.poss_intersezioni)
print(network.intersezioni)
print(network.frac_traces)
print(network.traces)
network.genfrac(2)
network.scrittura1('file3.txt')
network.scrittura2('file4.txt')
network.visual3D('2')
print(network.poss_intersezioni)
print(network.intersezioni)
print(network.frac_traces)
print(network.traces)