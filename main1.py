import numpy as np
import DFN as dfn

# Settiamo i parametri del costruttore del DFN

N = 4
Xmin = 0
Xmax = 3
Ymin = 0
Ymax = 3
Zmin = 0
Zmax = 5
alpha_pl = 2
radius_l = 1.5
radius_u = 2.5
k = 3
mode_vector = np.array([[0.], [0.], [1.]])
fixed_n_edges = 0

network = dfn.DiscreteFractureNetwork(N, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                                      alpha_pl, radius_l, radius_u, k, mode_vector, fixed_n_edges)
network.visual3D('1')
network.scrittura1()
network.scrittura2()
print(network.poss_intersezioni)
print(network.intersezioni)
print(network.frac_traces)
print(network.traces)
n2 = dfn.DiscreteFractureNetwork(2, Xmin, Xmax, Ymin, Ymax, Zmin, Zmax,
                                      alpha_pl, radius_l, radius_u, k, mode_vector, fixed_n_edges)
network.aggiungi(n2.fractures)
network.scrittura1('file3.txt')
network.scrittura2('file4.txt')
network.visual3D('2')
print(network.poss_intersezioni)
print(network.intersezioni)
print(network.frac_traces)
print(network.traces)