import numpy as np

area_file_name = "./solution/old_opt_sol/areas.dat"
node_file_name = "./Wings/W_315/data_nodes.dat"
elem_file_name = "./Wings/W_315/data_elems.dat"

areas = np.loadtxt(area_file_name)
nodes = np.loadtxt(node_file_name)
elems = np.loadtxt(elem_file_name).astype(int)

M = elems.shape[0]

L = np.zeros(M)

for i in range(M):

    ind1 = elems[i, 0]
    ind2 = elems[i, 1]
    
    node1 = nodes[ind1, :]
    node2 = nodes[ind2, :]

    L[i] = np.linalg.norm(node2 - node1)

volume = 0 
for i in range(M):

    dVolume = areas[i] * L[i]

    volume += dVolume

rho = 2767.0
mass = volume * rho

print("mass: ", mass, "kg")
