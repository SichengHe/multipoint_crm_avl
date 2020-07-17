import numpy as np
import copy

# load data
node_file_name = "./Wings/W_315/data_nodes.dat"
elem_file_name = "./Wings/W_315/data_elems.dat"
con_file_name = "./Wings/W_315/data_constraints.dat"
force_file_name = "./OUT_STRUCT/Wings/W_315/data_forces_Cl-0.244Ma0.64.dat" # -1g
# force_file_name = "./OUT_STRUCT/Wings/W_315/data_forces_Cl0.611Ma0.64.dat" # 2.5g
# force_file_name = "./OUT_STRUCT/Wings/W_315/data_forces_Cl0.403Ma0.86.dat" # gust
R_file_name = "./Wings/W_315/R.dat"

node = np.loadtxt(node_file_name)
elem = np.loadtxt(elem_file_name)
elem = elem.astype(int)
con = np.loadtxt(con_file_name)
con = con.astype(int)
force = np.loadtxt(force_file_name)
R = np.loadtxt(R_file_name) # w/o constraint nodes

force_wo_con = np.delete(force, con, 0)
node_wo_con = np.delete(node, con, 0)

# construct K
E = 68947.57293 * 1e6
h = 0.4 # thickest bar
A = h**2

N_elem = elem.shape[0]
EI_l = np.zeros(N_elem)
for i in range(N_elem):

    ind1 = elem[i, 0]
    ind2 = elem[i, 1]

    node1 = node[ind1, :]
    node2 = node[ind2, :]

    delta = np.linalg.norm(node2 - node1)

    EI_l[i] = E * A / delta

K = R.dot(np.diag(EI_l)).dot(np.transpose(R))
print("np.linalg.cond(K)", np.linalg.cond(K))
print("K", K)

# solve
u_wo_con = np.linalg.solve(K, force_wo_con.flatten()).reshape(-1, 3)
u = copy.deepcopy(u_wo_con)

for con_loc in con:
    u = np.insert(u, con_loc, 0.0, axis = 0)

nodes_d = node + u


# plot
x_center = (max(nodes_d[:, 0]) + min(nodes_d[:, 0])) / 2
x_radius = (max(nodes_d[:, 0]) - min(nodes_d[:, 0])) / 2
y_center = (max(nodes_d[:, 1]) + min(nodes_d[:, 1])) / 2
y_radius = (max(nodes_d[:, 1]) - min(nodes_d[:, 1])) / 2
z_center = (max(nodes_d[:, 2]) + min(nodes_d[:, 2])) / 2
z_radius = (max(nodes_d[:, 2]) - min(nodes_d[:, 2])) / 2

alpha = 1.3
x_lb = x_center - x_radius * alpha
x_ub = x_center + x_radius * alpha
y_lb = y_center - y_radius * alpha
y_ub = y_center + y_radius * alpha
z_lb = z_center - z_radius * alpha
z_ub = z_center + z_radius * alpha

x_lb = min(min(x_lb, y_lb), z_lb)
x_ub = max(max(x_ub, y_ub), z_ub)
y_lb = x_lb
z_lb = x_lb
y_ub = x_ub
z_ub = x_ub

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure(figsize=(16, 16))
ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
# ax.set_aspect('equal')
ax.pbaspect = [1,1,1]

ax.set_xlim([0, 30])
ax.set_ylim([0, 30])
ax.set_zlim([-15, 15])
ax.grid(False)
ax._axis3don = False

# jig shape
for i in range(N_elem):

    ind1 = elem[i, 0]
    ind2 = elem[i, 1]

    xyz1 = node[ind1, :]
    xyz2 = node[ind2, :]

    ax.plot([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], 'k')

# deformed
for i in range(N_elem):

    ind1 = elem[i, 0]
    ind2 = elem[i, 1]

    xyz1 = nodes_d[ind1, :]
    xyz2 = nodes_d[ind2, :]

    ax.plot([xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]], 'b')


ax.plot([x_lb, x_lb, x_lb, x_lb, x_ub, x_ub, x_ub, x_ub],
[y_lb, y_lb, y_ub, y_ub, y_ub, y_ub, y_lb, y_lb],
[z_lb, z_ub, z_ub, z_lb, z_lb, z_ub, z_ub, z_lb], alpha = 0)

ax.view_init(elev=0, azim=0)

# plt.show()
fig.savefig("n1g.pdf", bbox_inches='tight') 