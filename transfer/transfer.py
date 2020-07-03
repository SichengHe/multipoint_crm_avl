import numpy as np


# # load and node coord
# FE = np.loadtxt("FE.txt")
# node = node = FE[:, 0:3]
# load = FE[:, 3]


# # node coord to be interpolated with
# node_new = np.loadtxt("data_nodes.dat")



def getAssociation(node, node_new):
    
    # associate the node (we know load) with node_new (we dont know load)
    # For each "node" we pick the closest "node_new" and the load is transfered there
    
    num_node, _ = node.shape
    num_node_new, _ = node_new.shape

    node_new_2_node = []
    for i in range(num_node_new):
        node_new_2_node.append([])

    node_2_node_new = []

    for i in range(num_node):
        d = []
        for j in range(num_node_new):
            d_loc = node_new[j, :] - node[i, :]
            d_loc = np.sqrt(d_loc.dot(d_loc))

            d.append(d_loc)

        ind = d.index(min(d))

        node_new_2_node[ind].append(i)
        node_2_node_new.append(ind)

    return node_new_2_node, node_2_node_new



def getLoad_new(node, load, node_new):
    
    # INPUT: 
    # node coord, load - "node, load" from the old config to be interpolated with
    # node coord - "node_new" from the new config to be evaluated at

    # OUTPUT:
    # load for the new config

    # ASSUMPTION:
    # load is equally distributed on lower and upper point in the new config
    # assume that the "node_new" data has specific order:
    # odd:  [x, y, z+]
    # even: [x, y, z-]
    # -- they share the same x, y pair


    node_xy = node[:, 0:2]
    node_new_half_xy = node_new[:, 0:2][0::2]


    node_new_half_xy_2_node, node_2_node_new_half_xy = getAssociation(node_xy, node_new_half_xy)
    load_new_half_xy = np.zeros(node_new_half_xy.shape[0])



    for i in range(load_new_half_xy.shape[0]):
        
        child = node_new_half_xy_2_node[i]

        for j in range(len(child)):
            
            child_loc = child[j]
            load_new_half_xy[i] += load[child_loc]


    load_new = np.zeros((node_new.shape[0], 3))

    for i in range(load_new_half_xy.shape[0]):
        load_new[2*i, 2] = load_new_half_xy[i]/2
        load_new[2*i + 1, 2] = load_new_half_xy[i]/2

    return load_new


