print(elem_new.shape)

import matplotlib
import matplotlib.pyplot as plt

def plotWing(node, node_new, elem_new, node_2_node_new):
    
    node_xy = node[:, 0:2]
    node_new_xy = node_new[:, 0:2]

    num_elem, _ = elem_new.shape
    num_node, _ = node_xy.shape
    num_node_new, _ = node_new_xy.shape


    for i in xrange(num_elem):
        ind1 = elem_new[i, 0]
        ind2 = elem_new[i, 1]

        x1 = node_new_xy[ind1, 0]
        y1 = node_new_xy[ind1, 1]
        x2 = node_new_xy[ind2, 0]
        y2 = node_new_xy[ind2, 1]

        plt.plot([x1, x2], [y1, y2], '-b')

    alpha_array = np.linspace(0.1, 1, num_node_new)
    for i in xrange(num_node):
        alpha_loc = alpha_array[node_2_node_new[i]]
        plt.plot(node_xy[i, 0], node_xy[i, 1], 'bo', alpha=alpha_loc)


    plt.show()
    
plotWing(node, node_new, elem_new, node_2_node_new)
    