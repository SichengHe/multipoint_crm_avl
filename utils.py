from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import copy

def plotRod(xyz1, xyz2, r, color, ax):

    '''
        Plot a rod with end node coordinates "xyz1" and "xyz2",
        radius "r", with color "color" in the frame of "ax".
    '''

    # plot one rod
    dphi = np.pi/3 # use 6 triangular pyramid
    N = int(np.pi*2 / dphi)

    # get unit rod dir vec, v: (dx, dy, dz)
    x1 = xyz1[0]
    y1 = xyz1[1]
    z1 = xyz1[2]

    x2 = xyz2[0]
    y2 = xyz2[1]
    z2 = xyz2[2]

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    dl = np.sqrt(dx**2 + dy**2 + dz**2)

    dx = dx/dl
    dy = dy/dl
    dz = dz/dl

    # perpendicular to v pick a line v_list[0]: (dx1, dy1, dz1)
    v_list = []
    if dz !=0:
        dx1 = 0
        dy1 = np.sqrt(dz**2/(dz**2 + dy**2))
        dz1 = -dy*dy1/dz

        v_list.append(np.array([dx1, dy1, dz1]))
    else:
        print('ERROR: Divide by zero -- dz=0!')
        exit()
        

    # Euler-Rodrigues formula
    a = np.cos(dphi/2)
    b = dx * np.sin(dphi/2)
    c = dy * np.sin(dphi/2)
    d = dz * np.sin(dphi/2)

    ER_mat = np.zeros((3, 3))

    ER_mat[0, 0] = a**2 + b**2 - c**2 - d**2
    ER_mat[0, 1] = 2*(b*c - a*d)
    ER_mat[0, 2] = 2*(b*d + a*c)
    ER_mat[1, 0] = 2*(b*c + a*d)
    ER_mat[1, 1] = a**2 + c**2 - b**2 -d**2
    ER_mat[1, 2] = 2*(c*d - a*b)
    ER_mat[2, 0] = 2*(b*d - a*c)
    ER_mat[2, 1] = 2*(c*d + a*b)
    ER_mat[2, 2] = a**2 + d**2 - b**2 - c**2

    # generate the vecs
    for i in range(N - 1):

        v_loc = v_list[i]
        v_list.append(ER_mat.dot(v_loc))

    # for convinience ...
    v_list.append(copy.deepcopy(v_list[0]))

    # scale to r
    for i in range(len(v_list)):
        v_list[i] *= r

    # for i in range(len(v_list)):
    #     print(np.linalg.norm(v_list[i]))
 
    # plot
    # disk 1:
    for i in range(N):

        tri_x1 = x1
        tri_y1 = y1
        tri_z1 = z1

        tri_x2 = tri_x1 + v_list[i][0]
        tri_y2 = tri_y1 + v_list[i][1]
        tri_z2 = tri_z1 + v_list[i][2]

        tri_x3 = tri_x1 + v_list[i + 1][0]
        tri_y3 = tri_y1 + v_list[i + 1][1]
        tri_z3 = tri_z1 + v_list[i + 1][2]

        tri_x = np.array([tri_x1, tri_x2, tri_x3])
        tri_y = np.array([tri_y1, tri_y2, tri_y3])
        tri_z = np.array([tri_z1, tri_z2, tri_z3])

        verts = np.zeros((3, 3))
        verts[:, 0] = tri_x[:]
        verts[:, 1] = tri_y[:]
        verts[:, 2] = tri_z[:]

        poly3d = Poly3DCollection(verts)
        poly3d.set_facecolor(color)
        ax.add_collection3d(poly3d)

    # disk 2:
    for i in range(N):

        tri_x1 = x2
        tri_y1 = y2
        tri_z1 = z2

        tri_x2 = tri_x1 + v_list[i][0]
        tri_y2 = tri_y1 + v_list[i][1]
        tri_z2 = tri_z1 + v_list[i][2]

        tri_x3 = tri_x1 + v_list[i + 1][0]
        tri_y3 = tri_y1 + v_list[i + 1][1]
        tri_z3 = tri_z1 + v_list[i + 1][2]

        tri_x = np.array([tri_x1, tri_x2, tri_x3])
        tri_y = np.array([tri_y1, tri_y2, tri_y3])
        tri_z = np.array([tri_z1, tri_z2, tri_z3])

        verts = np.zeros((3, 3))
        verts[:, 0] = tri_x[:]
        verts[:, 1] = tri_y[:]
        verts[:, 2] = tri_z[:]

        poly3d = Poly3DCollection(verts)
        poly3d.set_facecolor(color)
        ax.add_collection3d(poly3d)

    # side
    for i in range(N):

        # # from disk 1
        # quad_x1 = x1 + v_list[i][0]
        # quad_y1 = y1 + v_list[i][1]
        # quad_z1 = z1 + v_list[i][2]

        # quad_x2 = x1 + v_list[i + 1][0]
        # quad_y2 = y1 + v_list[i + 1][1]
        # quad_z2 = z1 + v_list[i + 1][2]

        # # from disk 2
        # quad_x3 = x2 + v_list[i][0]
        # quad_y3 = y2 + v_list[i][1]
        # quad_z3 = z2 + v_list[i][2]

        # quad_x4 = x2 + v_list[i + 1][0]
        # quad_y4 = y2 + v_list[i + 1][1]
        # quad_z4 = z2 + v_list[i + 1][2]

        # quad_x = [quad_x1, quad_x2, quad_x4, quad_x3]
        # quad_y = [quad_y1, quad_y2, quad_y4, quad_y3]
        # quad_z = [quad_z1, quad_z2, quad_z4, quad_z3]

        # verts = np.zeros((4, 3))
        # verts[:, 0] = quad_x[:]
        # verts[:, 1] = quad_y[:]
        # verts[:, 2] = quad_z[:]


        # poly3d = Poly3DCollection(verts)
        # poly3d.set_facecolor(color)
        # ax.add_collection3d(poly3d)

        # from disk 1
        quad_x1 = x1 + v_list[i][0]
        quad_y1 = y1 + v_list[i][1]
        quad_z1 = z1 + v_list[i][2]

        quad_x2 = x1 + v_list[i + 1][0]
        quad_y2 = y1 + v_list[i + 1][1]
        quad_z2 = z1 + v_list[i + 1][2]

        # from disk 2
        quad_x3 = x2 + v_list[i][0]
        quad_y3 = y2 + v_list[i][1]
        quad_z3 = z2 + v_list[i][2]

        quad_x4 = x2 + v_list[i + 1][0]
        quad_y4 = y2 + v_list[i + 1][1]
        quad_z4 = z2 + v_list[i + 1][2]


        # break the quad into two triangulars
        # tri 1
        tri1_x = [quad_x1, quad_x2, quad_x4]
        tri1_y = [quad_y1, quad_y2, quad_y4]
        tri1_z = [quad_z1, quad_z2, quad_z4]

        verts1 = np.zeros((3, 3))
        verts1[:, 0] = tri1_x[:]
        verts1[:, 1] = tri1_y[:]
        verts1[:, 2] = tri1_z[:]

        poly3d = Poly3DCollection(verts1)
        poly3d.set_facecolor(color)
        ax.add_collection3d(poly3d)

        # tri 2
        tri2_x = [quad_x4, quad_x3, quad_x1]
        tri2_y = [quad_y4, quad_y3, quad_y1]
        tri2_z = [quad_z4, quad_z3, quad_z1]
 
        verts2 = np.zeros((3, 3))
        verts2[:, 0] = tri2_x[:]
        verts2[:, 1] = tri2_y[:]
        verts2[:, 2] = tri2_z[:]

        poly3d = Poly3DCollection(verts2)
        poly3d.set_facecolor(color)
        ax.add_collection3d(poly3d)






def preprocess(filename_node, filename_elem, \
filename_sol, filename_con, \
ind_beg_elem, ind_end_elem, ind_beg_node, ind_end_node, \
length_scale, yield_stress):

    '''
        Extract areas, stress, buckling state, and deformed truss coordinates.
    '''

    nodes = np.loadtxt(filename_node)
    nodes = nodes/length_scale
    elems = np.loadtxt(filename_elem).astype(int)
    nodes_d = copy.deepcopy(nodes)

    N_elem, __ = elems.shape

    areas = np.zeros(N_elem)
    stress = np.zeros(N_elem)
    bucklings = np.zeros(N_elem)

    i = 0
    j = 0
    with open(filename_sol) as f:
        for line in f:
            if (ind_beg_elem<=j and j<=ind_end_elem):
                line_split = line.split()

                areas_loc = float(line_split[1])/length_scale**2 # m^2
                stress_loc = float(line_split[4]) 
                buckling_loc = float(line_split[5]) 

                if buckling_loc>0:
                    buckling_loc = 0.0
                else:
                    buckling_loc = -buckling_loc

                areas[i] = areas_loc
                stress[i] = stress_loc
                bucklings[i] =  buckling_loc

                i += 1
                
            j += 1

    i = 0
    j = 0

    dys = np.zeros(ind_end_node - ind_beg_node + 1)
    with open(filename_sol) as f:
        for line in f:
            if (ind_beg_node<=j and j<=ind_end_node):
                line_split = line.split()

                dys[i] = float(line_split[1])

                i += 1
                
            j += 1

    cons = []
    with open(filename_con) as f:
        for line in f:
                line_split = line.split()

                cons.append(int(float(line_split[0])))


    dys /= length_scale
    dys = dys.reshape((-1, 3))
    nodes_d = nodes + dys


    abs_stress = abs(stress)
    rel_abs_stress = abs_stress/yield_stress
    rel_stress = -stress/yield_stress

    return nodes, nodes_d, elems, cons, areas, rel_stress, bucklings