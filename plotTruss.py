# Sicheng He hschsc@umich.edu
# 01/24/2018
# truss result visualization 
# area/stress/buckling

from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import copy
import utils as utils

# =================
# get data
# =================
class trussPlotter():

    '''
        Plot the truss.
        It can plot the undeformed and deformed truss with 
        stress and buckling field. Stress has been normalized 
        between -1, 1 by sigma_Y and buckling is normalized 
        with (gamma * x). In addition, for buckling, we take

            max( - sigma / (gamma x), 0 ).
        
        Inputs:
            nodes (w con): nodes, float
            elems (w con): elements, integer
            cons: constrained nodes, integer
            areas: crossectional areas (in m^2), float
            disps: displacement (in m), float
            rel_stress: sigma / sigma_Y, float
            bucklings: sigma / (gamma A), float

    '''

    def __init__(self, nodes, elems, cons, areas, disps, rel_stress, bucklings):

        self.nodes = nodes
        self.elems = elems
        self.cons = cons
        self.areas = areas
        self.disps = disps
        self.rel_stress = rel_stress
        self.bucklings = bucklings

        self.nodes_d = nodes + disps

        self.N_elem = self.elems.shape[0]

        # generated the coordinates of a box later used 
        # to create equal space figure
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

        self.x_lb = min(min(x_lb, y_lb), z_lb)
        self.x_ub = max(max(x_ub, y_ub), z_ub)
        self.y_lb = self.x_lb
        self.z_lb = self.x_lb
        self.y_ub = self.x_ub
        self.z_ub = self.x_ub

    def plot_stress(self, elev, azim, job_name):

        coolwarm_cmap = matplotlib.cm.get_cmap('coolwarm')
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)

        # view 1
        fig = plt.figure(figsize=(16, 16))
        ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')
        ax.pbaspect = [1,1,1]

        ax.set_xlim([0, 30])
        ax.set_ylim([0, 30])
        ax.set_zlim([-15, 15])
        ax.grid(False)
        ax._axis3don = False

        # jig shape
        for i in range(self.N_elem):

            ind1 = elems[i, 0]
            ind2 = elems[i, 1]

            xyz1 = nodes[ind1, :]
            xyz2 = nodes[ind2, :]

            r_loc = np.sqrt(areas[i]/np.pi)

            color = [192.0/255, 192.0/255, 192.0/255]
            utils.plotRod(xyz1, xyz2, r_loc, color, ax)

        # deformed
        for i in range(self.N_elem):

            ind1 = elems[i, 0]
            ind2 = elems[i, 1]

            xyz1 = nodes_d[ind1, :]
            xyz2 = nodes_d[ind2, :]

            r_loc = np.sqrt(areas[i]/np.pi)

            rgb_loc = matplotlib.colors.colorConverter.to_rgb(coolwarm_cmap(norm(-rel_stress[i])))
            color = rgb_loc
            utils.plotRod(xyz1, xyz2, r_loc, color, ax)

        # BC
        for ind in cons:

            x_loc, y_loc, z_loc = nodes_d[ind, :]

            ax.plot([x_loc, x_loc], [y_loc, y_loc], [z_loc, z_loc], 'ko', markersize = 6)

        ax.plot([self.x_lb, self.x_lb, self.x_lb, self.x_lb, self.x_ub, self.x_ub, self.x_ub, self.x_ub],
        [self.y_lb, self.y_lb, self.y_ub, self.y_ub, self.y_ub, self.y_ub, self.y_lb, self.y_lb],
        [self.z_lb, self.z_ub, self.z_ub, self.z_lb, self.z_lb, self.z_ub, self.z_ub, self.z_lb], alpha = 0)

        ax.view_init(elev=elev, azim=azim)

        fig.savefig(job_name + "_" + str(elev) + "_" + str(azim) + "_stress_" + ".pdf", bbox_inches='tight') 

    def plot_buckling(self, elev, azim, job_name):

        coolwarm_cmap = matplotlib.cm.get_cmap('coolwarm')
        norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
        
        fig = plt.figure(figsize=(16, 16))
        ax = fig.add_subplot(111, projection='3d', proj_type = 'ortho')

        ax.pbaspect = [1,1,1]

        ax.set_xlim([0, 30])
        ax.set_ylim([0, 30])
        ax.set_zlim([-15, 15])
        ax.grid(False)
        ax._axis3don = False

        # jig shape
        for i in range(self.N_elem):

            ind1 = elems[i, 0]
            ind2 = elems[i, 1]

            xyz1 = nodes[ind1, :]
            xyz2 = nodes[ind2, :]

            r_loc = np.sqrt(areas[i]/np.pi)

            color = [192.0/255, 192.0/255, 192.0/255]
            utils.plotRod(xyz1, xyz2, r_loc, color, ax)

        # deformed
        for i in range(self.N_elem):

            ind1 = elems[i, 0]
            ind2 = elems[i, 1]

            xyz1 = nodes_d[ind1, :]
            xyz2 = nodes_d[ind2, :]

            r_loc = np.sqrt(areas[i]/np.pi)

            rgb_loc = matplotlib.colors.colorConverter.to_rgb(coolwarm_cmap(norm(bucklings[i])))
            color = rgb_loc
            utils.plotRod(xyz1, xyz2, r_loc, color, ax)

        # BC
        for ind in cons:

            x_loc, y_loc, z_loc = nodes_d[ind, :]

            ax.plot([x_loc, x_loc], [y_loc, y_loc], [z_loc, z_loc], 'ko', markersize = 6)

        ax.plot([self.x_lb, self.x_lb, self.x_lb, self.x_lb, self.x_ub, self.x_ub, self.x_ub, self.x_ub],
        [self.y_lb, self.y_lb, self.y_ub, self.y_ub, self.y_ub, self.y_ub, self.y_lb, self.y_lb],
        [self.z_lb, self.z_ub, self.z_ub, self.z_lb, self.z_lb, self.z_ub, self.z_ub, self.z_lb], alpha = 0)

        ax.view_init(elev=elev, azim=azim)
        fig.savefig(job_name + "_" + str(elev) + "_" + str(azim) + "_buckling_" + ".pdf", bbox_inches='tight') 



node_file_name = "./Wings/W_315/data_nodes.dat"
elem_file_name = "./Wings/W_315/data_elems.dat"
con_file_name = "./Wings/W_315/data_constraints.dat"
area_file_name = "./solution/old_opt_sol/areas.dat"
disp_file_name_list = ["./solution/new_case_1/sol_disp.txt",
    "./solution/new_case_2/sol_disp.txt",
    "./solution/new_case_3/sol_disp.txt",
    "./solution/old_case_1/sol_disp.txt",
    "./solution/old_case_2/sol_disp.txt",
    "./solution/old_case_3/sol_disp.txt"]
stress_file_name_list = ["./solution/new_case_1/sol_stress.txt",
    "./solution/new_case_2/sol_stress.txt",
    "./solution/new_case_3/sol_stress.txt",
    "./solution/old_case_1/sol_stress.txt",
    "./solution/old_case_2/sol_stress.txt",
    "./solution/old_case_3/sol_stress.txt"]
buckling_file_name_list = ["./solution/new_case_1/sol_buckling.txt",
    "./solution/new_case_2/sol_buckling.txt",
    "./solution/new_case_3/sol_buckling.txt",
    "./solution/old_case_1/sol_buckling.txt",
    "./solution/old_case_2/sol_buckling.txt",
    "./solution/old_case_3/sol_buckling.txt"]


case_name_list = ["new_case_1",
"new_case_2",
"new_case_3",
"old_case_1",
"old_case_2",
"old_case_3",]

for i in range(len(case_name_list)):
    nodes = np.loadtxt(node_file_name)
    elems = np.loadtxt(elem_file_name)
    elems = elems.astype(int)
    cons = np.loadtxt(con_file_name)
    cons = cons.astype(int)
    areas = np.loadtxt(area_file_name)
    disps = np.loadtxt(disp_file_name_list[i])
    rel_stress = np.loadtxt(stress_file_name_list[i])
    bucklings = np.loadtxt(buckling_file_name_list[i])

    nodes_d = nodes + disps


    elev_list = [90, 0, 30]
    azim_list = [0, 0, -90]

    trussPlotter_obj = trussPlotter(nodes, elems, cons, areas, disps, rel_stress, bucklings)
    for j in range(len(elev_list)):

        elev = elev_list[j]
        azim = azim_list[j]

        file_name = "./solution/" + case_name_list[i] + "/" + case_name_list[i]
        trussPlotter_obj.plot_stress(elev, azim, file_name)
        trussPlotter_obj.plot_buckling(elev, azim, file_name)




