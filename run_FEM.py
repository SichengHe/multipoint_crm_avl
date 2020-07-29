import numpy as np
import copy

class FEM():

    '''
        A truss-based FEM solver.
        Takes in the nodes coordinates (nodes), topology (elems),
        constrained nodal indices (cons), load (forces), areas (areas).
        It can solve for the displacement. Afterwards, we can conduct
        postprocess to generate dimensionless yield stress 
        and buckling.
        
        Sizes:
        N: total node number (with constrained elements)
        M: total bar number
        N_con: constrained node number

        Dimensions:
        nodes: (N x 3), float
        elems: (M x 2), int
        cons: (N_con), int
        forces: (N x 3), float
        areas: (M), float

        Assumptions:
        1. nodes, elems, forces contain constrained nodes.
        2. R contain constrained nodes.
        3. elems contained nodes indexed starting from 0.
        4. cons order in monotonically increasing way.
    '''

    def __init__(self, nodes, elems, cons, forces, areas, E):

        self.nodes = nodes
        self.elems = elems
        self.cons = cons
        self.forces = forces
        self.areas = areas

        self.forces = np.delete(self.forces, self.cons, 0)
        self.nodes_wo_cons = np.delete(self.nodes, self.cons, 0)

        self.N = self.nodes.shape[0]
        self.M = self.elems.shape[0]
        self.N_con = self.cons.shape[0]

        self.E = E

        self.compute_R()

    def compute_R(self):

        e1 = np.array([1.0, 0, 0])
        e2 = np.array([0, 1.0, 0])
        e3 = np.array([0, 0, 1.0])

        # compute R
        self.R = np.zeros((self.M, 3 * self.N))
        for i in range(self.M):
            ind1 = self.elems[i, 0]
            ind2 = self.elems[i, 1]

            node1 = self.nodes[ind1, :]
            node2 = self.nodes[ind2, :]

            v12 = node2 - node1
            v12 /= np.linalg.norm(v12)
            v21 = - v12

            R11 = np.dot(v21, e1)
            R12 = np.dot(v21, e2)
            R13 = np.dot(v21, e3)
            R21 = np.dot(v12, e1)
            R22 = np.dot(v12, e2)
            R23 = np.dot(v12, e3)

            R11_ind = ind1 * 3
            R12_ind = ind1 * 3 + 1
            R13_ind = ind1 * 3 + 2
            R21_ind = ind2 * 3
            R22_ind = ind2 * 3 + 1
            R23_ind = ind2 * 3 + 2

            self.R[i, R11_ind] = R11
            self.R[i, R12_ind] = R12
            self.R[i, R13_ind] = R13

            self.R[i, R21_ind] = R21
            self.R[i, R22_ind] = R22
            self.R[i, R23_ind] = R23


        # expand cons to d.o.f index
        cons_expanded = np.zeros(self.N_con * 3, int)
        for i in range(self.N_con):
            cons_expanded[3 * i] = 3 * self.cons[i]
            cons_expanded[3 * i + 1] = 3 * self.cons[i] + 1
            cons_expanded[3 * i + 2] = 3 * self.cons[i] + 2

        self.R = np.delete(self.R, cons_expanded, 1)

    def compute_length(self):

        self.L = np.zeros(self.M)

        for i in range(self.M):

            ind1 = self.elems[i, 0]
            ind2 = self.elems[i, 1]

            node1 = self.nodes[ind1, :]
            node2 = self.nodes[ind2, :]

            self.L[i] = np.linalg.norm(node2 - node1)

    def __call__(self):

        self.compute_length()
        self.set_K()
        self.compute_u()

    def set_K(self):

        EI_l = np.zeros(self.M)
        for i in range(self.M):

            A = self.areas[i]
            EI_l[i] = self.E * A / self.L[i]

        # K w/o con
        self.K = self.R.T.dot(np.diag(EI_l)).dot(self.R)

    def compute_u(self):

        self.u_wo_cons = np.linalg.solve(self.K, self.forces.flatten()).reshape((-1, 3))

        self.u = copy.deepcopy(self.u_wo_cons)

        for con_loc in self.cons:
            self.u = np.insert(self.u, con_loc, 0.0, axis = 0)

    def get_u(self):

        return self.u

    def set_sigmaY(self, sigmaY):

        self.sigmaY = sigmaY

    def compute_stress(self):

        # computes sigma / sigmaY for each bar

        self.dimless_sigma = np.zeros(self.M)

        elongation = self.R.dot(self.u_wo_cons.flatten())
        for i in range(self.M):

            epsilon = elongation[i] / self.L[i]
            stress = self.E * epsilon

            self.dimless_sigma[i] = stress / self.sigmaY

    def get_stress(self):

        return self.dimless_sigma

    def compute_buckling(self):

        # computes max( - sigma / (gamma * area), 0)

        self.dimless_buckling = np.zeros(self.M)

        elongation = self.R.dot(self.u_wo_cons.flatten())
        for i in range(self.M):

            epsilon = elongation[i] / self.L[i]
            stress = self.E * epsilon

            gamma = np.pi * self.E / (4.0 * self.L[i]**2)

            self.dimless_buckling[i] = max(- stress / (gamma * self.areas[i]), 0)

    def get_buckling(self):

        return self.dimless_buckling

    def write_to_file(self, jobname):

        np.savetxt(jobname + "_" + "disp.txt",  self.u)
        np.savetxt(jobname + "_" + "stress.txt",  self.dimless_sigma)
        np.savetxt(jobname + "_" + "buckling.txt", self.dimless_buckling)

        
def FEM_wrapper(job_name, node_file_name, elem_file_name, con_file_name, force_file_name, area_file_name):

    node = np.loadtxt(node_file_name)
    elem = np.loadtxt(elem_file_name)
    elem = elem.astype(int)
    con = np.loadtxt(con_file_name)
    con = con.astype(int)
    force = np.loadtxt(force_file_name)
    area = np.loadtxt(area_file_name)

    force_wo_con = np.delete(force, con, 0)
    node_wo_con = np.delete(node, con, 0)

    # construct K
    # E = 68947.57293 * 1e6
    E = 69.0 * 1e9
    sigmaY = 270.0 * 1e6

    # construct and solve
    FEM_obj = FEM(node, elem, con, force, area, E)
    FEM_obj()

    # compute stress
    FEM_obj.set_sigmaY(sigmaY)
    FEM_obj.compute_stress()

    # compute buckling
    FEM_obj.compute_buckling()

    # write to file
    FEM_obj.write_to_file(job_name)

# load data
node_file_name = "./Wings/W_315/data_nodes.dat"
elem_file_name = "./Wings/W_315/data_elems.dat"
con_file_name = "./Wings/W_315/data_constraints.dat"
force_file_name_list = [ "./solution/new_case_1/data_forces_Cl0.611Ma0.64.dat",
    "./solution/new_case_2/data_forces_Cl-0.244Ma0.64.dat",
    "./solution/new_case_3/data_forces_Cl0.403Ma0.86.dat",
    "./solution/old_case_1/data_forces_0.dat",
    "./solution/old_case_2/data_forces_Cl0.262Ma0.64.dat",
    "./solution/old_case_3/data_forces_Cl0.500Ma0.85.dat"]
area_file_name = "./solution/old_opt_sol/areas.dat"

job_name_list = ["./solution/new_case_1/sol",
    "./solution/new_case_2/sol",
    "./solution/new_case_3/sol",
    "./solution/old_case_1/sol",
    "./solution/old_case_2/sol",
    "./solution/old_case_3/sol"]

for i in range(len(force_file_name_list)):

    force_file_name = force_file_name_list[i]
    job_name = job_name_list[i]

    FEM_wrapper(job_name, node_file_name, elem_file_name, con_file_name, force_file_name, area_file_name)