import numpy as np

force_file_name_list = [ "./solution/new_case_1/data_forces_Cl0.611Ma0.64.dat",
    "./solution/new_case_2/data_forces_Cl-0.244Ma0.64.dat",
    "./solution/new_case_3/data_forces_Cl0.403Ma0.86.dat",
    "./solution/old_case_1/data_forces_0.dat",
    "./solution/old_case_2/data_forces_Cl0.262Ma0.64.dat",
    "./solution/old_case_3/data_forces_Cl0.500Ma0.85.dat"]

for file_name in force_file_name_list:

    force = np.loadtxt(file_name)

    print(np.sum(force[:, -1]) / 9.8 / 1000)