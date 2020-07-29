import numpy as np

filename_sol = "solution_315_scenario_0.dat"

# start and end line of the bars
ind_beg_elem = 2
ind_end_elem = 316

length_scale = 100.0 # in cm convert to m

areas = np.zeros(ind_end_elem - ind_beg_elem + 1)

i = 0
j = 0
with open(filename_sol) as f:
    for line in f:
        if (ind_beg_elem<=j and j<=ind_end_elem):
            line_split = line.split()

            areas_loc = float(line_split[1])/length_scale**2 # m^2

            areas[i] = areas_loc

            i += 1
            
        j += 1

np.savetxt("areas.dat", areas)