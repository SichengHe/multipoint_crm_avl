import sys
sys.path.append("./aerodynamic")

from getLoad import getLoad
import numpy as np


# ================================
# Preprocessing
# ================================

# --------------------------------
# Compute geometry and CL_factor
# --------------------------------

raw_data = np.array(
[[22.98279158699809, 0.0],
[36.59655831739962, 0.0],
[25.08728179551122, 2.8336658354114803],
[36.957605985037404, 2.8336658354114803],
[31.07231920199501, 10.866832917705743],
[38.35411471321696, 10.866832917705743],
[45.137157107231914, 29.390523690773072],
[47.83042394014963, 29.390523690773072]])


raw_data = raw_data - raw_data[0, :]



# chord
chord0 = raw_data[1, 0] - raw_data[0, 0]
chord1 = raw_data[3, 0] - raw_data[2, 0]
chord2 = raw_data[5, 0] - raw_data[4, 0]
chord3 = raw_data[7, 0] - raw_data[6, 0]

dh_01 = raw_data[2, 1] - raw_data[0, 1]
dh_12 = raw_data[4, 1] - raw_data[2, 1]
dh_23 = raw_data[6, 1] - raw_data[4, 1]

print("raw_data", raw_data)
print("chord0", chord0)
print("chord1", chord1)
print("chord2", chord2)
print("chord3", chord3)
print("dh_01", dh_01)
print("dh_12", dh_12)
print("dh_23", dh_23)

# area
area0 = (chord0 + chord1)/2 * dh_01
area1 = (chord1 + chord2)/2 * dh_12
area2 = (chord2 + chord3)/2 * dh_23

print("area0", area0)
print("area1", area1)
print("area2", area2)

area_aero = (area1 + area2) * 2
area_ref = (area0 + area1 + area2) * 2
CL_factor = area_ref / area_aero

print("area_aero", area_aero)
print("area_ref", area_ref)
print("CL_factor", CL_factor)


# span
span = raw_data[-1, 1]*2
print("span", span)



# --------------------------------
# Compute 2.5 g CL
# --------------------------------

# parameters
MTOW = 297500 # kg
Area = 410.8 # m^2 reference area computed in ucrm9_geo.py
Ma = 0.64 
CL_factor = 1.21 # reference area / true area in ucrm9_geo.py

# height 0
T = 288.150
rho = 1.22500

gamma = 1.4
R = 287.058

# calculation
sos = np.sqrt(gamma*R*T)
velocity = sos*Ma
weight = MTOW*9.8
loadFactor = 2.5
CL_manuever = loadFactor * weight/(0.5*rho*velocity**2*Area)

print('CL_manuever', CL_manuever)




# ================================
#       Aerodynamic analysis
# ================================

# conduct aerodynamic analysis with avl 
# under multiple operation conditions
# NOTICE: THE WING AREAS NOT MATCH!!! (DID NOT INCLUDE THE ROOT PORTION)!!!

# ref
# https://arc.aiaa.org/doi/pdf/10.2514/6.2014-3274
# https://arc.aiaa.org/doi/pdfplus/10.2514/1.J056603 (Table 6)
# https://www.digitaldutch.com/atmoscalc/
# 1g, at 37000 ft        

# --------------------------------
# cruise
# --------------------------------

T_oper = 216.650
rho_oper = 0.348331

Cl_Ma_pair = [
    [0.5, 0.85],
    [0.475, 0.85],
    [0.525, 0.85],
    [0.512, 0.84],
    [0.488, 0.86]
]

sum_weight_list = []
for i in range(len(Cl_Ma_pair)):
        
    Ma_loc = Cl_Ma_pair[i][1]
    Cl_loc = Cl_Ma_pair[i][0]

    node, load = getLoad(Cl_loc, Ma_loc, rho_oper, T_oper, CL_factor = CL_factor)

    sum_weight_list.append(np.sum(load)/9.8*2)

# --------------------------------
# 2.5 g manuever
# --------------------------------

Ma = 0.64 
rho = 1.22500
T = 288.150
node, load = getLoad(CL_manuever, Ma, rho, T, CL_factor = CL_factor)

# Maximum take-off weight (MTOW) 297 500 kg
# Maximum landing weight (MLW) 213 180 kg
# Maximum zero fuel weight (MZFW) 195 040 kg
# Operational empty weight 138 100 kg

# sanity check
print(sum_weight_list)
print(np.sum(load)/9.8*2)

# THE RESULT IS SMALL!
#[188421.66188482844, 179142.77207199342, 197679.83922435757, 188356.79997978732, 188324.49259952435]


# ================================
#        transfer the load
# ================================
# get the directory to be interpolated with
import os
rootdir = 'Wings'

dir_list = []
dir_list_total = os.listdir(rootdir)

for i in range(len(dir_list_total)):
    dir_loc = dir_list_total[i]
    if dir_loc.startswith('W_'):
        dir_list.append(rootdir + '/' + dir_loc)

# transfer the load
import sys
sys.path.insert(0, 'transfer')
from transfer import getLoad_new

os.system('mkdir %s' % "./OUT_STRUCT")
os.system('mkdir %s' % "./OUT_STRUCT/Wings")
for i in range(len(Cl_Ma_pair)):
    
    Cl_loc = Cl_Ma_pair[i][0]
    Ma_loc = Cl_Ma_pair[i][1]
    dir = "OUT_AERO" + "/" + "Cl%.3fMa%.2f" %(Cl_loc, Ma_loc)

    node = np.loadtxt(dir+'/'+'node.txt')
    load = np.loadtxt(dir+'/'+'load.txt')

    # print("np.sum(load)", np.sum(load))

    for j in range(len(dir_list)):

        dir_loc = dir_list[j]

        os.system('mkdir %s' % "./OUT_STRUCT" + "/" + dir_loc)

        node_new = np.loadtxt(dir_loc + '/' + 'data_nodes.dat')

        load_new = getLoad_new(node, load, node_new)

        np.savetxt("./OUT_Struct" + "/" + dir_loc + '/' + "data_forces_Cl%.3fMa%.2f.dat" %(Cl_loc, Ma_loc), load_new)

        print("np.sum(load_new)", np.sum(load_new))

        
# 2.5g
dir = "OUT_AERO" + "/" + "Cl%.3fMa%.2f" %(CL_manuever, Ma)
node = np.loadtxt(dir+'/'+'node.txt')
load = np.loadtxt(dir+'/'+'load.txt')

for j in range(len(dir_list)):

    dir_loc = dir_list[j]

    node_new = np.loadtxt(dir_loc + '/' + 'data_nodes.dat')

    load_new = getLoad_new(node, load, node_new)

    np.savetxt("./OUT_Struct" + "/" + dir_loc + '/' + "data_forces_Cl%.3fMa%.2f.dat" %(CL_manuever, Ma), load_new)

    print("np.sum(load_new)", np.sum(load_new))