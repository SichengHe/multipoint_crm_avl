import os
from subprocess import Popen, PIPE
import numpy as np


def getLoad(Cl, Ma, rho, T, CL_factor = 1.0, outdirectory = "./OUT_AERO", indirectory = "./aerodynamic"):
    
    # Input:   in CL, Ma
    # Process: mkdir for all combinations of Cl, Ma;
    #          avl files for all Ma;
    #          run avl files for all Cl and get FE files;
    # Output: FE data

    gamma = 1.4
    R = 287.058
    vel = np.sqrt(gamma*R*T) * Ma
    dyn_pressure = 0.5*rho*vel**2


    # mkdir for all combinations of Cl, Ma
    os.system('mkdir %s' % outdirectory)
    dir = outdirectory + '/' + 'Cl%.3fMa%.2f' %(Cl, Ma)
    os.system('mkdir %s' % dir)


    # avl files for all Ma;
    filename = dir +'/'+ "wing_lehigh_Ma%.2f.avl" % Ma
    lines = open(indirectory + "/" + "wing_lehigh.avl", "rt")
    output = open(filename, "wt")

    ind = 0
    for line in lines:
        if ind == 1:
            output.write("%.2f" %Ma + " " + "Ma" + "\n")
        else:
            output.write(line+"\n")
        ind += 1

    lines.close()
    output.close()

    # run avl files for all Cl and get FE files
    rawFEname = dir+"/" +"FE_raw.txt"
    p = Popen(indirectory + "/" + "avl", stdin=PIPE) 
    p.communicate(os.linesep.join(["LOAD", indirectory + "/" + "wing_lehigh.avl",\
    "OPER",\
    "A", "C", '%.4f' %(Cl * CL_factor),\
    "X",\
    "FE", rawFEname, "O"]))
    # p = Popen(["./avl", "LOAD", "wing_lehigh.avl"]) 

    # get FE data
    FE = []

    ind_aft_strip = 0
    flag_inside_strip = False

    ind_aft_I = 0
    flag_inside_I = False 
    with open(rawFEname, "r") as ins:  
        for line in ins:
            line_loc = line.split()
            print(line_loc)

            # get the width
            if len(line_loc)>0:
                if line_loc[0] == 'Strip':
                    ind_aft_strip = 0

                    flag_inside_strip = True

            if ind_aft_strip == 2:
                width = np.float(line_loc[6])
                print("width", width)

                ind_aft_strip = 0
                flag_inside_strip = False

            if flag_inside_strip:
                ind_aft_strip += 1

            # get the x, y, z, area, cp
            if len(line_loc)>0:
                if line_loc[0] == 'I':
                    
                    ind_aft_I = 0
                    flag_inside_I = True

            if len(line_loc)==0 or line_loc[0][0]=='-':
                
                ind_aft_I = 0
                flag_inside_I = False

            
            if flag_inside_I and ind_aft_I>=1:
                
                x_loc = np.float(line_loc[1])
                y_loc = np.float(line_loc[2])
                z_loc = np.float(line_loc[3])
                dx = np.float(line_loc[4])
                cp = np.float(line_loc[6])

                area = dx*width
                force = area*dyn_pressure*cp

                FE_loc = [x_loc, y_loc, z_loc, force]

                FE.append(FE_loc)
            
            if flag_inside_I:
                
                ind_aft_I += 1

    FE = np.asarray(FE)

    node = FE[:, 0:3]
    load = FE[:, 3]

    nodename = dir + "/" + "node.txt"
    np.savetxt(nodename, node)
    loadname = dir + "/" + "load.txt"
    np.savetxt(loadname, load)

    return node, load


