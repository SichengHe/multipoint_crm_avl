# Purpose

In the repo, we generate the multi-scenario load for a wing structure (uCRM-9).
The truss model of the wing was developed in a separate repo and we saved here in this repo.
The aerodynamic load is computed using AVL. 
The load is first computed over an underlying aerodynamic grids within AVL. 
Then we transfer the load to the structure.

Following is a detailed explanations of functions of different files.

The overall file structure is

* `aerodynamic` stores files for aerodynamic load computation for all scenarios.
The inputs are flight conditions in `driver.py` file for the flight conditions.
The `wing_lehigh.avl` and `ucrm9_geo.py` store the geometry of the aerodynamic mesh.
The outputs are the load for the truss structure.
* `Wings` stores the structure of the wing with different fidelity (grid sizes).
* `transfer` transfers the aerodynamic load to the structure.
A very simple transfer method is applied here: 
we associate each aerodynamic node to the closest truss node.
After this, we sum the aerodynamic load contribution for each truss node.
In the process, we drop the moment contribution and only keep the force.
We also split the load equally for the upper and lower surface of the truss.

The calling relationship between the files are given as follows


|-- `driver.py`\
| &nbsp; &nbsp; &nbsp; |-- `getLoad.py`\
| &nbsp; &nbsp; &nbsp; | &nbsp; &nbsp; &nbsp; |-- `wing_lehigh.avl`\
| &nbsp; &nbsp; &nbsp; |-- `transfer.py`

* `driver.py` is the main file which is responsible to generate loads for all the cases.
* `getLoad.py` generate one load case using the `wing_lehigh.avl` template with given parameters.
* `wing_lehigh.avl` is the template files used for multi-scenario load computation.
* `transfer.py` transfer the load from the aerodynamic mesh to the truss.

Besides the aforementioned main structure to generate load, we have the following files

* `ucrm9_geo.py` is used to compute the area and span etc of the wing.
It also stores the geometry of the wing.
* `calcaulator_twohalfg.py` is used to compute $C_L$ for the $2.5~g$ manuever case.


# TODO

* The structure of the files can be improved.
The `driver.py` is better to be placed outside the three components and holds a position as the main function of the whole directory rather than the position it holds right now.
* The `transfer.py` class might be overly simple.