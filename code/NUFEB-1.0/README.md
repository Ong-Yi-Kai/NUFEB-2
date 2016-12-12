NUFEB is an open source tool for Individual Based model (IBm) simulation. 
The tool is implemented as a user package within LAMMPS - 
a molecular dynamics simulator offering basic functionalities for 
Discrete Element Method (DEM) simulations. NUFEB aims to improve those capability
with the goal to apply it to biological modelling.

NUFEB is a freely-available open-source code, distributed under the terms
of the GNU Public License.

NUFEB development has been funded by the EPSRC project An New
Frontier in Design: The Simulation of Open Engineered Biological Systems
(NUFEB).

### Building

NUFEB requires GCC/G++ for a successful build. 
It has been rigorously tested on Ubuntu-14.10 and Fedora-22.

To compile this code, go to The directory:
$ cd NUFEB/src/

then execute the following commands to compile code in the /STUBS directory:
$ cd STUBS/
$ make clean
$ make
$ cd ..

Now, install the NUFEB and granular packages with the following instruction:
$ make yes-USER-NUFEB
$ make yes-GRANULAR

Finally, execute the following
command to compile the NUFEB executable:
$ make serial

### Running

You can run NUFEB by going to /examples directory and run:
$  ../src/lmp_serial < Inputscript.lammps

