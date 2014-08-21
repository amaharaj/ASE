from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write
from ase.io import read
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, \
                                        Stationary, ZeroRotation
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.constraints import Hookean
from ase.constraints import FixAtoms

from ase import Atoms

from ase.io.trajectory import PickleTrajectory
from ase import *

import numpy as np

text_file = open("Output.txt", "w")
atoms = read('HH.xyz')

atoms.set_calculator(Dftb(label = 'HH.xyz',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_H='"s"',	
                         ))


# Constrain all atoms with an x value lower than -2.0
i = 0
array = []
array2 = []
OP = []

while i <= 1:
	
	if atoms.positions[i][0] <= 0.2:
		array.append(i)
		i += 1
	else:
		array2.append(i)
		i += 1
			
c = FixAtoms(indices = array)
print "array1: ", array, "array2: ", array2

# Apply Hookean Force to all atoms with an x value greater than -2.0
for a in array2:
#	D = atoms.get_distance(a0 = 0, a1 =1) 
#	c2 = Hookean(a1=0, a2 = 1, rt = D, k = 15)
	c2 = Hookean(a1=1, a2 = (1,0,0,0),  k=1)
#	atoms.set_constraint([c, c2])
	atoms.set_constraint([c2])

def print_distance(a = atoms):
#	if len(OP) > 22:


#		OP_f = [OP[i:i+2] for i in xrange(0,len(OP),2)]
#		OP_f_text = np.savetxt('out.tmp', OP_f)
#		OP_f_infile = open('out.tmp','r')
#		OP_f_text = OP_f_infile.read()
#		text_file.write(OP_f_text)
#		text_file.close()
#		exit()

	x_dist = a.get_distance(a0 = 0, a1 = 1)
	OP.append(x_dist)	
	epot = a.get_potential_energy() / len(a)
	ekin = a.get_kinetic_energy() / len(a)
	ETOT = epot+ekin
	OP.append(ETOT)

qn = QuasiNewton(atoms, trajectory='qn.traj')
qn.attach(print_distance)
qn.run(fmax=0.001)
qn.attach(print_distance)
write('qn.final.xyz', atoms)

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300*units.kB)
print 'Removing linear momentum and angular momentum'
Stationary(atoms) # zero linear momentum
ZeroRotation(atoms) # zero angular momentum

# We want to run MD using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 0.1*units.fs, trajectory='moldyn4.traj') # save trajectory.

#Function to print the potential, kinetic and total energy.
def printenergy(t=atoms):    #store a reference to atoms in the definition.
    epot = t.get_potential_energy() / len(t)
    ekin = t.get_kinetic_energy() / len(t)

    print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
           (epot, ekin, ekin/(1.5*units.kB), epot+ekin))
    ETOT = ekin + epot
    print ETOT
    x_dist = atoms.get_distance(a0 = 0, a1 = 1)
    forces = atoms.get_forces()
    OP.append(x_dist)
    OP.append(ekin+epot)
    
dyn.attach(printenergy, interval=10)

# Now run the dynamics
dyn.run(400)
OP_f = [OP[i:i+2] for i in xrange(0,len(OP),2)]
OP_f_text = np.savetxt('out.tmp', OP_f)         # reads array and saves into a file as a string
OP_f_infile = open('out.tmp','r')       # open file to write string array onto
OP_f_text = OP_f_infile.read()

text_file.write(OP_f_text)
text_file.close()
