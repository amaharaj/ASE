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
from ase.constraints import Filter
from ase import Atoms

import math

atoms = read('znorg.xyz')

atoms.set_calculator(Dftb(label = 'znorg.xyz',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
			 Hamiltonian_MaxAngularMomentum_Zn='"d"',
			 Hamiltonian_MaxAngularMomentum_O = '"p"',
                         ))


# Constrain all atoms with an x value lower than -2.0

i = 0
array = []
array2 = []
while i < 6:
	
	if atoms.positions[i][0] <= -2.0:
		array.append(i)
		i += 1
	else:
		array2.append(i)
		i += 1	

c = FixAtoms(indices = array)

# Constrain all atoms with an x value greater than -2.0 to oscillate between some theta range

f1 = Filter(atoms, indices = array)
print f1.get_positions()



#Apply Hookean Force to all atoms with an x value greater than -2.0
for a in array2:
	D = atoms.positions[a][0] - 0.5
	c2 = Hookean(a1=a, a2=(1, 0, 0, D), k=1.0)
	m = math.pi
	print atoms.positions
#	f1 = Filter(atoms, indices = [a])
#	for atom in atoms:
#		if atoms.positions[a][0] > -2.0:
#	atoms = atoms[a]
#	if atoms.positions[a][0] > -2.0:
	f1.rotate('y',m) 
	#rint atoms.positions, "!!"
	atoms.set_constraint([c, c2])

# relax with Quasi Newtonian
qn = QuasiNewton(atoms, trajectory='qn.traj')
qn.run(fmax=0.001)
write('qn.final.xyz', atoms)

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300*units.kB)
print 'Removing linear momentum and angular momentum'
Stationary(atoms) # zero linear momentum
ZeroRotation(atoms) # zero angular momentum

# We want to run MD using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 0.1*units.fs, trajectory='moldyn4.traj') # save trajectory.

#Function to print the potential, kinetic and total energy.
def printenergy(a=atoms):    #store a reference to atoms in the definition.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
           (epot, ekin, ekin/(1.5*units.kB), epot+ekin))
dyn.attach(printenergy, interval=10)

# Now run the dynamics
printenergy()
dyn.run(2000)
