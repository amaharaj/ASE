import math
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

atoms = read('h2o.xyz')

atoms.set_calculator(Dftb(label = 'h2o.xyz',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
			 Hamiltonian_MaxAngularMomentum_O = '"p"',
                         ))

a = math.pi/2
atoms.rotate('z',a)

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
