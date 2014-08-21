from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write, read
from ase.constraints import FixAtoms, Hookean

atoms = read('HH.xyz')
atoms.set_calculator(Dftb(label='hh',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))

c = Hookean(a1=0, a2=1, rt=0.1, k=50.)
atoms.set_constraint(c)
dyn = QuasiNewton(atoms, trajectory='atoms.traj')
dyn.run(steps=10)
print atoms.get_distance(a0=0,a1=1)
write('test.final.xyz', atoms)
