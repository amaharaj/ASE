from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write
from ase.io import read

from ase.structure import molecule
test = read('initial.xyz')
test.set_calculator(Dftb(label='h2o',
                         atoms=test,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_O='"p"',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))

dyn = QuasiNewton(test, trajectory='test.traj')
dyn.run(fmax=0.001)
write('test.final.xyz', test)
