from ase.calculators.dftb import Dftb
from ase.optimize import QuasiNewton
from ase.io import write, read
from ase.constraints import FixAtoms, Hookean_Always

import numpy as np

text_file = open("Output.txt", "w")

atoms = read('HH.xyz')
atoms.set_calculator(Dftb(label='hh',
                         atoms=atoms,
                         Hamiltonian_MaxAngularMomentum_='',
                         Hamiltonian_MaxAngularMomentum_H='"s"',
                         ))

c = Hookean_Always(a1=0, a2=1, k=0, rt=0.6)
atoms.set_constraint(c)
OP = []

def print_distance(a=atoms):
        distance = a.get_distance(a0=0,a1=1)
	OP.append(distance)
        epot = a.get_potential_energy() / len(a)
        ekin = a.get_kinetic_energy() / len(a)
        ETOT = epot+ekin
        OP.append(ETOT)

dyn = QuasiNewton(atoms, trajectory='atoms.traj')
dyn.attach(print_distance)
dyn.run(steps=10)
print atoms.get_distance(a0=0, a1=1)
write('test.final.xyz', atoms)

OP_f = [OP[i:i+2] for i in xrange(0,len(OP),2)]
OP_f_text = np.savetxt('out.tmp', OP_f)         # reads array and saves into a file as a string
OP_f_infile = open('out.tmp','r')       # open file to write string array onto
OP_f_text = OP_f_infile.read()

text_file.write(OP_f_text)
text_file.close()
                     
