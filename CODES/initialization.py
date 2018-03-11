import numpy
import scipy
import Bio
from pyrosetta import *
init()
from pyrosetta.toolbox import cleanATOM
cleanATOM("1yy8.pdb")
#from pyrosetta.toolbox import pose_from_rcsb
pose = pose_from_pdb("1YY8.clean.pdb")
print(pose)
print(pose.sequence())
print("Protein has", pose.total_residue(), "residues.")
print(pose.residue(500).name())
print(pose.pdb_info().chain(500))
print(pose.pdb_info().number(500))
print(pose.phi(5))
print(pose.psi(5))
print(pose.chi(1, 5))
R5N = AtomID(1, 5)
R5CA = AtomID(2, 5)
R5C = AtomID(3, 5)
print(pose.conformation().bond_length(R5N, R5CA))
print(pose.conformation().bond_length(R5CA, R5C))
from pyrosetta import PyMOLMover
pymol = PyMOLMover()
pymol.apply(pose)
# Calculating energy score
print("Calculating energy score of 6Q21 protein")
ras = pose_from_pdb("6q21.pdb")
print(ras)
print(ras.sequence())
from pyrosetta.teaching import *
scorefxn = get_fa_scorefxn()
print("Score function of 6Q21 protein is :")
print(scorefxn)
scorefxn2 = ScoreFunction()
scorefxn2.set_weight(fa_atr, 1.0)
scorefxn2.set_weight(fa_rep, 1.0)
print(scorefxn(ras))
scorefxn.show(ras)

print("==========================")
pose2 = pose_from_sequence('A'*10)
pmm = PyMOLMover()
pmm.apply(pose)

print(ras.energies().show(314))
r1 = ras.residue(24)
r2 = ras.residue(20)
a1 = r1.atom("N")
a2 = r2.atom("O")
print(etable_atom_pair_energies(a1, a2, scorefxn))
