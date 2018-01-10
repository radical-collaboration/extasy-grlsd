import mdtraj
import os

os.system('cat out*.gro >> tmp.gro')

file=mdtraj.load('tmp.gro')
selection1 = file.topology.select("symbol != 'H'")
file2=file.atom_slice(selection1, inplace=True)
file2.save_gro('tmpha.gro')

