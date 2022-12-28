import os, sys, platform

def spheres_2_pymol_script(spheres, resultFileName, increaseRadius):
    with open(resultFileName, 'w') as f:
        f.writelines('''
from pymol.cgo import *
from pymol import cmd

spheres = [
    COLOR, 1.0, 1.0, 1.0,
''')

        for x,y,z,r in spheres:
            f.writelines('''    SPHERE, {}, {}, {}, {},\n'''.format(x,y,z,r))

        f.writelines('''    ]

cmd.load_cgo(spheres, 'spheres_{}', 1)'''.format(increaseRadius))

libPath = os.getcwd()
pdbFilepath = os.getcwd()
resultFilepath = os.getcwd()

if platform.system() == 'Windows':
    libPath = libPath + str("\\lib\\")
    pdbFilepath += "\\data\\"
    resultFilepath += "\\results\\"

else:
    libPath = libPath + str("/lib/")
    pdbFilepath += "/data/"
    resultFilepath += "/results/"

sys.path.insert(0, libPath)
import PyMGOS

pdbFileName = pdbFilepath+'6g2j.pdb'
includeHETATM = False

if os.path.exists(pdbFileName.replace('.pdb', '.a.qtf')):
    os.remove(pdbFileName.replace('.pdb', '.a.qtf'))

MG = PyMGOS.MolecularGeometry()

if includeHETATM:
    MG.load(pdbFileName)
else:
    MG.load_except_PDB_HETATM(pdbFileName)

MG.preprocess()

increaseRadius = 0.0
increaseStep = 0.5

numAllAtoms = MG.get_all_atoms().size()
LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(increaseRadius)
print(increaseRadius, LRboudnaryAtomSet.size())

while LRboudnaryAtomSet.size() > (numAllAtoms/5):
    increaseRadius += increaseStep

    LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(increaseRadius)
    print(increaseRadius, LRboudnaryAtomSet.size())

LRboudnaryAtoms = []
for atom in LRboudnaryAtomSet:
    sphere = atom.geometry()

    x = sphere.center().x()
    y = sphere.center().y()
    z = sphere.center().z()
    r = sphere.radius()

    LRboudnaryAtoms.append((x, y, z, r+increaseRadius))

resultFileName = pdbFileName.replace('.pdb', '_spheres_{}.py'.format(increaseRadius))
spheres_2_pymol_script(LRboudnaryAtoms, resultFileName, increaseRadius)

    

    # LRVoids = MG.compute_voids_of_Lee_Richards_model(increaseRadius)
    # print(LRVoids.size())

    
pass
