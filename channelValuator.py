from protein import Protein

import os
import numpy as np

class ChannelValuator:
    def __init__(self):
        self._pdbFileName       = None
        self._protein           = Protein()
        self._includeHETATM     = False
        
        self._channelFilePaths  = None
        self._channels          = dict()
        
        self._min = np.array([np.inf, np.inf, np.inf]   , dtype=np.float32)
        self._max = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float32)

    def set_protein(self, pdbFileName, includeHETATM):
        from VDWradius import VDWRadius
        import os
        
        if not os.path.exists(pdbFileName):
            raise FileExistsError('No Protein File : {}'.format(pdbFileName))
        
        self._pdbFileName = pdbFileName
        self._includeHETATM = includeHETATM
        with open(self._pdbFileName, 'r') as f:
            # print("Protein file : {} Loaded.\n".format(proteinFilePath))

            proteinName = pdbFileName.split('/')[-1].replace('.pdb','')
            self._protein.set_name(proteinName)
            
            lines = f.readlines()
            if includeHETATM:
                for line in lines:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        line = line.rstrip()

                        atomName= line[-3:]
                        atomNum = int(line[6:12])
                        resName = line[17:20]

                        x       = np.float32(line[30:38])
                        y       = np.float32(line[38:46])
                        z       = np.float32(line[46:52])
                        
                        try:
                            r = VDWRadius[atomName.split()[-1].rjust(2,' ')]
                        except:
                            # raise BaseException('NO van der Waals radius : {}'.format(atomName))
                            parsedAtomName = ''
                            for char in atomName:
                                if char.isalpha():  
                                    parsedAtomName += char
                            r = VDWRadius[parsedAtomName.split()[-1].rjust(2,' ')]
                        
                        atom    = np.array([x,y,z,r], dtype=np.float32)
                        self._protein.add_atoms((atom, atomName, resName, atomNum))

                        # self._min = np.minimum(atom[:3] - r, self._min)
                        # self._max = np.maximum(atom[:3] + r, self._max)
            else:
                for line in lines:
                    if line.startswith("ATOM"):
                        line = line.rstrip()

                        atomName= line[-3:]
                        atomNum = int(line[6:12])
                        resName = line[17:20]

                        x       = np.float32(line[30:38])
                        y       = np.float32(line[38:46])
                        z       = np.float32(line[46:52])
                        
                        try:
                            r = VDWRadius[atomName.split()[-1].rjust(2,' ')]
                        except:
                            # raise BaseException('NO van der Waals radius : {}'.format(atomName))
                            parsedAtomName = ''
                            for char in atomName:
                                if char.isalpha():
                                    parsedAtomName += char
                            r = VDWRadius[parsedAtomName.split()[-1].rjust(2,' ')]
                        
                        atom    = np.array([x,y,z,r], dtype=np.float32)
                        self._protein.add_atoms((atom, atomName, resName, atomNum))

                        # self._min = np.minimum(atom[:3] - r, self._min)
                        # self._max = np.maximum(atom[:3] + r, self._max)

    def inflate_protein(self, cutoffRatio):
        import os, sys, platform

        pdbFileName = self._;pdbFileName


        libPath = os.getcwd() + "\\lib\\"
        if platform.system() != 'Windows':
            libPath = libPath.replace('\\', '/')

        
        import PyMGOS

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
        while LRboudnaryAtomSet.size() > (numAllAtoms/cutoffRatio):
            increaseRadius += increaseStep

            LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(increaseRadius)

        print(numAllAtoms, LRboudnaryAtomSet.size(), increaseRadius)

        LRboudnaryAtoms = []
        for atom in LRboudnaryAtomSet:
            sphere = atom.geometry()

            x = sphere.center().x()
            y = sphere.center().y()
            z = sphere.center().z()
            r = sphere.radius()

            LRboudnaryAtoms.append((x, y, z, r+increaseRadius))            


    def set_grid(self, gridSize):
        from grid import Grid
                
        self._gridSize = gridSize
        self._grid = Grid(self._min, self._max, gridSize)
        self._grid._split()

    def set_channels(self, channelFilePaths):
        from channel import Channel
        import os
        
        self._channelFilePaths = channelFilePaths
        for software in channelFilePaths.keys():
            channelFilePath = channelFilePaths[software]
            if os.path.exists(channelFilePath):
                # print("     Channel {} Loaded.".format(software))
                self._channels[software] = {
                    'channel'           : Channel(software),
                    'TotalVertices'     : 0,
                    'BuriedVertices'    : 0,
                    'RevealedVertices'  : 0
                }

                channelFunc = getattr(self, '_parse_channel_'+software)
                channelFunc(channelFilePath, self._channels[software]['channel'])
                
    def _parse_channel_MGOS(self, MGOSchannelFilePath, channels):
        with open(MGOSchannelFilePath) as f:
            lines = f.readlines()
            channel_check = False
            for line in lines:
                if line.startswith('channel'):
                    channel_check = True
                    continue
                if line.startswith('contribAtom'):
                    channel_check = False

                if channel_check is True and line.startswith('    SPHERE'):
                    _, x, y, z, r, _ = line.split(',')
                    x = np.float32(x)
                    y = np.float32(y)
                    z = np.float32(z)
                    r = np.float32(r)
                    sphere = (x,y,z,r)
                    channels.add_sphere(sphere)
                                    
    def _parse_channel_MOLE(self, MOLEchannelFilePath, channels):
        with open(MOLEchannelFilePath) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('def Tunnel'):
                    continue
                if line.startswith('  addAtom'):
                    sphereInfo = line.split('(')[1].replace(')', '')
                    _, _, r, x, y, z = sphereInfo.split(',')
                    x = np.float32(x)
                    y = np.float32(y)
                    z = np.float32(z)
                    r = np.float32(r)
                    sphere = (x,y,z,r)
                    channels.add_sphere(sphere)
    
    def _parse_channel_CAVER(self, CAVERchannelFilePath, channels):
        tunnelFiles = list()
        for (root, dirs, files) in os.walk(CAVERchannelFilePath):
            for filename in  files:
                if filename.endswith('.pdb'):
                    tunnelFiles.append(CAVERchannelFilePath + filename)
        for tunnelFile in tunnelFiles:
            with open(tunnelFile) as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('ATOM'):
                        splitted_line = line.split()
                        x       = np.float32(splitted_line[-4])
                        y       = np.float32(splitted_line[-3])
                        z       = np.float32(splitted_line[-2])
                        r       = np.float32(splitted_line[-1])
                        sphere  = (x,y,z,r)
                        channels.add_sphere(sphere)

    def _parse_channel_3V(self, VVVchannelFilePath, channels):
        pass
    
    def find_ground_truth_grid(self, pdbFileName, includeHETATM, cutoffRatio):
        self._define_surface_grids(pdbFileName, includeHETATM, cutoffRatio)
        

    def _define_surface_grids(self, pdbFileName, includeHETATM, cutoffRatio):
        import os, sys, platform

        libPath = os.getcwd() + "\\lib\\"
        if platform.system() != 'Windows':
            libPath = libPath.replace('\\', '/')

        sys.path.insert(0, libPath)
        import PyMGOS

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
        while LRboudnaryAtomSet.size() > (numAllAtoms/cutoffRatio):
            increaseRadius += increaseStep

            LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(increaseRadius)

        print(numAllAtoms, LRboudnaryAtomSet.size(), increaseRadius)

        LRboudnaryAtoms = []
        for atom in LRboudnaryAtomSet:
            sphere = atom.geometry()

            x = sphere.center().x()
            y = sphere.center().y()
            z = sphere.center().z()
            r = sphere.radius()

            LRboudnaryAtoms.append((x, y, z, r+increaseRadius))

        # resultFileName = pdbFileName.replace('.pdb', '_spheres_{}.py'.format(increaseRadius))
        # self._spheres_2_pymol_script(LRboudnaryAtoms, resultFileName, increaseRadius)

        for sphere in LRboudnaryAtoms:
            includingGridVertices = self._grid.get_including_gridVertices_of((sphere,))
            
            for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                vertex = ivertex.item()
                vertex.check_intersection_with_increased_sphere(sphere)

    def _define_interior_vertices(self):
        self._grid


    def _spheres_2_pymol_script(spheres, resultFileName, increaseRadius):
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

    def check_intersection_with_atoms(self):
        for atom in self._protein.atoms:
            includingGridVertices = self._grid.get_including_gridVertices_of(atom)
            
            for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                vertex = ivertex.item()
                vertex.check_intersection_with_atom(atom)

    def check_intersection_with_channels(self):
        for software in self._channels.keys():
            for sphere in self._channels[software]['channel']._spheres:
                includingGridVertices = self._grid.get_including_gridVertices_of((sphere,))
                
                if includingGridVertices.size == 0:
                    continue
                
                for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                    vertex = ivertex.item()

                    where = vertex.has_intersection_with_sphere(sphere)
                    if where == 'Boundary':
                        self._channels[software]['channel'].add_intersecting_vertex(vertex)
                        self._channels[software]['channel'].add_boundary_vertex(vertex)
                    elif where == 'Inside':
                        self._channels[software]['channel'].add_intersecting_vertex(vertex)
                    elif where == 'Outside':
                        pass
                

    def validate_channel(self, includeWater):
        for software in self._channels.keys():
            self._channels[software]['TotalVertices'] = len(self._channels[software]['channel']._intersectingVertices)

            for vertex in self._channels[software]['channel']._intersectingVertices:
                # vertex = ivertex.item()
                if len(vertex.intersectingAtoms) > 0:
                    self._channels[software]['BuriedVertices'] += 1
                else:
                    self._channels[software]['RevealedVertices'] += 1

    def write_results(self, resultFilePath):
        with open(resultFilePath, 'w') as f:
            for software in self._channels.keys():
                f.writelines('''
    Software        :   {}
    TotalVolume     :   {}
    BuriedVertices  :   {}
    BoundaryVertices:   {}
    RevealedVertices:   {}

        '''.format(
                software,
                self._channels[software]['TotalVertices'],
                self._channels[software]['BuriedVertices'],
                len(self._channels[software]['channel']._boundaryVertices),
                self._channels[software]['RevealedVertices'])
                )
            
           
    def write_PyMOL_script(self, pyMOLScriptFile, includeWater):
        radius = 0.15
        with open(pyMOLScriptFile, 'w') as f:
            f.writelines(
'''
from pymol.cgo import *
from pymol import cmd
import os
view = cmd.get_view()     
'''
            )
            for software in self._channels.keys():
                RevealedVertices    = ['RevealedVertices_{} = [\n'.format(software), '\tCOLOR, 1.000000, 1.000000, 1.000000,\n']
                BuriedVertices      = ['BuriedVertices_{} = [\n'.format(software),   '\tCOLOR, 0.000000, 0.000000, 1.000000,\n']
                for vertex in self._channels[software]['channel']._intersectingVertices:
                    if len(vertex.intersectingAtoms) > 0:
                        BuriedVertices.append('\tSPHERE, {}, {}, {}, {},\n'.format(vertex.point[0], vertex.point[1], vertex.point[2], radius))
                    else:
                        RevealedVertices.append('\tSPHERE, {}, {}, {}, {},\n'.format(vertex.point[0], vertex.point[1], vertex.point[2], radius*1.5))
                BuriedVertices.append(']\n')
                RevealedVertices.append(']\n')
                
                f.write(''.join(BuriedVertices))
                f.write("cmd.load_cgo(BuriedVertices_{}, 'BuriedVertices_{}')\n".format(software, software))
                f.write(''.join(RevealedVertices))  
                f.write("cmd.load_cgo(RevealedVertices_{}, 'RevealedVertices_{}')\n".format(software, software))
                
                f.write("cmd.group('{}','BuriedVertices_{}', 'add')\n".format(software, software))
                f.write("cmd.group('{}','RevealedVertices_{}', 'add')\n".format(software, software))