from symbol import atom
from protein import Protein

import os
import numpy as np

class ChannelValuator:
    def __init__(self):
        self._proteinFilePath   = None
        self._protein           = Protein()
        
        self._channelFilePaths  = None
        self._channels          = dict()
        
        self._emptyGrids    = list()
        self._fullGrids     = list()
        self._partialGrids  = list()
        
        self._min = np.array([np.inf, np.inf, np.inf]   , dtype=np.float32)
        self._max = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float32)
    

    def set_protein(self, proteinFilePath, includeHETATM):
        from VDWradius import VDWRadius
        import os
        
        if not os.path.exists(proteinFilePath):
            raise FileExistsError('No Protein File : {}'.format(proteinFilePath))
        
        self._proteinFilePath = proteinFilePath
        with open(self._proteinFilePath, 'r') as f:
            # print("Protein file : {} Loaded.\n".format(proteinFilePath))

            proteinName = proteinFilePath.split('/')[-1].replace('.pdb','')
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

                        self._min = np.minimum(atom[:3] - r, self._min)
                        self._max = np.maximum(atom[:3] + r, self._max)
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

                        self._min = np.minimum(atom[:3] - r, self._min)
                        self._max = np.maximum(atom[:3] + r, self._max)


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

    def _parse_channel_3V(self, VVVchannelFilePath,channels):
        pass
    

    def check_intersection_with_atoms(self):
        for atom in self._protein.atoms:
            includingGridVertices = self._grid.get_including_gridVertices_of(atom)

            isWater = False
            if atom[2] == 'HOH':
                isWater = True
            
            for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                vertex = ivertex.item()
                vertex.check_intersection_with_atom(atom, isWater)

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

                if includeWater:
                    if len(vertex.intersectingAtoms) > 0 or len(vertex.intersectingWaters) > 0:
                        self._channels[software]['BuriedVertices'] += 1
                    else:
                        self._channels[software]['RevealedVertices'] += 1
                else:
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
                
        # with open(resultFilePath.replace('.txt', '_details.txt'), 'w') as f:
        #     for software in self._channels.keys():
        #         for vertex in self._channels[software]['channel']._intersectingVertices:
        #             # vertex = ivertex.item()
        #             for atom in vertex.intersectingAtoms:
        #                 xyzr, atomName, resName, atomNum = atom
        #                 dist = np.linalg.norm(vertex.point - xyzr[0:3])
        #                 f.write('{}|{}|{}|{}|{}|{}|{}|{}\n'.format(
        #                     atomName, resName, atomNum,
        #                     vertex.point, xyzr[0:3],
        #                     xyzr[3], dist, xyzr[3]-dist)
        #                 )

        #         f.writelines('BoundaryVertices\n')
        #         for vertex in self._channels[software]['channel']._boundaryVertices:
        #             # vertex = ivertex.item()
        #             for atom in vertex.intersectingAtoms:
        #                 xyzr, atomName, resName, atomNum = atom
        #                 dist = np.linalg.norm(vertex.point - xyzr[0:3])
        #                 f.write('{}|{}|{}|{}|{}|{}|{}|{}\n'.format(
        #                     atomName, resName, atomNum,
        #                     vertex.point, xyzr[0:3],
        #                     xyzr[3], dist, xyzr[3]-dist)
        #                 )
           
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
                    if includeWater:
                        if len(vertex.intersectingAtoms) > 0 or len(vertex.intersectingWaters) > 0:
                            BuriedVertices.append('\tSPHERE, {}, {}, {}, {},\n'.format(vertex.point[0], vertex.point[1], vertex.point[2], radius))
                        else:
                            RevealedVertices.append('\tSPHERE, {}, {}, {}, {},\n'.format(vertex.point[0], vertex.point[1], vertex.point[2], radius*1.5))
                    else:
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

#     def write_grids_into_PyMOL_script(self, gridFilePath):
#         with open(gridFilePath, 'w') as f:
#             lineWidth = 1.5            
#             f.writelines(
# '''
# from pymol.cgo import *
# from pymol import cmd
# import os
# view = cmd.get_view()     
# '''
#             )
#             self._write_grids(lineWidth, f)
#             self._write_channel_grids(lineWidth, f)
    
#     def _write_grids(self, lineWidth, f):
#         f.writelines(
# '''
# EmptyGrids = [
#     LINEWIDTH, {},
#     BEGIN, LINES,
#     COLOR, 1.000000, 1.000000, 1.000000,
# '''.format(float(lineWidth)))
#         for cell in self._emptyGrids:
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))            
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
                
#         f.writelines('''    END
# ]
# cmd.load_cgo(EmptyGrids, 'EmptyGrids')


# ''')
#         f.writelines(
# '''
# FullGrids = [
#     LINEWIDTH, {},
#     BEGIN, LINES,
#     COLOR, 1.000000, 0.000000, 0.000000,
# '''.format(lineWidth)
#         )
#         for cell in self._fullGrids:
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
                
#         f.writelines('''    END
#     ]
# cmd.load_cgo(FullGrids, 'FullGrids')

            
# ''')
        
#         f.writelines(
# '''
# PartialGrids = [
#     LINEWIDTH, {},
#     BEGIN, LINES,
#     COLOR, 0.000000, 0.000000, 1.000000,
# '''.format(lineWidth)
#         )
#         for cell in self._partialGrids:
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))                
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#             f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
                
#         f.writelines('''    END
#     ]
# cmd.load_cgo(PartialGrids, 'PartialGrids')     
# ''')
#         f.writelines('''
# cmd.group('Grids', 'EmptyGrids', 'add')
# cmd.group('Grids', 'FullGrids', 'add')
# cmd.group('Grids', 'PartialGrids', 'add')


# ''')       
        
#     def _write_channel_grids(self, lineWidth, f):
#         for software in self._channels.keys():
            
#             f.writelines(
# '''
# {}_Grids = [
#     LINEWIDTH, {},
#     BEGIN, LINES,
#     COLOR, 1.000000, 1.000000, 1.000000,
# '''.format(software, lineWidth))
#             for cell in self._channels[software]['totalGrids']:
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][0].x, cell._vertices[0][0][0].y, cell._vertices[0][0][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][0][1].x, cell._vertices[0][0][1].y, cell._vertices[0][0][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))                
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][1].x, cell._vertices[0][1][1].y, cell._vertices[0][1][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[0][1][0].x, cell._vertices[0][1][0].y, cell._vertices[0][1][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][1].x, cell._vertices[1][0][1].y, cell._vertices[1][0][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))                
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][1].x, cell._vertices[1][1][1].y, cell._vertices[1][1][1].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))                
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][1][0].x, cell._vertices[1][1][0].y, cell._vertices[1][1][0].z))
#                 f.writelines('''    VERTEX, {}, {}, {},\n'''.format(cell._vertices[1][0][0].x, cell._vertices[1][0][0].y, cell._vertices[1][0][0].z))
                
#             f.writelines('''    END
#     ]
# cmd.load_cgo({}_Grids, '{}_Grids')     
# '''.format(software, software))
            
#     def write_results(self, resultFilePath):
#         with open(resultFilePath, 'w') as f:
#             totalGrids      = self._grid._Vertices.size
#             emptyGrids      = len(self._emptyGrids)
#             partialGrids    = len(self._partialGrids)
#             fullGrids       = len(self._fullGrids)
            
#             unitVolume      = self._grid._gridSize**3
            
#             f.writelines('''
# Protein             :   {}
#     Total Grids     :   {}  Angstrom^3
#     Full Grids      :   {}  Angstrom^3
#     Partial Grids   :   {}  Angstrom^3
#     Empty Grids     :   {}  Angstrom^3
# '''.format(self._proteinName, totalGrids*unitVolume, fullGrids*unitVolume, partialGrids*unitVolume, emptyGrids*unitVolume))
            
#             for software in self._channels.keys():
#                 software_totalGrids     = len(self._channels[software]['totalGrids'])
#                 software_fullGrids      = self._channels[software]['fullGrids']
#                 software_partialGrids   = self._channels[software]['partialGrids']
#                 software_emptyGrids     = self._channels[software]['emptyGrids']
                
#                 f.writelines('''
# Software            :   {}
#     Total Grids     :   {}  Angstrom^3
#     Full Grids      :   {}  Angstrom^3
#     Partial Grids   :   {}  Angstrom^3
#     Empty Grids     :   {}  Angstrom^3
# '''.format(software, software_totalGrids*unitVolume, software_fullGrids*unitVolume, software_partialGrids*unitVolume, software_emptyGrids*unitVolume))