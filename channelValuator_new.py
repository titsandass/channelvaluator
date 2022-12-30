from protein import Protein

import numpy as np
import os

class ChannelValuator:
    def __init__(self):
        self._pdbFileName       = None
        self._protein           = None
        self._inflatedAtoms     = None
        self._includeHETATM     = False
        self._channelFilePaths  = None
        self._channels          = dict()

        self._min               = np.array([np.inf, np.inf, np.inf]   , dtype=np.float16)
        self._max               = np.array([-np.inf, -np.inf, -np.inf], dtype=np.float16)

###################################################################    
    def set_protein(self, pdbFileName, includeHETATM):
        if not os.path.exists(pdbFileName):
            print('File Not Exists : {}'.format(pdbFileName))
            return

        self._pdbFileName = pdbFileName
        self._includeHETATM = includeHETATM
        proteinName = pdbFileName.split('/')[-1].replace('.pdb','')
        self._protein = Protein(proteinName)

        with open(pdbFileName, 'r') as pdbFile:
            lines = pdbFile.readlines()
            if includeHETATM:
                for line in lines:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atomInfo = self._parse_line_into_atom(line)
                        self._protein.add_atoms(atomInfo)
            else:
                for line in lines:
                    if line.startswith("ATOM"):
                        atomInfo = self._parse_line_into_atom(line)
                        self._protein.add_atoms(atomInfo)

    def _parse_line_into_atom(self, line):
        from VDWradius import VDWRadius

        atomName= line[-3:]
        atomNum = int(line[6:12])
        resName = line[17:20]

        x       = np.float32(line[30:38])
        y       = np.float32(line[38:46])
        z       = np.float32(line[46:52])

        try:
            r = VDWRadius[atomName.split()[-1].rjust(2,' ')]
        except:
            parsedAtomName = ''
            for char in atomName:
                if char.isalpha():  
                    parsedAtomName += char
            r = VDWRadius[parsedAtomName.split()[-1].rjust(2,' ')]

        atom    = np.array([x,y,z,r], dtype=np.float32)
        return (atom, atomName, resName, atomNum)

###################################################################
    def set_channels(self, channelFilePaths):
        from channel import Channel

        self._channelFilePaths = channelFilePaths
        for software in channelFilePaths.keys():
            channelFilePath = channelFilePaths[software]
            if os.path.exists(channelFilePath):
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


###################################################################
    def inflate_protein(self, cutoffRatio, initialIncrement, step):
        import sys, platform
        libPath = os.getcwd() + "\\lib\\"
        if platform.system() != 'Windows':
            libPath = libPath.replace('\\', '/')
        sys.path.insert(0, libPath)
        import PyMGOS
        
        if os.path.exists(self._pdbFileName.replace('.pdb', '.a.qtf')):
            os.remove(self._pdbFileName.replace('.pdb', '.a.qtf'))

        MG = PyMGOS.MolecularGeometry()
        if self._includeHETATM:
            MG.load(self._pdbFileName)
        else:
            MG.load_except_PDB_HETATM(self._pdbFileName)            
        MG.preprocess()

        numAllAtoms = MG.get_all_atoms().size()
        LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(initialIncrement)
        while LRboudnaryAtomSet.size() > (numAllAtoms/cutoffRatio):
            initialIncrement += step
            LRboudnaryAtomSet = MG.find_boundary_atoms_in_Lee_Richards_model(initialIncrement)
        print(numAllAtoms, LRboudnaryAtomSet.size(), initialIncrement)

        self._inflatedAtoms = Protein()
        for atom in LRboudnaryAtomSet:
            sphere = atom.geometry()
            x = sphere.center().x()
            y = sphere.center().y()
            z = sphere.center().z()
            r = sphere.radius() + initialIncrement
            atom = np.array([x,y,z,r], dtype=np.float16)

            self._inflatedAtoms.add_atoms((atom,))

            self._min = np.minimum(atom[:3] - r, self._min)
            self._max = np.maximum(atom[:3] + r, self._max)

###################################################################        
    def set_grid(self, gridSize):
        from grid import Grid

        self._gridSize = gridSize
        self._grid = Grid(self._min, self._max, gridSize)
        self._grid._split()

###################################################################
    def set_ground_truth(self):
        self._set_boundary_gridVertices()
        self._determine_interior_gridVertices()

    def _set_boundary_gridVertices(self):
        for inflatedAtom in self._inflatedAtoms:
            includingGridVertices = self._grid.get_including_gridVertices_of(inflatedAtom)
            for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                vertex = ivertex.item()
                if vertex.is_boundary():
                    continue
                vertex.check_intersection_with_inflated_atom(inflatedAtom)

    def _determine_interior_gridVertices(self):
        # x_planes = [:, 0, :], [ :, -1,  :]
        # y_planes = [:, :, 0], [ :,  :, -1]
        # z_planes = [0, :, :], [-1,  :,  :]

        imax, jmax, kmax = self._grid._Vertices.shape

        #x_planes
        self.__propagate_plane(range(imax), range(kmax), range(jmax))
        self.__propagate_plane(range(imax), range(kmax), reversed(range(jmax)))

        #y_planes
        self.__propagate_plane(range(imax), range(jmax), range(kmax))
        self.__propagate_plane(range(imax), range(jmax), reversed(range(kmax)))

        #z_planes
        self.__propagate_plane(range(jmax), range(kmax), range(imax))
        self.__propagate_plane(range(jmax), range(kmax), reversed(range(imax)))

    def __propagate_plane(self, i_range, j_range, k_range):
        for i in i_range:
            for j in j_range:
                for k in k_range:
                    vertex = self._grid._Vertices[i,j,k].item()
                    if vertex.is_boundary():
                        break
                    else:
                        if not vertex.is_interior():
                            continue
                        else:
                            vertex.set_interior(False)

###################################################################
    def verify_atom_overlapping_vertices(self):
        for atomInfo in self._protein.atoms:
            includingGridVertices = self._grid.get_including_gridVertices_of(atomInfo)
            for ivertex in np.nditer(includingGridVertices, flags=['refs_ok']):
                vertex = ivertex.item()
                if not vertex.is_boundary():
                    continue

                if vertex.is_checked_atom():
                    continue
                else:
                    vertex.check_intersection_with_atom(atomInfo)

###################################################################
    def verify_channel_overlapping_vertices(self):
        for software in self._channels.keys():
            for sphere in self._channels[software]['channel']._spheres:
                includingGridVertices = self._grid.get_including_gridVertices_of((sphere,))
                
                if includingGridVertices.size == 0:
                    raise NotImplementedError("Channel Sphere is not in the Grid")
                
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


