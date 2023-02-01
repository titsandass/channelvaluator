import numpy as np

class Grid:
    def __init__(self, minPt, maxPt, gridSize):
        self._gridSize  = gridSize
        
        self._minPt     = minPt
        self._maxPt     = maxPt
        
        self._bottomleft= np.floor(self._minPt / gridSize) * gridSize
        self._topright  = np.ceil( self._maxPt / gridSize) * gridSize
        
        self._numX, self._numY, self._numZ = ((self._topright - self._bottomleft) / gridSize).astype(np.int32) + 1

        self._Vertices  = np.empty(( self._numX, self._numY, self._numZ ), dtype=object)
    
    def _split(self):
        #gridVertex
        for i in range(self._numX):
            for j in range(self._numY):
                for k in range(self._numZ):
                    point = self._bottomleft + (i*self._gridSize, j*self._gridSize, k*self._gridSize)
                    self._Vertices[i][j][k] = self._GridVertex(point)

    def get_including_gridVertices_of(self, atomORSphere):
        center  = atomORSphere[0][:3]
        radii   = atomORSphere[0][3]

        bottomleft  = np.array([np.floor((i-radii)/self._gridSize) * self._gridSize for i in center], dtype=np.float16)
        topright    = np.array([np.ceil( (i+radii)/self._gridSize) * self._gridSize for i in center], dtype=np.float16)

        bottomleft_gridIndex    = np.array((bottomleft - self._bottomleft)  / self._gridSize, dtype=np.uint16)
        topright_gridIndex      = np.array((topright   - self._bottomleft)  / self._gridSize, dtype=np.uint16)

        bottomleft_gridIndex    = np.maximum(bottomleft_gridIndex   , [0, 0, 0])
        topright_gridIndex      = np.minimum(topright_gridIndex     , [self._numZ, self._numX, self._numY])

        gridVertices = self._Vertices[
            bottomleft_gridIndex[0]:topright_gridIndex[0]+1,
            bottomleft_gridIndex[1]:topright_gridIndex[1]+1,
            bottomleft_gridIndex[2]:topright_gridIndex[2]+1
        ]
        
        return gridVertices

    class _GridVertex:
        def __init__(self, point):
            self.point              = point

            self.checked_atom       = False
            self.checked_channel    = False

            self.boundary           = False
            self.interior           = True       

        def __eq__(self, other):
            return id(self) == id(other)

        def __hash__(self):
            return hash(id(self))

        def check_intersection_with_atom(self, atom):
            center  = np.round(atom[0][:3], 3)
            radii   = np.round(atom[0][3], 3)

            dist = np.linalg.norm(self.point - center)
            if dist < radii:
                self.checked_atom = True
                # self.intersectingAtoms.append(atom)

        def check_intersection_with_inflated_atom(self, inflatedAtom):
            center  = inflatedAtom[0][:3]
            radii   = inflatedAtom[0][3]

            dist = np.linalg.norm(self.point - center)
            if dist < radii:
                self.boundary = True            

        def has_intersection_with_sphere(self, sphere):
            center  = sphere[:3]
            radii   = sphere[3]

            dist = np.linalg.norm(self.point - center)
            if dist < radii:
                return True
            else:
                return False



        def set_checked_atom(self, checked:bool):
            self.checked_atom = checked
        
        def set_checked_channel(self, checked:bool):
            self.checked_channel = checked

        def set_boundary(self, checked:bool):
            self.boundary = checked

        def set_interior(self, checked:bool):
            self.interior = checked
        


        def is_checked_atom(self):
            return self.checked_atom

        def is_checked_channel(self):
            return self.checked_channel

        def is_boundary(self):
            return self.boundary

        def is_interior(self):
            return self.interior