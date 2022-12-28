import numpy as np

class Grid:
    def __init__(self, minPt, maxPt, gridSize):
        self._gridSize  = gridSize
        
        self._minPt     = np.round(minPt, 4)
        self._maxPt     = np.round(maxPt, 4)
        
        self._bottomleft= np.round(np.floor(self._minPt / gridSize) * gridSize, 2)
        self._topright  = np.round(np.ceil(self._maxPt / gridSize) * gridSize, 2)
        
        self._width, self._depth, self._height = np.round(self._topright - self._bottomleft, 2)

        self._numWidth  = np.int32(self._width / gridSize)
        self._numDepth  = np.int32(self._depth / gridSize)
        self._numHeight = np.int32(self._height / gridSize)

        self._Vertices  = np.empty((self._numWidth+1, self._numDepth+1  , self._numHeight+1 ),  dtype=object)
    

    def _split(self):
        #gridVertex
        for i in range(self._numWidth+1):
            for j in range(self._numDepth+1):
                for k in range(self._numHeight+1):
                    point = self._bottomleft + (i*self._gridSize, j*self._gridSize, k*self._gridSize)
                    self._Vertices[i][j][k] = self._GridVertex(np.round(point,2))

    def get_including_gridVertices_of(self, atomORSphere):
        center  = atomORSphere[0][:3]
        radii   = atomORSphere[0][3]

        bottomleft  = np.array([np.floor((i-radii)/self._gridSize)  *self._gridSize for i in center], dtype=np.float32)
        topright    = np.array([np.ceil((i+radii)/self._gridSize)   *self._gridSize for i in center], dtype=np.float32)

        bottomleft_gridIndex    = np.array((bottomleft - self._bottomleft)  / self._gridSize, dtype=np.uint16)
        topright_gridIndex      = np.array((topright - self._bottomleft)    / self._gridSize, dtype=np.uint16)

        bottomleft_gridIndex    = np.maximum(bottomleft_gridIndex   , [0, 0, 0])
        topright_gridIndex      = np.minimum(topright_gridIndex     , [self._numWidth+1, self._numDepth+1, self._numHeight+1])

        gridVertices = self._Vertices[
            bottomleft_gridIndex[0]:topright_gridIndex[0]+1,
            bottomleft_gridIndex[1]:topright_gridIndex[1]+1,
            bottomleft_gridIndex[2]:topright_gridIndex[2]+1
        ]
        
        return gridVertices

    class _GridVertex:
        def __init__(self, point):
            self.point              = point
            self.checked            = False
            self.checked_channel    = False

            self.intersectingAtoms  = list()
            self.intersectingWaters = list()
            
        @property
        def x(self):
            return self.point[0]

        @property
        def y(self):
            return self.point[1]
            
        @property
        def z(self):
            return self.point[2]        

        def __eq__(self, other):
            return id(self) == id(other)

        def __hash__(self):
            return hash(id(self))

        def check_intersection_with_atom(self, atom, isWater):
            center  = np.round(atom[0][:3], 3)
            radii   = np.round(atom[0][3], 3)

            dist = np.linalg.norm(self.point - center)
            if dist < radii:
                self.checked = True
                if isWater:
                    self.intersectingWaters.append(atom)
                else:
                    self.intersectingAtoms.append(atom)

        def has_intersection_with_sphere(self, sphere):
            center  = np.round(sphere[:3], 4)
            radii   = np.round(sphere[3], 3)

            dist = np.linalg.norm(self.point - center)            
            if np.absolute(dist - radii) < 1e-6:
                return 'Boundary'
            elif dist < radii:
                return 'Inside'
            else:
                return 'Outside'

        def set_checked(self, checked:bool):
            self.checked = checked
        
        def set_checked_channel(self, checked:bool):
            self.checked_channel = checked
        
        def is_checked(self):
            return self.checked

        def is_checked_channel(self):
            return self.checked_channel