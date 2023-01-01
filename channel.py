class Channel:
    def __init__(self, software):
        self._software = software
        self._spheres = set()
        self._intersectingVertices = set()
        self._buriedVertices = set()

    def add_sphere(self, sphere):
        self._spheres.add(sphere)

    def add_intersecting_vertex(self, ivertex):
        self._intersectingVertices.add(ivertex)
        
    def add_buried_vertex(self, ivertex):
        self._buriedVertices.add(ivertex)