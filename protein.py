class Protein: 
    def __init__(self, name=None):
        self.name = name
        self.atoms = list()
        # np.array([x,y,z,r], dtype=np.float32)
        # (atom, atomName, resName, atomNum)

    def set_name(self, name):
        self.name = name

    def add_atoms(self, atom):
        self.atoms.append(atom)

    