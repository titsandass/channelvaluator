import numpy as np

class Protein: 
    def __init__(self):
        self.name = None
        self.atoms = list()

    def set_name(self, name):
        self.name = name

    def add_atoms(self, atom):
        self.atoms.append(atom)

    