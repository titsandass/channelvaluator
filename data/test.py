from pymol.cgo import *
from pymol import cmd

view = cmd.get_view()
a = [
    COLOR, 1.000000, 0.000000, 0.000000,
    SPHERE, 5.899, -38.704, -35.697, 0.5,
    SPHERE, 5.899, -38.704, 6.806, 0.5,
    SPHERE, 5.899, 0.287, -35.697, 0.5,
    SPHERE, 5.899, 0.287, 6.806, 0.5,
    
    COLOR, 1.000000, 1.000000, .000000,
    SPHERE, 45.225, -38.704, -35.697, 0.5,
    SPHERE, 45.225, -38.704, 6.806, 0.5,
    SPHERE, 45.225, 0.287, -35.697, 0.5,
    SPHERE, 45.225, 0.287, 6.806, 0.5,
    
    COLOR, 1.000000, 1.000000, 1.00000,
    SPHERE, 45.225, -38.704, 6.806, 0.6,
    
    COLOR, 1.000000, 1.000000, 1.00000,
    SPHERE, 45.225, 0.287, 6.806, 0.6,
    
    COLOR, 0.000000, 1.000000, 1.00000,
    SPHERE, 45.225, 0.287, 6.806, 0.7,
]
cmd.load_cgo(a, 'a')

cmd.fetch('2OBI')
cmd.remove('HETATM')

cmd.set_view(view)