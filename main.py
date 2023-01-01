import os, sys, platform

from channelValuator import ChannelValuator

pdbFilepath = os.getcwd() + "\\data\\"
resultFilepath = os.getcwd() + "\\results\\"
if platform.system() != 'Windows':
    pdbFilepath = pdbFilepath.replace('\\', '/')
    resultFilepath = resultFilepath.replace('\\', '/')

pdbFileName = pdbFilepath+'2OBI.pdb'

channelFilePaths    = {
    'MGOS'  : pdbFileName.replace('.pdb','_channel_MGOS.py'),
    'MOLE'  : pdbFileName.replace('.pdb','_channel_MOLE.py'),
    'CAVER' : pdbFileName.replace('.pdb','/tunnels/'),
    '3V'    : pdbFileName.replace('.pdb','_channel_3V.py')
}

includeHETATM       = False
cutoffRatio         = 3
initialIncrement    = 0
step                = 0.5
gridSize            = 0.5

CV = ChannelValuator()

print('1.set_protein')
CV.set_protein(pdbFileName, includeHETATM)
print('2.set_channels')
CV.set_channels(channelFilePaths)

print('3.inflate_protein')
CV.inflate_protein(cutoffRatio, initialIncrement, step)
print('4.set_grid')
CV.set_grid(gridSize)

print('5.set_ground_truth')
CV.set_ground_truth()
print('6.verify_atom_overlapping_vertices')
CV.verify_atom_overlapping_vertices()
print('7.verify_channel_overlapping_vertices')
CV.verify_channel_overlapping_vertices()

print('8.write_result_in_PyMOL_script')
CV.write_result_in_PyMOL_script(pdbFileName.replace('.pdb', '_CV_result.py'))
print('Done')