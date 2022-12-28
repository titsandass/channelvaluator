import os, sys, platform

from channelValuator import ChannelValuator

pdbFilepath = os.getcwd() + "\\data\\"
resultFilepath = os.getcwd() + "\\results\\"
if platform.system() != 'Windows':
    pdbFilepath = pdbFilepath.replace('\\', '/')
    resultFilepath = resultFilepath.replace('\\', '/')

pdbFileName = pdbFilepath+'6g2j.pdb'

includeHETATM = False
gridSize = 0.5
cutoffRatio = 5

channelFilePaths    = {
    'MGOS'  : pdbFileName.replace('.pdb','_channel_MGOS.py'),
    'MOLE'  : pdbFileName.replace('.pdb','_channel_MOLE.py'),
    'CAVER' : pdbFileName.replace('.pdb','/tunnels/'),
    '3V'    : pdbFileName.replace('.pdb','_channel_3V.py')
}

CV = ChannelValuator()

CV.set_protein(pdbFileName, includeHETATM)
CV.set_grid(gridSize)
CV.set_channels(channelFilePaths)

CV.find_ground_truth_grid(cutoffRatio)
CV.check_intersection_with_channels()
CV.validate_channel()

CV.write_results(None)

pass
