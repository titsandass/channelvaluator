import os, sys, platform, time, pickle

libPath         = os.getcwd() + "\\lib\\"
pdbFilepath     = os.getcwd() + "\\data\\"
resultFilepath  = os.getcwd() + "\\channels\\"
logFilepath     = os.getcwd() + "\\channels\\"
if platform.system() != 'Windows':
    libPath = libPath.replace('\\', '/')
    pdbFilepath = pdbFilepath.replace('\\', '/')
    resultFilepath = resultFilepath.replace('\\', '/')
    logFilepath = logFilepath.replace('\\', '/')
sys.path.insert(0, libPath)

import PyMGOS

pdbFiles = list()
for (root, dirs, files) in os.walk(pdbFilepath):
    for file in files:
        if file.endswith('.pdb'):
            pdbFiles.append(root+file)
    break

# pdbFiles = ['C:\\Users\\hwkim\\Desktop\\workspaces\\ChannelValuator\\data\\4bpc.pdb']

includeHETATM       = False   
solventProbeRadius  = 1.4
gateProbeRadius     = 3

calc_times = list()
for pdbFile in pdbFiles:
    # if os.path.exists(pdbFile.replace('.pdb', '.a.qtf')):
    #     os.remove(pdbFile.replace('.pdb', '.a.qtf'))
    
    start = time.time()

    MG = PyMGOS.MolecularGeometry()
    if includeHETATM:
        MG.load(pdbFile)
    else:
        MG.load_except_PDB_HETATM(pdbFile)            
    MG.preprocess()

    preprocess = time.time()

    channels = MG.compute_channels(solventProbeRadius, gateProbeRadius)

    compute_channel = time.time()

    MG.write_PyMOL_script(channels, pdbFile.replace('.pdb', '_channel_MGOS.py'))

    if platform.system() != 'Windows':
        calc_times.append({'pdb' : pdbFile.split('/')[-1].replace('.pdb',''), 'preprocess' : preprocess - start, 'compute_channel' : compute_channel - preprocess})
        print('{} Done'.format(pdbFile.split('/')[-1].replace('.pdb','')))
    else:
        calc_times.append({'pdb' : pdbFile.split('\\')[-1].replace('.pdb',''), 'preprocess' : preprocess - start, 'compute_channel' : compute_channel - preprocess})
        print('{} Done'.format(pdbFile.split('\\')[-1].replace('.pdb','')))

with open(logFilepath+'MGOS_calctime.pickle', 'wb') as f:
    pickle.dump(calc_times, f)

# MolecularChannelSet MG.compute_channels(const double& solventProbeRadius, const double& gateSize)
# MolecularVoidSet    MG.compute_voids_of_Lee_Richards_model(const double& solventProbeRadius)
# MolecularPocketSet  MG.compute_pockets(const double& ligandSize, const double& solventProbeRadius)

# void write_PyMOL_script(MolecularChannelSet& channels, const string& fileName)
