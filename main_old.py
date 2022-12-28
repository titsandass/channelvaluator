from channelValuator import ChannelValuator

import time
import os

oneFile = True

proteinFiles    = []
if oneFile:
    proteinFiles = ["./data/2OBI.pdb"]
else:
    proteinFilesDir     = "./data/pdb_files/"
    for (root, dirs, files) in os.walk(proteinFilesDir):
        for filename in  files:
            if filename.endswith('.pdb'):
                proteinFiles.append(proteinFilesDir + filename)
        break

done = ["./data/pdb_files/6g2j.pdb",
        "./data/pdb_files/6qvc.pdb",
        "./data/pdb_files/6tjv.pdb",
        "./data/pdb_files/6vx7.pdb",
        "./data/pdb_files/7rpi.pdb",
        "./data/pdb_files/7ose.pdb",
        "./data/pdb_files/6zkc.pdb",
        ]
dead = [
        # "./data/pdb_files/7dgz.pdb"
        ]

for proteinFilePath in reversed(proteinFiles):
    if proteinFilePath in done or proteinFilePath in dead:
            continue

    channelMGOSFilePath = proteinFilePath.replace('.pdb','_channel_MGOS.py')
    channelMOLEFilePath = proteinFilePath.replace('.pdb','_channel_MOLE.py')
    channelCAVERFilePath= proteinFilePath.replace('.pdb','/tunnels/')
    channel3VFilePath   = proteinFilePath.replace('.pdb','_channel_3V.py')

    channelFilePaths    = {
        'MGOS'  : channelMGOSFilePath,
        'MOLE'  : channelMOLEFilePath,
        'CAVER' : channelCAVERFilePath,
        '3V'    : channel3VFilePath
    }

    proteinFilename     = proteinFilePath.split('/')[-1].replace('.pdb','')
    logFilePath         = './results/{}_log.txt'.format(proteinFilename)
    resultDir           = './results/'

    gridSize            = 0.3
    includeWater        = False
    includeHETATM       = False

    with open(logFilePath, 'w') as f:
        f.write('''grid Size	                        :	{}\n'''.format(gridSize))
        print('''{}'''.format(proteinFilename))
        #preprocess
        CV = ChannelValuator()

        print('''set_protein''')
        now = time.time()
        CV.set_protein(proteinFilePath, includeHETATM)

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''set_protein                        :	{} s'''.format(dt))
        f.write('''set_protein                          :	{} s\n'''.format(dt))

        print('''set_grid''')
        now = time.time()
        CV.set_grid(gridSize)

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''set_grid                           :	{} s'''.format(dt))
        f.write('''set_grid	                            :	{} s\n'''.format(dt))

        print('''set_channels''')
        now = time.time()
        CV.set_channels(channelFilePaths)

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''set_channels                       :	{} s'''.format(dt))
        f.write('''set_channels                         :	{} s\n'''.format(dt))        

        #intersection check
        print('''check_intersection_with_atoms''')
        now = time.time()
        CV.check_intersection_with_atoms()

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''check_intersection_with_atoms      :	{} s'''.format(dt))
        f.write('''check_intersection_with_atoms        :	{} s\n'''.format(dt))        

        print('''check_intersection_with_channels''')
        now = time.time()
        CV.check_intersection_with_channels()

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''check_intersection_with_channels   :	{} s'''.format(dt))
        f.write('''check_intersection_with_channels     :	{} s\n'''.format(dt))        

        #validation
        print('''validate_channel''')
        now = time.time()
        CV.validate_channel(includeWater)

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''validate_channel                   :	{} s'''.format(dt))
        f.write('''validate_channel                     :	{} s\n'''.format(dt))        

        # write
        print('''write_results''')
        now = time.time()
        CV.write_results('{}{}_result.txt'.format(resultDir,proteinFilename))

        dt = time.time() - now
        print ("\033[A                             \033[A")
        print('''write_results                      :	{} s'''.format(dt))
        f.write('''write_results                        :	{} s\n'''.format(dt))
        
        # print('''write_PyMOL_script''')
        # now = time.time()
        # CV.write_PyMOL_script('{}{}_vertices.py'.format(resultDir,proteinFilename), includeWater)

        # dt = time.time() - now
        # print ("\033[A                             \033[A")
        # print('''write_PyMOL_script                  :	{} s'''.format(dt))
        # f.write('''write_PyMOL_script                    :	{} s\n'''.format(dt))

