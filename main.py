def CV(pdbFileName):
    print('0.Protein {}'.format(pdbFileName))

    channelFilePaths    = {
        'MGOS'  : pdbFileName.replace('.pdb','_channel_parsed_MGOS.py'),
        'MOLE'  : pdbFileName.replace('.pdb','_channel_MOLE.py'),
        'CAVER' : pdbFileName.replace('.pdb','/tunnels/'),
        '3V'    : pdbFileName.replace('.pdb','_channel_3V.py')
    }

    start = time.time()

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

    done = time.time()

print('8.write_channel_statistics')
CV.write_channel_statistics(pdbFileName.replace('.pdb', '_CV_stats.py'), done-start)
print('9.write_result_in_PyMOL_script')
CV.write_result_in_PyMOL_script(pdbFileName.replace('.pdb', '_CV_result.py'))

    print('Done')

if __name__ == '__main__':
    import os, sys, platform
    import time

    from channelValuator import ChannelValuator

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

    includeHETATM       = False
    cutoffRatio         = 3
    initialIncrement    = 0
    step                = 0.2
    gridSize            = 0.5

    # print(*sys.argv)
    CV(sys.argv[1])