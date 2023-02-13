import os

currdir = os.getcwd() + '\\data\\'
MGOSChannels = list()
for (root, dirs, files) in os.walk(currdir):
    for file in files:
        if file.endswith('channel_MGOS.py'):
            MGOSChannels.append(currdir+file)


for channelFile in MGOSChannels:
    lines = None

    deleteIndices = list()
    with open(channelFile, 'r') as f:
        lines = f.readlines()
        delete = False
        for i, line in enumerate(lines):
            if line.startswith('contribAtom'):
                deleteIndices.append(i)
                delete = True
            elif line.startswith('cmd.load_cgo(spine'):
                deleteIndices.append(i)
                delete = False
            elif delete == True:
                deleteIndices.append(i)
    
        for idx in sorted(deleteIndices, reverse=True):
            del lines[idx]

    with open(channelFile.replace('_MGOS.py','_parsed_MGOS.py'), 'w') as f:
        f.writelines(lines)