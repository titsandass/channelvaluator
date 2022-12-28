def test_mgos_cavity_api(pdbFilename, includeHETATM, solventRadius, gateRadius):
    import PyMGOS

    MG = PyMGOS.MolecularGeometry()

    if includeHETATM:
        print("Include HETATM")
        MG.load(pdbFilename)
        resultChannelName = pdbFilename.replace(".pdb", "_channel_MGOS_includeHETATM.py")

    else:
        print("exclude HETATM")
        MG.load_except_PDB_HETATM(pdbFilename)
        resultChannelName = pdbFilename.replace(".pdb", "_channel_MGOS_excludeHETATM.py")

    print("PDB loaded")
    MG.preprocess()

    channels = MG.compute_channels(solventRadius, gateRadius)
    MG.write_PyMOL_script(channels, resultChannelName)
    print("# channels : ", channels.number_of_channels())
