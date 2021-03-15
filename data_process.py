"""
SPC's scripts to process gmx trajectories using gmxapi.
"""


def process_trajectory(
    files,
    path=".",
    begin=0,
    end=1.0e20,
    skip=1,
    output_group="all",
    center_group="protein",
    align=False,
    name_base="",
    name=None,
    ndx_path="index.ndx",
    tpr_path="topol.tpr",
    cat_overwrite=True,
):

    import os

    import gmxapi as gmx

    def did_it_run(command):
        did_it_run = True
        if command.output.returncode.result() != 0:
            print(command.output.erroroutput.result())
            exit
            did_it_run = False
        return did_it_run

    """
    Center structure, align structure, do pbc, skip frames and reduce output of
    given structure file or files.
      Parameters
      ----------
      files: structure file or files.(xtc,trr,gro,pdb,cpt)
      path: path to structure file/s.
      skip: skip every skip frame.
      output_group: output group of index.ndx.
      center_group: center/align group of index.ndx.
      name_base: use this name to start the name output.
      name: use this name for the output
      ndx_path: index file path.
      tpr_path: tpr file path.
      Returns
      -------
    """
    # Asertions on input
    assert os.path.exists(path), "Path does not exists."
    assert isinstance(files, str) or isinstance(
        files, list
    ), "File/s should be a string or list"
    if isinstance(files, str):
        assert os.path.isfile(path + f"/{files}"), (
            "File" + files + " does not exist."
        )
    else:
        for f in files:
            assert isinstance(f, str), "File content of files must be str."
            assert os.path.isfile(path + f"/{f}"), (
                "File" + f + " does not exist."
            )

    for f in [ndx_path, tpr_path]:
        assert os.path.isfile(path + f"/{f}"), "File" + f + " does not exist."
    assert isinstance(skip, int), "Skip should be int."

    for st in [path, output_group, center_group]:
        assert isinstance(st, str), f"{st} should be a string."
    if name is not None:
        assert isinstance(name, str), f"{name} should be a string"
    assert isinstance(name_base, str), f"{name_base} should be a string"

    cdminus = os.getcwd()
    os.chdir(path)

    if isinstance(files, list):
        if cat_overwrite:
            args = ["trjcat"]
        else:
            args = ["trjcat -nooverwrite"]

        catcomm = gmx.commandline_operation(
            "gmx",
            arguments=["trjcat"],
            input_files={"-f": files},
            output_files={"-o": "cat.xtc"},
        )
        catcomm.run()
        assert did_it_run(catcomm), "Concating failed failed"
        file = "cat.xtc"
        extension = "xtc"
        print(files)
    else:
        base, extension = files.split(".")
        file = files

    if extension == "cpt":
        extension = "gro"

        # Output name
    if name is None:
        if name_base != "":
            name_base = name_base + "_"
        if extension == "gro" or extension == "pdb":
            name = f"{name_base}{output_group}_pbc.{extension}"
            name_al = f"{name_base}{output_group}_pbc_al.{extension}"
        else:
            name = f"{name_base}{output_group}_sk{skip}_pbc.{extension}"
            name_al = f"{name_base}{output_group}_sk{skip}_pbc_al.{extension}"
    else:
        name_al = name
    begin_end = ["-b", str(begin), "-e", str(end)]
    for f in [name, name_al, "kk.xtc"]:
        if os.path.isfile(f):
            os.remove(f)

    # Cluster pbc
    trjconv0 = gmx.commandline_operation(
        "gmx",
        arguments=["trjconv", "-skip", str(skip), "-pbc", "cluster"]
        + begin_end,
        input_files={"-f": file, "-n": ndx_path, "-s": tpr_path},
        stdin=f"{center_group} {output_group}",
        output_files={"-o": f"kk.{extension}"},
    )
    trjconv0.run()
    assert did_it_run(trjconv0), "Clustering failed"

    # Center and pbc
    trjconv1 = gmx.commandline_operation(
        "gmx",
        arguments=["trjconv", "-pbc", "mol", "-center"],
        input_files={"-f": f"kk.{extension}", "-n": ndx_path, "-s": tpr_path},
        stdin=f"{center_group} {output_group}",
        output_files={"-o": name},
    )
    trjconv1.run()
    os.remove(f"kk.{extension}")
    assert did_it_run(trjconv1), "Pbc failed"

    # Align
    if align:
        trjconv2 = gmx.commandline_operation(
            "gmx",
            arguments=["trjconv", "-fit", "rot+trans"],
            input_files={"-f": name, "-n": ndx_path, "-s": tpr_path},
            stdin=f"{center_group} {output_group}",
            output_files={"-o": name_al},
        )
        trjconv2.run()
        assert did_it_run(trjconv2), "Align failed"

    os.chdir(cdminus)
    return


def get_dssp(path, path_top, simplified=True):
    """
    Calculate and plot the secondary structure as a function of time and split per subunit.

    Parameters
    ----------
    path: path to structure file.
    path_top: path to topology file
    simplified: Should the structure types be simplified.

    Returns
    -------
    fig: matplotlib figure object.
    ax: axes of the figures.
    """
    import os

    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    import mdtraj as md

    assert isinstance(simplified, bool), "Simplified should be boulean."
    assert os.path.isfile(path), "Structure path does not exisist."
    assert os.path.isfile(path_top), "Topology path does not exisist."

    traj = mda.Universe(path_top, path)
    u = md.load(
        path,
        top=path_top,
        atom_indices=traj.select_atoms("backbone").atoms.indices,
    )
    dssp = md.compute_dssp(u, simplified)
    time = np.linspace(0, traj.trajectory[-1].time, traj.trajectory.n_frames)
    time = time / 1000
    n = traj.select_atoms("protein").n_segments
    resid = np.split(traj.select_atoms("name CA and protein").resids, n)[0]
    if simplified:
        ss = ["H", "E", "C"]
        ss_name = ["H", "B", "C"]
    else:
        ss = ["H", "B", "E", "I", "T", "S", ""]
        ss_name = [
            r"$\alpha$",
            r"$\beta$ Bridge",
            r"$\beta$ Ladder",
            "3/10",
            r"$\pi$",
            "Turn",
            "Bend",
            "Coil",
        ]
    dssp0 = np.zeros(dssp.shape)
    for i, code in enumerate(ss):
        dssp0[dssp == code] = i
    fig, ax = plt.subplots(n, 1, figsize=(17, 5 * n))
    k = []
    for segid, a in zip(np.split(dssp0.T, 4), ax.flat):
        k.append(
            a.pcolor(
                time, resid, segid, cmap=plt.cm.get_cmap("tab10", len(ss))
            )
        )
        formatter = plt.FuncFormatter(lambda val, loc: ss_name[val])
        a.set_xlabel("t (ns)")
        a.set_ylabel("resid")
    fig.colorbar(k[-1], ax=ax, ticks=range(len(ss)), format=formatter)
    return fig, ax
