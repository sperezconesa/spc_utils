'''
This module contains many scripts to do analysis on trajectories. 
'''
def mda_janin_with_CSTV(u,selection):
    '''
    Takes in a MDAnalysis universe and a selection of that universe that contains
    one or several aminoacids sharing resid and resname. It calculates the regular
    Janin angles but if it is a CSTV it gives out only chi1 if it is a GAP it returns
    none.
      Parameters
      ----------
      u: MDAnalysis universe.
      selection: selection of the universe containing one or several aminoacids
                sharing resid and resname.
      Returns
      -------
      output_angles: numpy array of angles, shape[n_frames, n_dihedrals, n_subunits]
    '''
    from MDAnalysis.analysis.dihedrals import Janin
    from MDAnalysis import Universe
    import numpy as np

    substitute = {
        'CYS': 'SG',
        'SER': 'OG',
        'THR': 'CG2',
        'VAL': 'CG1',
    }
    resname = np.unique(selection.resnames)
    resid = np.unique(selection.resids)
    assert resname.shape[0] != 1 or resid.shape[0] == 1, 'more than one resid'
    assert isinstance(u, Universe) , 'u should be a MDAnlaysis universe.'
    if resname[0] in ['CYS', 'SER', 'THR', 'VAL']:
        my_list = []
        for res in selection.residues:
            my_list0 =[]
            for ts in u.trajectory:
                chi1 = res.chi1_selection(cg_name=substitute[resname[0]]
                                   ).dihedral.value()
                my_list0.append(chi1)
            my_list.append(my_list0)
        output_angles = np.array(my_list)
        output_angles = output_angles.reshape([output_angles.shape[0],output_angles.shape[1],1])
        output_angles = output_angles + 180.
    else:
        try:
            output_angles = Janin(selection).run().angles
        except:
            output_angles = None
    return output_angles
def get_rmsd(u, ref=None, sel_str='name CA and protein',sel_str_al=None, write_trj = False, trj_path= 'align.xtc' ):
    '''
    Takes in a MDAnalysis universe and a selection of that universe and calculate the RMSD
    with a provided reference or the  initial snapshot.
      Parameters
      ----------
      u: MDAnalysis universe.
      ref: MDAnalysis snapshot used as reference.
      sel_str: selection string of the universe containing the atoms to calculate the RMSD.
      sel_str_al: selection string of the universe containing the atoms to align the trajectory.
      write_trj: Do you want to write the trajectory.
      trj_path: Path to write trajectory.
      Returns
      -------
      rmsd: numpy of time versus RMSD.
    '''
    import MDAnalysis as mda
    import numpy as np
    import os
    from MDAnalysis.analysis import rms
    from MDAnalysis.analysis.align import AlignTraj

    assert isinstance(u, mda.Universe) , 'u should be a MDAnlaysis universe.'
    assert isinstance(ref, mda.Universe) or ref == None , 'ref should be a MDAnlaysis universe or None.'
    assert isinstance(sel_str, str), 'sel_str should be string.'
    assert isinstance(sel_str_al, str) or sel_str_al == None, 'sel_str_al should be string.'
    assert isinstance(write_trj, bool), 'write_trj should be a bool.'
    assert isinstance(trj_path, str), 'trj_path should be a str.'



    if sel_str_al == None:
        sel_str_al = sel_str
    if ref == None:
        ref_CA=u.select_atoms(sel_str)
        ref = u.copy()
    else:
        ref_CA=ref.select_atoms(sel_str)

    t = []

    if write_trj == True:
        rmsd = AlignTraj(u, ref, select=sel_str_al, filename=trj_path).run().rmsd
        u.trajectory[0]
    else:
        rmsd = AlignTraj(u, ref, select=sel_str_al).run().rmsd
        u.trajectory[0]


    for i, ts in enumerate(u.trajectory):
        t.append(ts.dt*i)

    return np.vstack([np.array(t).flatten(),rmsd]).T
def get_dssp(path,path_top,simplified=True):
    '''
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
    '''
    import mdtraj as md
    import MDAnalysis as mda
    import os
    import matplotlib.pyplot as plt
    
    assert isinstance(simplified,bool),'Simplified should be boulean.'
    assert os.path.isfile(path),'Structure path does not exisist.'
    assert os.path.isfile(path_top),'Topology path does not exisist.'
    
    traj=mda.Universe(path_top,path)
    u=md.load(path,top=path_top,atom_indices=traj.select_atoms("backbone").atoms.indices)
    dssp=md.compute_dssp(u,simplified)
    time =np.linspace(0,traj.trajectory[-1].time,traj.trajectory.n_frames)
    time = time/1000
    n=traj.select_atoms("protein").n_segments
    resid=np.split(traj.select_atoms("name CA and protein").resids,n)[0]
    if simplified:
        ss=["H","E","C"]
        ss_name=["H","B","C"]
    else:
        ss=["H","B","E","I",'T','S','']
        ss_name=[r'$\alpha$',r'$\beta$ Bridge',r'$\beta$ Ladder','3/10',r'$\pi$','Turn','Bend','Coil']     
    dssp0=np.zeros(dssp.shape)
    for i,code in enumerate(ss):
        dssp0[dssp==code]=i
    fig,ax = plt.subplots(n,1,figsize = (17,5*n))
    k=[]
    for segid,a in zip(np.split(dssp0.T,4),ax.flat):
        k.append(a.pcolor(time,resid,segid,cmap=plt.cm.get_cmap("tab10", len(ss))))
        formatter = plt.FuncFormatter(lambda val, loc: ss_name[val])
        a.set_xlabel("t (ns)")
        a.set_ylabel("resid")
    fig.colorbar(k[-1],ax=ax,ticks=range(len(ss)),format=formatter)
    return fig,ax
def distance_atoms(u,sel1,sel2):
    '''
    Calculate the distance between two atoms (sel1, sel2) as a function of time in the trajectory trj.

    Parameters
    ----------
    u: MDA universe to analyz trajectory to analyze. 
    sel1: MDA selection containing 1 atom.
    sel2: MDA selection containing 1 atom.
      
    Returns
    -------
    d: matplotlib figure object.
    '''
    from MDAnalysis.analysis.distances import dist
    from MDAnalysis import Universe
    from MDAnalysis import AtomGroup
    from numpy import array

    assert isinstance(u, Universe) , 'u should be a MDAnlaysis universe.'
    assert isinstance(sel1, Universe) , 'sel1 should be a MDAnlaysis universe.'
    assert isinstance(sel2, Universe) , 'sel2 should be a MDAnlaysis universe.'

    d = []
    g1 = u.select_atoms(sel1)
    g2 = u.select_atoms(sel2)
    for ts in u.trajectory:
        d.append([ts.time/1000,dist(g1, g2, box=ts.dimensions)[2,0]])
    return array(d)
