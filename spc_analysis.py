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
def get_rmsd(u, ref=None, sel_str='name CA and protein', skip=1,sel_str_al=None, in_memory = False ):
    '''
    Takes in a MDAnalysis universe and a selection of that universe and calculate the RMSD
    with a provided reference or the  initial snapshot.
      Parameters
      ----------
      u: MDAnalysis universe.
      ref: MDAnalysis snapshot used as reference.
      sel_str: selection string of the universe containing the atoms to calculate the RMSD.
      sel_str_al: selection string of the universe containing the atoms to align the trajectory.
      in_memory: Do you want to print a temp.xtc or do everything in memory.
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
    assert isinstance(skip, int) and skip > 0 , 'Skip should be int.'
    assert isinstance(sel_str_al, str) or sel_str_al == None, 'sel_str_al should be string.'
    assert isinstance(in_memory, bool), 'in_memory should be a bool.'

    if in_memory:
        print('This is not yet implemented')
        return

    if sel_str_al == None:
        sel_str_al = sel_str
    if ref == None:
        ref_CA=u.select_atoms(sel_str)
        ref = u.copy()
    else:
        ref_CA=ref.select_atoms(sel_str)
    l0 = []
    t0 = []
    
    if in_memory == True:
        a = AlignTraj(u, ref, select=sel_str_al,in_memory=True).run()
    else:
        AlignTraj(u, ref, select=sel_str_al, filename='tmp.xtc').run()
        u.trajectory[0]
        u.select_atoms('all').write('tmp.pdb')
        u = mda.Universe('tmp.pdb','tmp.xtc')
        
    trj_CA=u.select_atoms(sel_str)
    for i, ts in enumerate(u.trajectory[::skip]):
        l0.append([ts.dt*i,rms.rmsd(trj_CA.positions,ref_CA.positions,superposition=False)])
        t0.append([ts.dt])
    rmsd = np.array(l0)
    if in_memory == False:
        os.remove('tmp.pdb')
        os.remove('tmp.xtc')
    return rmsd,np.array(t0)
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
