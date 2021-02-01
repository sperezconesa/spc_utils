'''
SPC's scripts to process gmx trajectories using gmxapi.
'''


def process_trajectory(file, path='.', begin=0,end=1.e+20,skip=1, output_group='all', center_group='protein', align=False, name_base='', name=None,  ndx_path='index.ndx'):
    
    
    import gmxapi as gmx
    import os
    def did_it_run(command):
        did_it_run = True
        if command.output.returncode.result() != 0:
            print(command.output.erroroutput.result())
            exit
            did_it_run = False
        return did_it_run

    '''
    Center structure, align structure, do pbc, skip frames and reduce output of
    given structure file.
      Parameters
      ----------
      path: path to structure file.
      file: structure file.(xtc,trr,gro,pdb,cpt)
      skip: skip every skip frame.
      output_group: output group of index.ndx.
      center_group: center/align group of index.ndx.
      name_base: use this name to start the name output.
      name: use this name for the output
      ndx_path: index file path. 
      Returns
      -------
    '''
    # Asertions on input
    assert os.path.exists(path), 'Path does not exists.'
    for f in [file, ndx_path, 'topol.tpr']:
        assert os.path.isfile(path + f'/{f}'), f'File {f} does not exist.'
    assert isinstance(skip, int), 'Skip should be int.'

    for st in [path,file,output_group,center_group]:
        assert isinstance(st, str), f'{st} should be a string.'
    if name != None:
        assert isinstance(name, str), f'{name} should be a string'
    assert isinstance(name_base, str), f'{name_base} should be a string'
    

    

    cdminus = os.getcwd()
    os.chdir(path)

    base, extension = file.split('.')
    
    if extension == 'cpt':
        extension = 'gro'
    
        # Output name
    if name == None:
        if name_base != '':
            name_base = name_base + '_'
        if extension == 'gro' or extension == 'pdb':
            name = f'{name_base}{output_group}_pbc.{extension}'
            name_al = f'{name_base}{output_group}_pbc_al.{extension}'
        else:
            name = f'{name_base}{output_group}_sk{skip}_pbc.{extension}'
            name_al = f'{name_base}{output_group}_sk{skip}_pbc_al.{extension}'
    else:
        name_al = name
    begin_end = ['-b', str(begin), '-e', str(end)]
    for files in [name, name_al, 'kk.xtc']:
        if os.path.isfile(files):
            os.remove(files)

    #Cluster pbc
    trjconv0 = gmx.commandline_operation('gmx',arguments=['trjconv', '-skip', str(skip), '-pbc', 'cluster']+begin_end,
                                    input_files = {'-f': file,'-n': ndx_path},
                                    stdin = f"{center_group} {output_group}",
                                    output_files = {'-o': f'kk.{extension}'})
    trjconv0.run()
    assert did_it_run(trjconv0),'Clustering failed' 

    #Center and pbc
    trjconv1 = gmx.commandline_operation('gmx',arguments=['trjconv', '-pbc', 'mol','-center'],
                                    input_files={'-f': f'kk.{extension}','-n': ndx_path},
                                    stdin = f'{center_group} {output_group}',
                                    output_files = {'-o': name})
    trjconv1.run()
    os.remove(f'kk.{extension}')
    assert did_it_run(trjconv1), 'Pbc failed'

    #Align
    if align:
        trjconv2 = gmx.commandline_operation('gmx',arguments=['trjconv', '-fit', 'rot+trans'],
                                        input_files = {'-f': name,'-n': ndx_path},
                                        stdin = f'{center_group} {output_group}',
                                        output_files = {'-o': name_al})
        trjconv2.run()
        assert did_it_run(trjconv2), 'Align failed'

    os.chdir(cdminus)
    return

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
