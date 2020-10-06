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
    substitute = {
        'CYS': 'SG',
        'SER': 'OG',
        'THR': 'CG2',
        'VAL': 'CB1',
    }
    resname = np.unique(selection.resnames)
    resid = np.unique(selection.resids)
    assert resname.shape[0] != 1 or resid.shape[0] == 1, 'more than one resid'
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