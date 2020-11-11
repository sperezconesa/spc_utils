'''
SPC's module with visualization utils.
'''
def fasta(filename):
    def fasta_in(data=''):
        from IPython.display import display
        bundle = {}
        bundle['application/vnd.fasta.fasta'] = data
        bundle['text/plain'] = data
        display(bundle, raw=True) 
    '''
    Visualize a fasta file.
    
    Parameters
    ----------
    filename: path and filename.  
    Returns
    -------
    '''
    import os
    
    assert os.path.exists(filename), 'File does not exist.'
    extension = filename.split('.')[-1]
    assert extension == 'fasta', 'File is not fasta.'
    
    file = open(filename,'r')
    data = file.read()
    file.close()
    fasta_in(data)    

def clustal(filename):
    def clustal_in(data=''):
        from IPython.display import display
        bundle = {}
        bundle['application/vnd.clustal.clustal'] = data
        bundle['text/plain'] = data
        display(bundle, raw=True) 
    '''
    Visualize a clustal file.
    
    Parameters
    ----------
    filename: path and filename.  
    Returns
    -------
    '''
    import os
    
    assert os.path.exists(filename), 'File does not exist.'
    extension = filename.split('.')[-1]
    assert extension == 'clustal', 'File is not clustal.'
    
    file = open(filename,'r')
    data = file.read()
    file.close()
    clustal_in(data)  
def AnnotateResidues(path, ax, x_lims, fontfamily=None, fontsize=None, H_color='b', H_hatch='///', L_color=[0.7, 0.7, 0.7],
                     B_color='r', SF_color='g', alpha=1):
    """
      Annotate the x-axis of a matplotlib axes with secondary structure of residues.
      The structure is read from a txt like this KcsA example:
      _____________________
      # Type start end name
      H  26 50 TM1
      L  50 62
      H   62 73 PH
      L 73 75
      SF 75 80
      L 80 88
      H 88 121 TM2
      _______________________
      This is to label resid segment 3 10 as helix (H) with name M2 (name optional). This can be done
      with beta sheets (B), loops (L) or selectivity filters (SF).
      Parameters
      ----------
      Path : str.
          Path to the annotation txt.
      ax : plt ax object.
      lim : list.
           List with xlims of the axes.
      fontsize: int
      fontfamily: str
      H_color: str
      H_hatch: str
      L_color: str
      B_color: str
      SF_color: str
      alpha: float

      Returns
      -------
      ax:
    """
    import matplotlib.patches as mpatches
    import matplotlib.transforms as mtransforms
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams['hatch.linewidth'] = 4.0

    y0=-0.1
    if H_hatch == "":
        fill = True
    else:
        fill = False
    with open(path) as file:
        for line in file:
            if line[0] != "#":
                l=line.split()
                trans = mtransforms.blended_transform_factory(ax.transData,ax.transAxes)
                if ((x_lims[0] > int(l[1])) and ( x_lims[0]>int(l[2]))) or ((x_lims[1] < int(l[1])) and ( x_lims[1]<int(l[2]))):
                    continue
                else:
                    x_tail = max(int(l[1]), x_lims[0])
                    x_head = min(int(l[2]), x_lims[1])
                    if l[0] == "H":
                        rect_style="simple,tail_width=25"
                        rect_backgr = mpatches.FancyArrowPatch((x_tail,y0), (x_head,y0), arrowstyle=rect_style,
                                                               color='white', zorder=0, transform=trans,shrinkA=0,shrinkB=0)
                        rect_backgr.set_clip_on(False)
                        rect = mpatches.FancyArrowPatch((x_tail,y0), (x_head,y0), arrowstyle=rect_style, fill=fill,
                                                        color=H_color, hatch=H_hatch, transform=trans,
                                                        snap=True, shrinkA=0, shrinkB=0, alpha=alpha)
                        if len(l)>3:
                            ax.text((x_head-x_tail)/2+x_tail, y0, l[3], ha="center", va="center", family=fontfamily,
                                    size=fontsize,transform=trans, bbox=dict(facecolor='white', alpha=1, edgecolor=H_color))
                        rect.set_clip_on(False)
                        ax.add_patch(rect_backgr)
                        ax.add_patch(rect)
                    elif l[0] == "B":
                        arrow_style="simple,head_length=15,head_width=30,tail_width=10"
                        arrow = mpatches.FancyArrowPatch((x_tail,y0), (x_head,y0), arrowstyle=arrow_style,
                                                         transform=trans,shrinkA=0,shrinkB=0, color=B_color)
                        arrow.set_clip_on(False)
                        ax.add_patch(arrow)
                    elif l[0] == "L":
                        arrow_style="simple,tail_width=6"
                        arrow = mpatches.FancyArrowPatch((x_tail,y0), (x_head,y0), arrowstyle=arrow_style,
                                                         transform=trans,edgecolor=None,shrinkA=0,shrinkB=0, color=L_color)
                        arrow.set_clip_on(False)
                        ax.add_patch(arrow)
                    elif l[0] == "SF":
                        arrow_style='simple,tail_width=6'
                        arrow = mpatches.FancyArrowPatch((x_tail,y0), (x_head,y0), arrowstyle=arrow_style, transform=trans,
                                                         linestyle=':', edgecolor=None,  color=SF_color, shrinkA=0,
                                                         shrinkB=0, alpha=alpha)

                        ax.text((x_head-x_tail)/2+x_tail, y0+0.025, "SF", ha="center", va="center", family=fontfamily,
                                size=fontsize,transform=trans,
                                bbox=dict(facecolor='white', alpha=0, edgecolor="White"),color=SF_color)

                        arrow.set_clip_on(False)
                        ax.add_patch(arrow)

    ax.set_xlim([x_lims[0], x_lims[1]])
    return ax
