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