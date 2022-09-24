import sys, os
import copy as cp
from pmx import *
from pmx import geometry
from pmx import ndx
from pmx.options import *
import numpy as np

# put these commands into pmx geometry.py
#def center_vector( v ):
#    vout = _p.center_vec( v )
#    return( vout )

#def calc_fit_R( cs1, cs2, m ):
#    R = _p.calc_fit_R(cs1, cs2, m)
#    return(R)
#

def find_last_protein_atom( m ):
    prev = ''
    chID = 0
    resID = 0
    for atom in m.atoms:
        if ('SOL' in atom.resname) or ('HOH' in atom.resname) \
           or ('NA' in atom.resname) or ('CL' in atom.resname) \
           or ('Na' in atom.resname) or ('Cl' in atom.resname) \
           or ('Mg' in atom.resname) or ('MG' in atom.resname) \
           or ('Zn' in atom.resname) or ('ZN' in atom.resname):
           chID = prev.chain_id
           resID = prev.resnr
           break
        prev = atom
    return(chID,resID)

def get_cs( m, ndx=[] ):
    if len(ndx)>0:
        #in python 3 map is an iterable, not a list. The evaluation of the lambda is defered.
        #using list() to force evaluation as a work around
        cs = list(map(lambda i: m.atoms[i].x, ndx))
    else:
        cs = m.coords()
    return(cs)


def fit(model1, model2, ndx1=[],ndx2=[] ):

    cs1 = get_cs( model1, ndx1 )
    cs2 = get_cs( model2, ndx2 )

    m = list(map(lambda x: 1., cs1)) # dummy array
    assert( len(cs1) == len(cs2) )

    # 1. remove COM
    v1 = geometry.center_vector( cs1 )
    v2 = geometry.center_vector( cs2 )
    model1.translate( [-v1[0], -v1[1], -v1[2] ] )
    model2.translate( [-v2[0], -v2[1], -v2[2] ] )
    cs1 = get_cs( model1, ndx1 )
    cs2 = get_cs( model2, ndx2 )

    # 2. rotate
    R = geometry.calc_fit_R(cs1, cs2, m)
    geometry.apply_fit_R( model2.atoms, R)

    # 3. translate back to the model1 position
    model1.translate( v1 )
    model2.translate( v1 )

    return(v1,v2,R)

def rotate_coords_R( atoms, R):
    for atom in atoms:
        x_old = list(map(lambda x: x, atom.x))
        for r in range(3):
            atom.x[r] = 0
            for c in range(3):
                atom.x[r]+=R[r][c]*x_old[c]

def rotate_velocities_R( m, R):
    for atom in m.atoms:
        v_old = list(map(lambda v: v, atom.v))
        for r in range(3):
            atom.v[r] = 0
            for c in range(3):
                atom.v[r]+=R[r][c]*v_old[c]

def select_ndx( fname="index.ndx", message=False ):
    if not os.path.isfile( fname ):
        sys.stdout.write("Could not find index file: %s\n" % fname )
        return(False)
    ndx_file = ndx.IndexFile( fname )
    names = ndx_file.names
    i = 0
    ndxDict = {}
    for name in names:
        atomNum = len(ndx_file[name].ids)
        sys.stdout.write('%d %s: %d atoms\n' % (i,name,atomNum) )
        ndxDict[i] = ndx_file[name].ids
        i+=1
    sys.stdout.write('\n')

    if message==False:
        sys.stdout.write('Select a group for analysis:\n')
    else:
        sys.stdout.write(message+'\n')

    ndxNum = -1
    while ndxNum==-1:
        ndxNum = input() #was raw_input() in python 2
        if ndxNum.isdigit()==False:
            sys.stdout.write('Wrong index group number selected (use an integer number)\n')
            ndxNum = -1
            continue

        ndxNum = int(ndxNum)
        if (ndxNum >= i) or (ndxNum < 0):
            sys.stdout.write('Wrong index group number selected\n')
            ndxNum = -1
    sys.stdout.write('Selected group %d\n\n' % ndxNum)

    res = []
    res = np.asarray(ndxDict[ndxNum])-1 # starting from 0
    return(res)


def main(argv):

    desc=('Go over multiple frames of a) coupled prot+lig; b) ligand in vacuo; c) apo protein.',
          'Rotate 1) apo protein onto coupled prot+lig; 2) ligand in vacuo onto coupled prot+lig.',
          'If .gro files contain velocities, they will be rotated accordingly',
          '',)

# define input/output files

    files= [
        FileOption("-protlig", "r/m",["pdb","gro"],"protlig.pdb", "prot+lig coupled system"),
        FileOption("-prot", "r/m",["pdb","gro"],"apo_protein.pdb", "protein structure (apo simulation)"),
        FileOption("-lig", "r/m",["pdb","gro"],"ligand.gro", "ligand structure (vacuum simulation)"),
        FileOption("-nprotlig", "r",["ndx"],"index.ndx", "index for -xray structure"),
        FileOption("-nprot", "r",["ndx"],"index.ndx", "index for -prot structure"),
        FileOption("-nlig", "r",["ndx"],"index.ndx", "index for -lig structure"),
        FileOption("-o", "w",["pdb","gro"],"fit.gro", "output"),
        ]

# define options

    options=[
        ]

    help_text = ('Go over multiple frames of a) coupled prot+lig; b) ligand in vacuo; c) apo protein.',
          'Rotate 1) apo protein onto coupled prot+lig; 2) ligand in vacuo onto coupled prot+lig.',
          'If .gro files contain velocities, they will be rotated accordingly',
          '',)

# pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options, fileoptions = files, program_desc = help_text, check_for_existing_files = False, version = "0.0" )

    mprotligName = cmdl['-protlig']
    mprotName = cmdl['-prot']
    mligName = cmdl['-lig']


    nprotlig_prot = select_ndx(cmdl['-nprotlig'],message='Select index group for protein from prot+lig structure:\n')
    nprotlig_lig = select_ndx(cmdl['-nprotlig'],message='Select index group for ligand from prot+lig structure:\n')
    nprot = select_ndx(cmdl['-nprot'],message='Select index group for protein from structure -prot:\n')
    nlig = select_ndx(cmdl['-nlig'],message='Select index group for ligand from structure -lig:\n')

    for mplName,mpName,mlName in zip(mprotligName,mprotName,mligName):
        #TODO: check if all files exist. If not, choose the missing file randomly

        # read
        mpl = Model(mplName,bPDBTER=True)
        mp = Model(mpName)
        ml = Model(mlName,bPDBTER=True)

        # all in nanometers
        mpl.a2nm()
        mp.a2nm()
        ml.a2nm()

####         # step1: fit apo protein onto the protein from prot+lig structure
####        (v1,v2,R) = fit( mpl, mp, nprotlig_prot, nprot )
####        # rotate velocities
####        rotate_velocities_R( mp, R )
####        mp.write('fitmp.pdb')
         # step1: fit prot+lig onto apo protein
        (v1,v2,R) = fit( mp, mpl, nprotlig_prot, nprot )
        # rotate velocities

        # step2: ligand in vacuo onto the ligand from prot+lig structure
        (v1,v2,R) = fit( mpl, ml, nprotlig_lig, nlig )
        # rotate velocities
        rotate_velocities_R( ml, R )

        # output
        mout = cp.deepcopy(mp)
        mout.unity = ml.unity
        # combine apo protein and ligand in vacuum
#        chID = mout.atoms[-1].chain_id
#        resID = mout.atoms[-1].resnr
        chID,resID = find_last_protein_atom( mout )
#        print chID,resID,ml.residues
        #mout.insert_residue(resID,ml.residues[0],chID)#,newResNum=True)
        #instead of adding a residue, replace the existing one
        molres=mout.residues[resID-1] #resID starts indexing at 1
        print(molres)
        mout.replace_residue(molres, ml.residues[0], bKeepResNum=True)
#change here        mout.atoms.append( list(map(lambda i: ml.atoms[i], nlig)) )
        mout.write(cmdl['-o'])
#        sys.exit(0)


if __name__ == '__main__':
    main( sys.argv )

