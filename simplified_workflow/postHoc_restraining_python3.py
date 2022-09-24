import sys, os
import copy as cp
from pmx import *
from pmx import geometry
from pmx import ndx
from pmx.options import FileOption, Option, Commandline
from pmx.model import Model
import numpy as np
from scipy.stats import kstest

def write_dg( dg, outfile ):
    fp = open(outfile,'w')
    dgkcal = np.round(dg/4.184,2)
    sys.stdout.write('\nRestraint contribution to free energy (according to B-K): %3.4f kJ/mol\n' % dg )
    sys.stdout.write('Restraint contribution to free energy (according to B-K): %3.4f kcal/mol\n' % dgkcal )
    sys.stdout.write('\n')
    fp.write('Restraint contribution to free energy (according to B-K): %3.4f kJ/mol\n' % dg )
    fp.write('Restraint contribution to free energy (according to B-K): %3.4f kcal/mol\n' % dgkcal )
    fp.write('\n')
    fp.close()

def write_ii_atoms( ii,outiiFile ):
    fp = open(outiiFile, 'w')
    fp.write('\n [ intermolecular_interactions ]\n')
    # bonds
    if 'bonds' in ii.keys():
       fp.write(' [ bonds ]\n')
       for b in ii['bonds']:
           fp.write('%6d %6d %6d' % ( b[0].id, b[1].id, b[2] ))
           if len(b)>3:
               for x in b[3]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    # angles
    if 'angles' in ii.keys():
       fp.write(' [ angles ]\n')
       for ang in ii['angles']:
           fp.write('%6d %6d %6d %6d' % ( ang[0].id, ang[1].id, ang[2].id, ang[3] ))
           if len(ang)>4:
               for x in ang[4]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    # dihedrals
    if 'dihedrals' in ii.keys():
       fp.write(' [ dihedrals ]\n')
       for dih in ii['dihedrals']:
           fp.write('%6d %6d %6d %6d %6d' % ( dih[0].id, dih[1].id, dih[2].id, dih[3].id, dih[4] ))
           if len(dih)>5:
               for x in dih[5]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    fp.close()

def write_ii( ii,outiiFile ):
    fp = open(outiiFile, 'w')
    fp.write('\n [ intermolecular_interactions ]\n')
    # bonds
    if 'bonds' in ii.keys():
       fp.write(' [ bonds ]\n')
       for b in ii['bonds']:
           fp.write('%6d %6d %6d' % ( b[0], b[1], b[2] ))
           if len(b)>3:
               for x in b[3]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    # angles
    if 'angles' in ii.keys():
       fp.write(' [ angles ]\n')
       for ang in ii['angles']:
           fp.write('%6d %6d %6d %6d' % ( ang[0], ang[1], ang[2], ang[3] ))
           if len(ang)>4:
               for x in ang[4]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    # dihedrals
    if 'dihedrals' in ii.keys():
       fp.write(' [ dihedrals ]\n')
       for dih in ii['dihedrals']:
           fp.write('%6d %6d %6d %6d %6d' % ( dih[0], dih[1], dih[2], dih[3], dih[4] ))
           if len(dih)>5:
               for x in dih[5]:
                   fp.write(' %14.6f' % x)
           fp.write('\n')
       fp.write('\n')
    fp.close()

def calc_dih_atoms( m, ind1, ind2, ind3, ind4 ):
    dih = m.atoms[ind1].dihedral(m.atoms[ind2],m.atoms[ind3],m.atoms[ind4],degree=True)
    return(dih)
def calc_dih( m, ind1, ind2, ind3, ind4, bDegree=True ):
    a = m[ind1]
    b = m[ind2]
    c = m[ind3]
    d = m[ind4]
    v_ba = subtract_vecs( a,b )
    v_cb = subtract_vecs( b,c )
    v_dc = subtract_vecs( c,d )
    c1 = vector_prod(v_ba,v_cb)
    c2 = vector_prod(v_cb,v_dc)
    cosdih = np.dot(c1,c2)/np.sqrt(np.dot(c1,c1)*np.dot(c2,c2))
    dih = np.arccos(cosdih)
    direction = np.dot(v_ba,c2)
    sign = 1.0
    if direction>0.0:
        sign = -1.0
    dih = sign*dih
    if bDegree==True:
        return( dih*180.0/np.pi )
    else:
        return(dih)

def calc_dist_atoms( m, ligID, protID ):
    d = m.atoms[ligID] - m.atoms[protID]
    return(d)
def calc_dist( m, ligID, protID ):
    d = np.sqrt( np.power(m[ligID][0]-m[protID][0],2) + \
                 np.power(m[ligID][1]-m[protID][1],2) + \
                 np.power(m[ligID][2]-m[protID][2],2) )
    return(d)

def calc_angle_atoms( m, ind1, ind2, ind3 ):
    ang = m.atoms[ind2].angle(m.atoms[ind1],m.atoms[ind3],degree=True)
    return(ang)
def calc_angle( m, ind1, ind2, ind3, bDegree=True ):
    a = m[ind1]
    b = m[ind2]
    c = m[ind3]
    v1 = subtract_vecs( b,a )
    v2 = subtract_vecs( b,c )
    angle = np.arccos( np.dot(v1,v2)/np.sqrt(np.dot(v1,v1)*np.dot(v2,v2)) )
    if bDegree==True:
        return( angle*180.0/np.pi )
    else:
        return( angle )

#  rvec_sub(b,a,ba);
#  rvec_sub(b,c,bc);
#  return acos (cos_angle(ba,bc) );

def subtract_vecs( a,b ):
    x = a[0] - b[0]
    y = a[1] - b[1]
    z = a[2] - b[2]
    return(np.array([x,y,z]))

def vector_prod( a,b ):
    x = a[1]*b[2] - a[2]*b[1]
    y = a[2]*b[0] - a[0]*b[2]
    z = a[0]*b[1] - a[1]*b[0]
    return(np.array([x,y,z]))


def check_angle( allModels, RT, ind1, ind2, ind3, alphaLevel ):
    angle = []
    for key in allModels.keys():
#    for m in allModels:
        m = allModels[key]
        angle.append( calc_angle( m, ind1, ind2, ind3  ) )
    angle = np.array(angle)
    m = np.mean(angle)
    s = np.sqrt( np.var(angle) )
    k = RT/(np.power(s/180.0*np.pi,2.0))
#    print m,s,k
#    print 0.5*k*np.power( (m-0.0)/180.0*np.pi,2 )
#    print 0.5*k*np.power( (m-180.0)/180.0*np.pi,2 )

    if 0.5*k*np.power( (m-0.0)/180.0*np.pi,2 ) / RT <5.0:
        return(False)
    if 0.5*k*np.power( (m-180.0)/180.0*np.pi,2 ) / RT <5.0:
        return(False)
    if kstest( (angle-m)/s, 'norm')[1] < alphaLevel:
        return(False)

    return(True)

def check_dist(allModels,ind1,ind2,alphaLevel):
    dist = []
    for key in allModels.keys():
#    for m in allModels:
        m = allModels[key]
        dist.append( calc_dist(m,ind1,ind2) )
    m = np.mean(dist)
    s = np.sqrt( np.var(dist) )

    pval=kstest( (dist-m)/s, 'norm')[1]
    if pval < alphaLevel:
        sys.stdout.write('\t\check_dist fail:{}-{}={}\n'.format(ind1, ind2, pval))
        return(False)

#    if kstest( (dist-m)/s, 'norm')[1] < alphaLevel:
#        return(False)
    return(True)

def check_dih(allModels,ind1,ind2,ind3,ind4,alphaLevel):
    dih = []
    for key in allModels.keys():
#    for m in allModels:
        m = allModels[key]
        dih.append( calc_dih(m,ind1,ind2,ind3,ind4) )
    m = np.mean(dih)
    s = np.sqrt( np.var(dih) )

    if kstest( (dih-m)/s, 'norm')[1] < alphaLevel:
        return(False)
    return(True)

def identify_atom_pairs( distMatVar, n, N, allModels, RT, ligAtomDict, protAtomDict, alphaLevel ):
    ligList = []
    protList = []
#    forbiddenLigList = [] # list of atoms in case several iterations needed
#    forbiddenProtList = [] # list of atoms in case several iterations needed
    forbiddenLigProtList = [] # list of ligand_atom pairs already identified in an unsuccessful iteration
    minVarList = []
    backupDistMatVar = cp.deepcopy(distMatVar)

    sys.stdout.write('\nStarting anchor search with alphaLevel=%f\n'%alphaLevel)
    reset_counter=0;

    found = 0
    counter = 0
    while found<3:
        if counter==np.shape(distMatVar)[0]:
#            sys.stdout.write('\nCould not identify atom pairs. Exiting...\n')
            sys.stdout.write('\nCould not identify atom pairs. Trying again... alphaLevel=%f\treset_counter=%f\n'%(alphaLevel,reset_counter))
            sys.stdout.flush()
            ligList = []
            protList = []
            counter = 0
            found = 0
            distMatVar = cp.deepcopy(backupDistMatVar)
            reset_counter+=1
            #if(reset_counter>2):
            #    alphaLevel*=0.8;
            if(reset_counter>10):
                sys.exit(1)
        counter+=1
        ind = np.argmin(distMatVar)
        ligInd = np.divmod( ind, N )[0]
        protInd = np.divmod( ind, N)[1]
        ligProtInd = str(ligInd)+'_'+str(protInd)
        if (ligInd not in ligList) and (protInd not in protList) and (found!=0 or ligProtInd not in forbiddenLigProtList):
#           and (ligInd not in forbiddenLigList) and (protInd not in forbiddenProtList):
            if found==0:
                # check dist
                if check_dist(allModels,ligAtomDict[ligInd],protAtomDict[protInd],alphaLevel)==False:
                    distMatVar[ind] = 99999.99
                    #sys.stdout.write('\t{} failed found=0\n'.format(ligProtInd))
                    #sys.stdout.flush()
                    continue
                else:
                     forbiddenLigProtList.append(ligProtInd)
#                    forbiddenLigList.append(ligInd)
#                    forbiddenProtList.append(protInd)
            if found==1:
                # check angle1 and dih2
                if check_angle(allModels,RT,ligAtomDict[ligInd],ligAtomDict[ligList[0]],protAtomDict[protList[0]],alphaLevel)==False or check_dih(allModels,ligAtomDict[ligInd],ligAtomDict[ligList[0]],protAtomDict[protList[0]],protAtomDict[protInd],alphaLevel)==False:
                    distMatVar[ind] = 99999.99
                    #sys.stdout.write('\t{} failed found=1\n'.format(ligProtInd))
                    continue
#                else:
#                    forbiddenLigList.append(ligInd)
#                    forbiddenProtList.append(protInd)
            if found==2:
                # check angle2 and dih1 and dih3
                if check_angle(allModels,RT,ligAtomDict[ligList[0]],protAtomDict[protList[0]],protInd,alphaLevel)==False or check_dih(allModels,ligAtomDict[ligInd],ligAtomDict[ligList[1]],ligAtomDict[ligList[0]],protAtomDict[protList[0]],alphaLevel)==False or check_dih(allModels,ligAtomDict[ligList[0]],protAtomDict[protList[0]],protAtomDict[protList[1]],protAtomDict[protInd],alphaLevel)==False:
                    distMatVar[ind] = 99999.99
                    #sys.stdout.write('\t{} failed found=2\n'.format(ligProtInd))
                    continue

            ligList.append(ligInd)
            protList.append(protInd)
            minVarList.append( distMatVar[ind] )
            found+=1
            #sys.stdout.write('\tfound={}: {}\t{}\n'.format(found,ligList, protList))
            #sys.stdout.flush()
        distMatVar[ind] = 99999.99

#    print ligList
#    print protList
#    print minVarList
    return(ligList,protList)

def calc_distances( arrLigx, arrLigy, arrLigz, arrProtx, arrProty, arrProtz ):
    t = np.shape(arrLigx)[0]
    n = np.shape(arrLigx)[1]
    N = np.shape(arrProtx)[1]

    out = np.zeros((t,n*N))
#    distances from one ligand atom to each protein atom
    indLig = np.repeat( range(0,np.shape(arrLigx)[1]), np.shape(arrProtx)[1] )
    indProt = np.tile( range(0,np.shape(arrProtx)[1]), np.shape(arrLigx)[1] )

    for row in range(0,t):
        distx = np.power(arrLigx[row,indLig] - arrProtx[row,indProt],2)
        disty = np.power(arrLigy[row,indLig] - arrProty[row,indProt],2)
        distz = np.power(arrLigz[row,indLig] - arrProtz[row,indProt],2)
        out[row,] = distx+disty+distz

    return(out)

def build_atom_dict( ngroup ):
    out = {}
    i = 0
    for n in ngroup:
        out[i] = n
        i+=1
    return(out)

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
        ndxNum = input()
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

    desc=('In a generated ensemble of structures with protein and ligand translated+rotated',
          'into the binding site, identify optimal restraints',
          '',)

# define input/output files

    files= [
        FileOption("-f", "r/m",["pdb","gro"],"protlig.pdb", "prot+lig system"),
        FileOption("-n", "r",["ndx"],"index.ndx", "index file"),
        FileOption("-oii", "w",["itp"],"ii.itp", "output restraints"),
        FileOption("-odg", "w",["dat"],"dg.dat", "output dg"),
        ]

# define options

    options=[
       Option( "-T", "float", 298.0, "temperature"),
       Option( "-alpha", "float", 0.05, "alpha level"),
        ]

    help_text = ('In a generated ensemble of structures with protein and ligand translated+rotated',
          'into the binding site, identify optimal restraints',
          '',)

# pass options, files and the command line to pymacs

    cmdl = Commandline( argv, options = options, fileoptions = files, program_desc = help_text, check_for_existing_files = False, version = "0.0" )

    mNames = cmdl['-f']
    alphaLevel = cmdl['-alpha']

    nprot = select_ndx(cmdl['-n'],message='Select index group for protein:\n')
    nlig = select_ndx(cmdl['-n'],message='Select index group for ligand:\n')

    # build dictionaries for protein and ligand
    ligAtomDict = build_atom_dict( nlig )  # ligAtomDict[0...nlig] = ligAtomID
    protAtomDict = build_atom_dict( nprot )  # protAtomDict[0...nprot] = protAtomID

    arrLigx = []
    arrLigy = []
    arrLigz = []
    arrProtx = []
    arrProty = []
    arrProtz = []
    allModels = {}
    counter = 0
    for mName in mNames:
        # read
        sys.stdout.write('\r Reading: %s\n' % mName)
        sys.stdout.flush()
        m = Model(mName,bPDBTER=True)
#        allModels.append(m)
        #TODO: save only the protein and ligand atoms
        allModels[counter] = list(map(lambda i: m.atoms[i].x, range(0,np.shape(m.atoms)[0])))
#        nprotlig = np.hstack((nprot,nlig))
#        mProtLig = Model()
#        mProtLig.atoms = map(lambda i: m.atoms[i], nprotlig)
#        foo = cp.deepcopy(mProtLig)
#        allModels.append(foo)

        # all in nanometers
#        m.a2nm()

        # extract atoms into arrays
        fooLigx = list(map(lambda i: m.atoms[i].x[0], nlig))
        fooLigy = list(map(lambda i: m.atoms[i].x[1], nlig))
        fooLigz = list(map(lambda i: m.atoms[i].x[2], nlig))
        fooProtx = list(map(lambda i: m.atoms[i].x[0], nprot))
        fooProty = list(map(lambda i: m.atoms[i].x[1], nprot))
        fooProtz = list(map(lambda i: m.atoms[i].x[2], nprot))
        if counter==0:
            arrLigx = fooLigx
            arrLigy = fooLigy
            arrLigz = fooLigz
            arrProtx = fooProtx
            arrProty = fooProty
            arrProtz = fooProtz
        else:
            arrLigx = np.vstack([arrLigx,fooLigx])
            arrLigy = np.vstack([arrLigy,fooLigy])
            arrLigz = np.vstack([arrLigz,fooLigz])
            arrProtx = np.vstack([arrProtx,fooProtx])
            arrProty = np.vstack([arrProty,fooProty])
            arrProtz = np.vstack([arrProtz,fooProtz])

        counter+=1
        
    sys.stdout.write('\r Finished reading %d frames\n' % counter)
    sys.stdout.flush()

    # calculate distances
    distMat = calc_distances( arrLigx, arrLigy, arrLigz, arrProtx, arrProty, arrProtz )

    # calculate variances in distances
    distMatVar = np.var(distMat,axis=0)
    t = np.shape(arrLigx)[0]
    n = np.shape(arrLigx)[1]
    N = np.shape(arrProtx)[1]
#    distMatVar = distMatVar.reshape((n,N)) # rows - ligand atoms, cols - protein atoms

    # identify optimal atom pairs (combine with angle calculations)
    R = 8.31445985*0.001  # Gas constant in kJ/mol/K
    RT = cmdl['-T']*R
    sys.stdout.write('Identifying atom pairs...\n')
    sys.stdout.flush()
    ligList,protList = identify_atom_pairs( distMatVar, n, N, allModels, RT, ligAtomDict, protAtomDict, alphaLevel )
    print("ligList: {}\nprotList: {}\n".format(ligList,protList))
    sys.stdout.flush()

    # calculate restraints
    sys.stdout.write('Calculating restraints...\n')
    dist = []
    angle1 = []
    angle2 = []
    dih1 = []
    dih2 = []
    dih3 = []
    for key in allModels.keys():
        m = allModels[key]
        # dist
        dist.append( calc_dist( m, ligAtomDict[ligList[0]], protAtomDict[protList[0]] ) )
        # angles
        angle1.append( calc_angle( m, ligAtomDict[ligList[1]], ligAtomDict[ligList[0]], protAtomDict[protList[0]] ) )
        angle2.append( calc_angle( m, ligAtomDict[ligList[0]], protAtomDict[protList[0]], protAtomDict[protList[1]] ) )
        # dihedrals
        dih1.append( calc_dih( m, ligAtomDict[ligList[2]], ligAtomDict[ligList[1]], ligAtomDict[ligList[0]], protAtomDict[protList[0]] ) )
        dih2.append( calc_dih( m, ligAtomDict[ligList[1]], ligAtomDict[ligList[0]], protAtomDict[protList[0]], protAtomDict[protList[1]] ) )
        dih3.append( calc_dih( m, ligAtomDict[ligList[0]], protAtomDict[protList[0]], protAtomDict[protList[1]], protAtomDict[protList[2]] ) )
    dist = np.array(dist)
    angle1 = np.array(angle1)
    angle2 = np.array(angle2)
    dih1 = np.array(dih1)
    dih2 = np.array(dih2)
    dih3 = np.array(dih3)

    # perform KS test
    Ddist,pdist = kstest( (dist-np.mean(dist))/np.sqrt(np.var(dist)), 'norm')
    Dangle1,pangle1 = kstest( (angle1-np.mean(angle1))/np.sqrt(np.var(angle1)), 'norm')
    Dangle2,pangle2 = kstest( (angle2-np.mean(angle2))/np.sqrt(np.var(angle2)), 'norm')
    Ddih1,pdih1 = kstest( (dih1-np.mean(dih1))/np.sqrt(np.var(dih1)), 'norm')
    Ddih2,pdih2 = kstest( (dih2-np.mean(dih2))/np.sqrt(np.var(dih2)), 'norm')
    Ddih3,pdih3 = kstest( (dih3-np.mean(dih3))/np.sqrt(np.var(dih3)), 'norm')
#    print pdist,pangle1,pangle2,pdih1,pdih2,pdih3

    # fit force
    # distance
    mdist = np.mean(dist)
    kdist = RT/(np.var(dist))
    # angles
    mangle1 = np.mean(angle1)
    kangle1 = RT/(np.var(angle1/180.0*np.pi))
    mangle2 = np.mean(angle2)
    kangle2 = RT/(np.var(angle2/180.0*np.pi))
    # dihedrals
    mdih1 = np.mean(dih1)
    kdih1 = RT/(np.var(dih1/180.0*np.pi))
    mdih2 = np.mean(dih2)
    kdih2 = RT/(np.var(dih2/180.0*np.pi))
    mdih3 = np.mean(dih3)
    kdih3 = RT/(np.var(dih3/180.0*np.pi))

    # output restraints
    ii = {}
#    lig1 = allModels[0].atoms[ligAtomDict[ligList[0]]]
#    lig2 = allModels[0].atoms[ligAtomDict[ligList[1]]]
#    lig3 = allModels[0].atoms[ligAtomDict[ligList[2]]]
#    prot1 = allModels[0].atoms[protAtomDict[protList[0]]]
#    prot2 = allModels[0].atoms[protAtomDict[protList[1]]]
#    prot3 = allModels[0].atoms[protAtomDict[protList[2]]]
    lig1 = ligAtomDict[ligList[0]]+1
    lig2 = ligAtomDict[ligList[1]]+1
    lig3 = ligAtomDict[ligList[2]]+1
    prot1 = protAtomDict[protList[0]]+1
    prot2 = protAtomDict[protList[1]]+1
    prot3 = protAtomDict[protList[2]]+1

    bond1 = [ lig1, prot1, 6, [mdist, 0.0, mdist, kdist] ]
    ii['bonds'] = [bond1]
    angle1 = [ lig2, lig1, prot1, 1, [mangle1, 0.0, mangle1, kangle1] ]
    angle2 = [ lig1, prot1, prot2, 1, [mangle2, 0.0, mangle2, kangle2] ]
    ii['angles'] = [ angle1, angle2 ]
    dihedral1 = [ lig3, lig2, lig1, prot1, 2, [mdih1, 0.0, mdih1, kdih1] ]
    dihedral2 = [ lig2, lig1, prot1, prot2, 2, [mdih2, 0.0, mdih2, kdih2] ]
    dihedral3 = [ lig1, prot1, prot2, prot3, 2, [mdih3, 0.0, mdih3, kdih3] ]
    ii['dihedrals'] = [ dihedral1, dihedral2, dihedral3 ]
    write_ii( ii, cmdl['-oii'] )

    # calculate dG contribution
    V0 = 1.66            # standard volume in nm^3
    dgPrefactor = ( 8.0*np.power(np.pi,2.0)*V0/(np.power(mdist,2.0)*np.sin(mangle1*np.pi/180.0)*np.sin(mangle2*np.pi/180.0)) )
    dgForceConstants = np.sqrt(kdist*kangle1*kangle2*kdih1*kdih2*kdih3)/np.power(2.0*np.pi*RT,3.0)
    dg = -RT * np.log(dgPrefactor*dgForceConstants)
    dg = np.round(dg,2)

    # output dG contribution
    write_dg( dg, cmdl['-odg'] )

if __name__ == '__main__':
    main( sys.argv )

