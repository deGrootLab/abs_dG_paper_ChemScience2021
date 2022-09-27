#!/usr/bin/env python

import numpy as np
import os
import os.path
from pmx import ndx
from pmx.model import Model
from utils import read_from_mdp, readii_util
from pmx.xtc import Trajectory
from pmx.geometry import Rotation
from postHoc_restraining_python3 import calc_dist, calc_angle, calc_dih, vector_prod, subtract_vecs

# ==============================================================================
#                               Helper Functions
# ==============================================================================
def readii(fii):
    #Keep this function around becasue other external scripts (eg energetic decorelation) import it
    lig, pro, means, ks = readii_util(fii)
    return(np.array(lig), np.array(pro), np.array(means), np.array(ks)) #1-indexed becasue bynum takes that

def update_anchors(model, idx_lig, idx_pro):
        L=[model.atoms[idx_lig[a]-1].x for a in range(3)]
        L.reverse()
        P=[model.atoms[idx_pro[a]-1].x for a in range(3)]
        return(L+P) #ligand first: [l2 l1 l0 p0 p2 p3]

def print_cur_rot(anchors, global_idcs, goal):
    order=[ #type, mean&sigma id, anchor indeces
            ["dih", 5, [5,4,3,2], "dih_C"], #dih_C [P2, P1, P0, L2]
            ["ang", 2,   [4,3,2], "ang_B"], #ang_B [P1, P0, P2]
            ["dist",0,     [3,2], "dist"],  #dist [P0, L2]
            ["dih", 4, [4,3,2,1], "dih_B"], #dih_B [P1, P0, L2, L1]
            ["ang", 1,   [3,2,1], "ang_A"], #ang_A [P0, L2, L1]
            ["dih", 3, [3,2,1,0], "dih_A"]  #dih_A [P0, L2, L1, L0]
          ]
    for op in order:
        cur_val=0
        indeces=op[2]
        if(op[0] == "dih"):
            cur_val = calc_dih(anchors,*indeces,False)
        elif(op[0] == "ang"):
            cur_val = calc_angle(anchors,*indeces,False)
        elif(op[0] == "dist"):
            cur_disp = subtract_vecs(anchors[indeces[1]],
                                     anchors[indeces[0]])
            cur_val = np.linalg.norm(cur_disp)

        print("\t{}:\t{:>7.4}  {:>7.4}\t{}".format(op[3], cur_val, goal[op[1]], [global_idcs[a] for a in indeces]))

def rotate_and_translate_Lig_to(g, lig, model, idx_lig, idx_pro, order=None, debug=False):
    goal=g

    if(not order):
        order=[ #type, mean&sigma id, anchor indeces
                ["dih", 5, [5,4,3,2], "dih_C"], #dih_C [P2, P1, P0, L2]
                ["ang", 2,   [4,3,2], "ang_B"], #ang_B [P1, P0, P2]
                ["dist",0,     [3,2], "dist"],  #dist [P0, L2]
                ["dih", 4, [4,3,2,1], "dih_B"], #dih_B [P1, P0, L2, L1]
                ["ang", 1,   [3,2,1], "ang_A"], #ang_A [P0, L2, L1]
                ["dih", 3, [3,2,1,0], "dih_A"]  #dih_A [P0, L2, L1, L0]
              ]

    opid=0
    for op in order:
        indeces=op[2]
        anchors=update_anchors(model, idx_lig, idx_pro)
        if(op[0] == "dih"):
            cur_val = calc_dih(anchors,*indeces,False)
            delta = goal[op[1]] - cur_val
            rot = Rotation(anchors[indeces[1]],anchors[indeces[2]])
            for num, a in enumerate(lig): #rot needs to be around 0,0,0, same as start of rot axis
                a.x=rot.apply(a.x, delta)
                a.v=np.array(rot.apply(np.array(a.v)+np.array(anchors[indeces[2]]), delta))-np.array(anchors[indeces[2]])
            if(debug):
                print("\tafter op #", opid, op[3])
                print("\t",lig[0].x, '\t', lig[0].v, "len vel:", np.linalg.norm(np.array(lig[0].v)))
        elif(op[0] == "ang"):
            cur_val = calc_angle(anchors,*indeces,False)
            delta = goal[op[1]] - cur_val
            normal = vector_prod(
                subtract_vecs(anchors[indeces[2]], anchors[indeces[1]]),
                subtract_vecs(anchors[indeces[1]], anchors[indeces[0]])
                )
            endp = np.array(anchors[indeces[1]])+np.array(normal)
            rot = Rotation(anchors[indeces[1]], endp)
            #for a in lig: #rot needs to be around 0,0,0, same as start of rot axis
            for num, a in enumerate(lig):
                a.x=rot.apply(a.x, delta)
                #velocity needs to be roatated around an axis starting at [0,0,0], these shidts achieve that
                a.v=np.array(rot.apply(np.array(a.v)+endp, delta))-endp
            if(debug):
                print("after op #", opid, op[3])
                print("\t",lig[0].x, '\t', lig[0].v, "len vel:", np.linalg.norm(np.array(lig[0].v)))
        elif(op[0] == "dist"):
            cur_disp = subtract_vecs(anchors[indeces[1]],
                                     anchors[indeces[0]])
            cur_val = np.linalg.norm(cur_disp)
            final_disp = cur_disp * goal[op[1]]/cur_val
            delta = final_disp - cur_disp
            for a in lig:
                for ax in range(3):
                    a.x[ax] = a.x[ax] + delta[ax]
                    #no change for velocity here

            if(debug):
                print("after op #", opid, op[3])
                print("\t",lig[0].x, '\t', lig[0].v, "len vel:", np.linalg.norm(np.array(lig[0].v)))
        else:
            raise(Exception("Unrecognized restraint type %s"%op[0]))

        opid+=1


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_decorelate_alg(SGETunedJobTask):
    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')
    m = luigi.IntParameter(description='Sampling sim number')
    sTI = luigi.Parameter(description='Coupling state for TI')

    folder_path = luigi.Parameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=False,
                 visibility=ParameterVisibility.HIDDEN,
                 description='Dict of study stettings '
                 'used to propagate settings to dependencies')

    restr_scheme = luigi.Parameter(significant=True,
                 description='Restraint scheme to use. '
                 'Aligned, Fitted or Fixed')

    T = luigi.FloatParameter(significant=False,
                 default=298.0, #K
                 description='Temperature in Kelvin. '
                 'Used to build distribution from restraint strengths.')

    stage="morphes"

    #request 1 cores
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    #avoid Prameter not a string warnings
    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{sTI}{i}_{m}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    #debug output
    debug = luigi.BoolParameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default=False,
        description="show debug output in a log.")

    extra_packages=[np]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._setupState()

        #set variables
        self.base_path = self.study_settings['base_path']
        self.sim_path = self.folder_path+"/state%s/repeat%d/%s%d"%(
            self.sTI, self.i, self.stage, self.m)
        self.mdp = self.study_settings['mdp_path'] +\
            "/protein/eq_npt_{0}.mdp".format(
                self.study_settings['TIstates'][self.sTI])

    def _setupState(self):
        if(self.sTI != "C"):
            raise(ValueError("Aligning morphes for TI state{}. "
                     "{} should only be done on TI stateC.".format(
                         self.sTI, self.__class__.__name__)))
            exit(1);
        # no need to set self.s as it is never used
        # else:
        #     self.s="B" # TI stateC depends on NPT sim in stateB



def work(args):
    #make the C state

    trj_A = Trajectory(args.f) #P+L

    m_A = Model(args.c,bPDBTER=False)
    m_A.a2nm()

    trj_out = Trajectory(args.o, mode='Out',
                            atomNum=len(m_A.atoms)) #aligned output

    ndx_file_A = ndx.IndexFile(args.i, verbose=False)
    if(not "MOL" in ndx_file_A):
        raise(Exception(f"index file {args.i} does not contain a [ MOL ] group!"))
    l_ndx = np.asarray(ndx_file_A["MOL"].ids)-1

    #select ligand
    lig=list(map(lambda i: m_A.atoms[i], l_ndx))


    idx_lig, idx_pro, means, ks = readii(args.ii)
    conv_ks=ks
    for k in range(1,6):
        conv_ks[k]=conv_ks[k]*(np.pi/180.0)**2
    kT = 8.31445985*0.001*self.T
    sigmas = np.sqrt(kT/conv_ks)
    #means and sigmas in deg and nm

    cov_mat=np.zeros((6,6))
    for j in range(6):
        if(j>0): #not the distance
            means[j]*=np.pi/180.0 #convert to rad
            sigmas[j]*=np.pi/180.0 #convert to rad
        cov_mat[j,j]=sigmas[j]*sigmas[j]

    if(args.debug):
        print("debug: starting on trajectory: {}".format(trj_A_src))

    #output frame files
    ods=os.path.splitext(args.od)

    #Frames are not acessible individually, just in sequence.
    #pmx.xtc.Trajectory is based on __iter__, so we need a custom
    #"for" loop to simultaneously go through both trajectories.
    #Based on https://www.programiz.com/python-programming/iterator
    iter_A = iter(trj_A)
    fridx=0
    while True:
        try:
            frame_A = next(iter_A)
        except StopIteration:
            break

        od=f"{ods[0]}{fridx}{ods[1]}"
        if(not os.path.isfile(od)):
            frame_A.update(m_A, uv=True, uf=True)
            if(self.debug):
                print("debug: \tframe {}".format(fridx))

            #draw restraint coords from independent multivariate distribution
            sample = np.random.multivariate_normal(means, cov_mat) # nm & rad

            #wrap sample so that dihedrals are in (-180,180). Angles shouldn't need this.
            for k in range(3,len(sample)):
                for k in range(3,len(sample)):
                    #if(k>=3): #dihedrals
                    if(sample[k]<-np.pi):
                        sample[k]+=2*np.pi
                    elif(sample[k]>np.pi):
                        sample[k]-=2*np.pi

            #rotate ligand to satisfy drawn restraint coords
            rotate_and_translate_Lig_to(sample, lig, m_A, idx_lig, idx_pro, debug=False)

            # output
            if(args.write_aligned_trj):
                x = np.zeros(len(m_A.atoms)*3)
                v = np.zeros(len(m_A.atoms)*3)
                for i, atom in enumerate(m_A.atoms):
                    x[i*3:(i+1)*3]=atom.x
                    v[i*3:(i+1)*3]=atom.v
                # v=None
                trj_out.write_xtc_frame(step=frame_A.step, time=frame_A.time,
                                        lam=1.0, box=frame_A.box, x=x, v=v,
                                        units=m_A.unity, bTrr=True )
            m_A.write(od)
        else:
            if(args.debug):
                print("debug: \tframe {} exists".format(fridx))
        fridx+=1


    trj_out.close()
    if(args.debug):
        print("debug: done with trajectory")





################################################################################
#start of execution
if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Rotates the ligand in the aligned trajectory to satisfy the uncorrelated restraints by transtalion/rotation of the ligand.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', dest='c', default="frame0.gro",
                        type=str, help='a frame of the aligned structure')
    parser.add_argument('-f', dest='f', default="aligned.trr",
                        type=str, help='aligned trajectory')
    parser.add_argument('-i', dest='i', default="index_prot_mol.ndx",
                        type=str, help='Index file with a [ MOL ] group for the ligand to be translated/rotated')
    parser.add_argument('-ii', dest='ii', default="ii.itp",
                        type=str, help='restraints file')
    parser.add_argument('-o', dest='o', default="decorrelated.trr",
                        type=str, help='decorelated trajectory')
    parser.add_argument('-od', dest='od', default="start.gro",
                        type=str, help='individual frames of the decorelated trajectory. A frame number will be added before the file extention.')
    #parser.add_argument('--no_reuse_existing_pbc_fixes', dest='reuse_existing_pbc_fixes', action='store_false')
    #parser.add_argument('-b', dest='b', default=2256.0,
                        #type=float, help='trj start time (ps)')
    #parser.add_argument('--no_write_aligned_trj', dest='write_aligned_trj', action='store_false')
    parser.add_argument('--first_frame_only', dest='first_frame_only', action='store_true')
    parser.add_argument('--debug', dest='debug', action='store_true')


    args = parser.parse_args()

    work(args)
