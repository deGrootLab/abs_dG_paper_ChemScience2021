#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:31:13 2021

@author: ykhalak
"""


import argparse
import numpy as np
import os
from pmx import ndx
from pmx.model import Model
from pmx.scripts.workflows.utils import check_file_ready
from pmx.scripts.workflows.fit_ligs_multiframes_python3 import fit,rotate_velocities_R, find_last_protein_atom
from pmx.scripts.workflows.utils import read_from_mdp
from pmx.xtc import Trajectory

def is_file_non_zero(fpath):
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

#helper function for .gro completeness check
def _check_gro_finished(fn):
    """Checks if a gro file is complete.

    Parameters
    ----------
    fn: filename

    Returns
    -------
    Boolean: True for finished.
    """
    ret=False
    with open(fn, 'rb') as f:
        #lines = f.read().splitlines()
        #box_size_line = lines[-2] # empty line after this

        #Faster version based on https://openwritings.net/pg/python/python-read-last-line-file
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n': # find start of last line
            f.seek(-2, os.SEEK_CUR)
        while f.read(1) != b'\n': # find start of second to last line
            f.seek(-2, os.SEEK_CUR)
        box_size_line=f.readline().decode()

        if(len(box_size_line.split())==9): # box info line should have 9 columns; atom lines only have 8
            ret=True

    return(ret)

def gen_ndx_w_chains(struct_path, ndx_path, nchains):
    #create default
    os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
              ">> align.log 2>&1".format(
                gro=struct_path, out=ndx_path ) )
    #add chains
    n = ndx.IndexFile(ndx_path)
    m = Model(struct_path, bPDBTER=True)
    for i in range(nchains):
        ch = m.chains[i]
        grp = ndx.make_index_group(ch.atoms, "chain_" + ch.id)
        n.add_group(grp)
    n.write(ndx_path)


def gen_ndx_common_Calpha(apo_struct_path, holo_struct_path, apo_ndx_path, holo_ndx_path):
    """Adds an index group of common C-alpha atoms to ndx files of both apo and holo states.

    Args:
        apo_struct_path  (str): apo structure file.
        holo_struct_path (str): holo structure file.
        apo_ndx_path     (str): apo index file.
        holo_ndx_path    (str): holo index file.

    """

    n_a = ndx.IndexFile(apo_ndx_path)
    n_h = ndx.IndexFile(holo_ndx_path)
    m_a = Model(apo_struct_path, bPDBTER=True)
    m_h = Model(holo_struct_path, bPDBTER=True)

    ready_a = "C-alpha_common" in n_a.names
    ready_h = "C-alpha_common" in n_h.names
    if(ready_a and ready_h):
        return; #this was alredy done before, so skip

    common_a=[]
    common_h=[]

    res_arr_a=[r_a for r_a in m_a.residues if r_a.moltype=='protein']
    res_arr_h=[r_h for r_h in m_h.residues if r_h.moltype=='protein']

    for r_a in res_arr_a:
        for r_h in res_arr_h:
            if(r_a.orig_id == r_h.orig_id and r_a.resname == r_h.resname): #residue matches
                for a_a in r_a.atoms:
                    for a_h in r_h.atoms:
                        if(a_a.name == a_h.name and (a_a.id in n_a["C-alpha"].ids) and (a_h.id in n_h["C-alpha"].ids)):
                            common_a.append(a_a.id)
                            common_h.append(a_h.id)

    common_a.sort()
    common_h.sort()

    #add the index groups and write to files
    if(not ready_a):
        g_a = ndx.IndexGroup(ids=common_a, name="C-alpha_common")
        n_a.add_group(g_a)
        n_a.write(apo_ndx_path)
    else:
        if common_a != n_a["C-alpha_common"].ids:
            raise(Exception("Generated common protein index group does not match the one that already exists: [C-alpha_common] in %s".format(apo_ndx_path)))

    if(not ready_h):
        g_h = ndx.IndexGroup(ids=common_h, name="C-alpha_common")
        n_h.add_group(g_h)
        n_h.write(holo_ndx_path)
    else:
        if common_h != n_h["C-alpha_common"].ids:
            raise(Exception("Generated common protein index group does not match the one that already exists: [C-alpha_common] in %s".format(holo_ndx_path)))

################################################################################
#start of execution
if __name__== "__main__":
    parser = argparse.ArgumentParser(description='Computes dG from Zwanzig formula.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', dest='c', default="../../../init.pdb",
                        type=str, help='holo structure file used for identifying protein chains. Has to be a pdb.')
    parser.add_argument('-fH', dest='fH', default="../../../stateA/repeat0/npt0/",
                        type=str, help='holo NPT folder')
    parser.add_argument('-fA', dest='fA', default="../../../../apoP/repeat0/npt0/",
                        type=str, help='apo NPT folder')
    parser.add_argument('-fW', dest='fW', default="../../../../water/stateB/repeat0/npt0/",
                        type=str, help='decoupled water NPT folder')
    parser.add_argument('--no_reuse_existing_pbc_fixes', dest='reuse_existing_pbc_fixes', action='store_false')
    parser.add_argument('-b', dest='b', default=2256.0,
                        type=float, help='trj start time (ps)')
    parser.add_argument('--no_write_aligned_trj', dest='write_aligned_trj', action='store_false')
    parser.add_argument('--first_frame_only', dest='first_frame_only', action='store_true')
    parser.add_argument('-r', dest='r', type=int, help='repeat number')
    parser.add_argument('-m', dest='m', type=int, help='indep. equil. sim. number', default=0)


    args = parser.parse_args()


    if(args.c[-4:]!=".pdb"):
        print(f"Protein structure file has to be a pdb to identify chains, not '{args.c}'")
        exit(1);

    #find number of protein chains
    m_init = Model(args.c, bPDBTER=True)
    n_prot_chains=0 #count number of protein chains. Assume that they are all the ones before the ligand
    for c in m_init.chains:
        if("-MOL-" in c.get_sequence()):
            break;
        else:
            n_prot_chains+=1
    if(n_prot_chains<1):
        raise(RuntimeError("There should be at least one protein chain before the ligand! Check if the ligand has the same chain id as the protein!"))

    if(os.path.isfile("align.log")): #clean old partial log if present
        os.unlink("align.log")


    #make the ndxs for the chains
    if(not is_file_non_zero("PL_w_chains.ndx")):
        #gen_ndx_w_chains(args.fH+"/confout.gro", "PL_w_chains.ndx", n_prot_chains) #P+L
        gen_ndx_w_chains(args.fH+f"/../../../ions{args.r}_{args.m}.pdb", "PL_w_chains.ndx", n_prot_chains) #P+L
        check_file_ready("PL_w_chains.ndx")
    if(not is_file_non_zero("P_w_chains.ndx")):
        #gen_ndx_w_chains(args.fA+"/confout.gro", "P_w_chains.ndx", n_prot_chains) #P
        gen_ndx_w_chains(args.fA+f"/../../ions{args.r}_{args.m}.pdb", "P_w_chains.ndx", n_prot_chains) #P
        check_file_ready("P_w_chains.ndx")

    #Update the PL and P ndx files with a group containing common C-alpha atoms.
    #Used for alignment of structures with missing residues.
    gen_ndx_common_Calpha(
        args.fA+"/confout.gro",
        args.fH+"/confout.gro",
        "P_w_chains.ndx", "PL_w_chains.ndx")

    #LW index
    if(not is_file_non_zero("LW.ndx")):
        os.system("echo 'q' | gmx make_ndx -f {gro} -o {out} "
                  ">> align.log 2>&1".format(
                  gro=args.fW+"/confout.gro",
                  out="LW.ndx" ) )
        check_file_ready("LW.ndx")

    #Cut the begining off of trjs and center them
    names = ["A","B","C"]
    sources = [args.fH, args.fA, args.fW]
    ndxs = ["-n PL_w_chains.ndx", "-n P_w_chains.ndx", ""]
    for i in range(3):
        src=sources[i]
        if(args.reuse_existing_pbc_fixes and os.path.exists("trj_{}.trr".format(names[i]))):
            continue;
        elif(args.reuse_existing_pbc_fixes):
            print("Trying to reuse trj_{}.trr but it does not exist".format(names[i]))


        #wrap mol centers
        os.system("echo System | gmx trjconv -s {tpr} -f {trj} -o {out} "
              "-b {b} -ur compact -pbc mol"
              ">> align.log 2>&1".format(
                  tpr=src+"tpr.tpr", trj=src+"traj.trr",
                  out="trj_{}_temp_cut_pbc.trr".format(names[i]),
                  b=args.b) )

        #center on chain_A and rewrap mol centers
        c = "chain_A"
        if(i==2): #LW
            c = "Other"
        os.system("echo {c} System | gmx trjconv -s {tpr} -f {trj} -o {out} "
              "{ndx} -ur compact -center -pbc mol"
              ">> align.log 2>&1".format(
                  tpr=src+"tpr.tpr", ndx=ndxs[i], c=c,
                  trj="trj_{}_temp_cut_pbc.trr".format(names[i]),
                  out="trj_{}.trr".format(names[i]) ) )

        check_file_ready("trj_{}.trr".format(names[i]))

        #clean temp
        os.unlink("trj_{}_temp_cut_pbc.trr".format(names[i]))

    #make the C state
    m_A = Model(args.fH+"/confout.gro")
    m_B = Model(args.fA+"/confout.gro") #apoP
    m_C = Model(args.fW+"/confout.gro") #vacL
    m_A.a2nm()
    m_B.a2nm()
    m_C.a2nm()


    trj_A = Trajectory("trj_A.trr") #P+L
    trj_B = Trajectory("trj_B.trr") #apoP
    trj_C = Trajectory("trj_C.trr") #vacL



    ndx_file_A = ndx.IndexFile("PL_w_chains.ndx", verbose=False)
    ndx_file_B = ndx.IndexFile("P_w_chains.ndx", verbose=False)
    ndx_file_C = ndx.IndexFile("LW.ndx", verbose=False)
    pA_ndx = np.asarray(ndx_file_A["C-alpha_common"].ids)-1 # as in Vytas' alignment script
    pB_ndx = np.asarray(ndx_file_B["C-alpha_common"].ids)-1 # as in Vytas' alignment script
    linA_ndx = np.asarray(ndx_file_A["MOL"].ids)-1
    l_ndx = np.asarray(ndx_file_C["MOL"].ids)-1


    #find chain and resID of the last residue of the protein
    mol_first_atom = m_A.atoms[linA_ndx[0]]
    resID = mol_first_atom.resnr #internal numbering (residue.id, not residue.orig_id)
    A_mol_res_index=-1;
    A_prot_end_res_index=-1;
    for i,r in enumerate(m_A.residues):
        if(r.moltype=='protein'):
            A_prot_end_res_index = m_A.residues.index(r)
        if(r.id==resID):
            A_mol_res_index = m_A.residues.index(r)
            break;
    if(A_mol_res_index<0):
        raise("Could not find residue with resID %d in protein+ligand."%(resID))
    if(A_prot_end_res_index<0):
        raise("Could not find the last protein residue in protein+ligand.")

    mol_res_index_shift=A_mol_res_index-A_prot_end_res_index

    B_prot_end_res_index=-1;
    for i,r in enumerate(m_B.residues):
        if(r.moltype=='protein'):
            B_prot_end_res_index = m_B.residues.index(r)
    if(B_prot_end_res_index<0):
        raise("Could not find the last protein residue in ApoP.")
    B_mol_res_index=B_prot_end_res_index+mol_res_index_shift #this is where the ligand will go



    num_aligned_atoms = len(m_B.atoms) + l_ndx.shape[0]

    if(args.write_aligned_trj):
        trj_out = Trajectory("aligned.trr", mode='Out',
                             atomNum = num_aligned_atoms) #aligned output


    #Frames are not acessible individually, just in sequence.
    #pmx.xtc.Trajectory is based on __iter__, so we need a custom
    #"for" loop to simultaneously go through both trajectories.
    #Based on https://www.programiz.com/python-programming/iterator
    iter_A = iter(trj_A)
    iter_B = iter(trj_B)
    iter_C = iter(trj_C)
    fridx=0
    while True:
        try:
            frame_A = next(iter_A)
            frame_B = next(iter_B)
            frame_C = next(iter_C)
        except StopIteration:
            break

        #don't want frames while still equilibrating
        if(frame_A.time<args.b):
            continue


        frame_A.update(m_A)
        #m_b needs to be reloaded to have correct # of atoms next iteration
        m_B = Model(args.fA+"/confout.gro") #apoP
        m_B.a2nm()
        frame_B.update(m_B, uv=True)
        frame_C.update(m_C, uv=True)

        # step1: fit prot from prot+lig onto apo protein
        (v1,v2,R) = fit( m_B, m_A, pB_ndx, pA_ndx )
        # rotate velocities
        # not needed. We aren't saving m_A

        # step2: ligand onto the ligand from prot+lig structure
        (v1,v2,R) = fit( m_A, m_C, linA_ndx, l_ndx )
        # rotate velocities
        rotate_velocities_R( m_C, R )

        #insert vac ligand into B
        #do the insertion explicitly without relying on chains
        mol = m_C.residues[0]
        mol.model = m_B
        m_B.residues.insert(B_mol_res_index, mol)
        #don't add to chain, it doesn't matter
        m_B.al_from_resl()
        m_B.renumber_atoms()
        m_B.al_from_resl()

        # output
        if(fridx==0 or not args.first_frame_only): # don't output frames after first?
            frame_fn="frame%d.gro"%fridx
            if(not os.path.isfile(frame_fn) or not _check_gro_finished(frame_fn)):
                m_B.write(frame_fn)

        x = np.zeros(len(m_B.atoms)*3)
        v = np.zeros(len(m_B.atoms)*3)
        for i, atom in enumerate(m_B.atoms):
            x[i*3:(i+1)*3]=atom.x
            v[i*3:(i+1)*3]=atom.v

        if(args.write_aligned_trj):
            trj_out.write_xtc_frame(step=frame_B.step, time=frame_B.time,
                                    lam=1.0, box=frame_B.box, x=x, v=v,
                                    units=m_B.unity, bTrr=True )
        else:
            if(args.first_frame_only): # stop going through the trajectories after first frame?
                break;

        fridx+=1

    if(args.write_aligned_trj):
        trj_out.close()
    os.system("rm -f \\#*")
