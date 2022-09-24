#!/usr/bin/env python

import argparse
import logging
import os
import re
import shutil as sh
import sys
from pmx.scripts.cli import check_unknown_cmd


# ==============================================================================
#                             HELPER CLASSES
# ==============================================================================
class NoMissingModuleFilter(logging.Filter):
    def filter(self, record):
        keywords=["module without the python package"]
        return not any(s in record.getMessage() for s in keywords)

# ==============================================================================
#                            HELPER FUNCTIONS
# ==============================================================================
def check_file_ready(fname):
    """Checks if a file was sucessfully created
    and gives an informative error if not.

    Parameters
    ----------
    fname : str
        Path to the file beeing checked

    Raises:
    ----------
    OSError:
        If fname is not found.
    """
    if(not os.path.isfile(fname)):
        raise OSError("Failed creating "+fname)

def copy_if_missing(src, trg):
    """Checks if trg file exists and copies it from src if not.

    Parameters
    ----------
    src : str
        Path to source file
    trg : str
        Path to target file

    Raises:
    ----------
    OSError:
        If trg still doesn't exist failed.
    """
    if(not os.path.isfile(trg)):
        sh.copy(src,trg)
        check_file_ready(trg)

def read_from_mdp(fname):
    """
    Parse mdp file for expected end time and time between saved frames.

    Returns:
    ----------
    end_time: float
        time of last frame (ps).
    dtframe: float
        time between saved frames (ps).
    """
    nsteps = 0
    tinit = 0.0
    dt = 0.0
    nstxout = 0
    with open(fname,"r") as f:
        lineList = f.readlines()
        for line in lineList:
            if("nsteps" in line):
                matchObj = re.match( r'nsteps\s*=\s*(\d+)', line, re.M|re.I)
                nsteps = int(matchObj.group(1))
            elif("tinit" in line):
                matchObj = re.match( r'tinit\s*=\s*([-+]?[0-9]*\.?[0-9]+)', line, re.M|re.I)
                tinit = float(matchObj.group(1))
            elif("dt" in line):
                matchObj = re.match( r'dt\s*=\s*([-+]?[0-9]*\.?[0-9]+)', line, re.M|re.I)
                dt = float(matchObj.group(1))
            elif("nstxout" in line):
                #print(fname,":\t",line)
                matchObj = re.match( r'nstxout\s*=\s*(\d+)', line, re.M|re.I)
                nstxout = int(matchObj.group(1))
    end_time=tinit+(nsteps*dt)
    dtframe=nstxout*dt
    return(end_time, dtframe)


# raw_input returns the empty string for "enter"
def confirm_defNO(msg):
    yes = {'y', 'yes', 'ye'}
    no = {'n','no', ''}
    sys.stderr.write(msg + "\t(y/N):\n")
    choice = input().lower()
    if choice in yes:
        return True
    elif choice in no:
        return False
    else:
        sys.stderr.write("Unexpected input string. Assuming no.\n")
        return False
    
    
    
def readii_util(fii):
    lig=[-1,-1,-1]
    pro=[-1,-1,-1]
    ligfirst=True;
    nang=0
    ndih=0
    means=[.0,.0,.0,.0,.0,.0]
    ks=[.0,.0,.0,.0,.0,.0]
    with open(fii,'r') as f:
        block=""
        for cnt, line in enumerate(f):
            l=line.strip()
            if(not l): #empty line
                continue
            if('['in l): #read block name
                s=l.split()
                block=s[1]
                continue
            s=l.split()
            # if(block=="bonds"):
                # if(int(s[0])<int(s[1])):
                    # ligfirst=False

            #assume lig is first, we'll flip in the end if needed
            if(block=="bonds"):
                lig[0]=int(s[0])
                pro[0]=int(s[1])
                means[0]=float(s[3])
                ks[0]=float(s[-1])

            elif(block=="angles"):
                means[1+nang]=float(s[4])
                ks[1+nang]=float(s[-1])
                nang+=1

            elif(block=="dihedrals"):
                if(ndih==0):
                    lig=[int(i) for i in s[0:3]] #reverse order
                    lig=lig[::-1]
                elif(ndih==2):
                    pro=[int(i) for i in s[1:4]]

                means[3+ndih]=float(s[5])
                ks[3+ndih]=float(s[-1])
                ndih+=1

        # #flip lig & pro if not ligfirst
        #if(not ligfirst):
            #lig,pro=pro,lig

    return(lig, pro, means, ks) #1-indexed becasue bynum takes that



def writeii_util(fii, lig, pro, means, ks):
    fp = open(fii, 'w')
    fp.write('\n [ intermolecular_interactions ]\n')
    
    # bonds
    fp.write(' [ bonds ]\n')
    fp.write('%6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[0], pro[0], 6, means[0], 0.0, means[0], ks[0]))
    fp.write('\n')
    
    # angles
    fp.write(' [ angles ]\n')
    fp.write('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[1], lig[0], pro[0], 1, means[1], 0.0, means[1], ks[1]))
    fp.write('%6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[0], pro[0], pro[1], 1, means[2], 0.0, means[2], ks[2]))
    fp.write('\n')

    # dihedrals
    fp.write(' [ dihedrals ]\n')
    fp.write('%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[2], lig[1], lig[0], pro[0], 2,
                                                              means[3], 0.0, means[3], ks[3]))
    fp.write('%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[1], lig[0], pro[0], pro[1], 2,
                                                              means[4], 0.0, means[4], ks[4]))
    fp.write('%6d %6d %6d %6d %6d %14.6f %14.6f %14.6f %14.6f\n' % ( lig[0], pro[0], pro[1], pro[2], 2,
                                                              means[5], 0.0, means[5], ks[5]))
    fp.write('\n')
    
    fp.close()
    

# ==============================================================================
#                      COMMAND LINE OPTIONS AND MAIN
# ==============================================================================
def parse_options(SGE=False):
    """Parse cmd-line options.

    Returns:
    ----------
    args: argparse.Namespace
        The processed command line arguments.
    """

    parser = argparse.ArgumentParser(description='Runs the whole workflow '\
            'for calculation of free energy of ligand binding '\
            'using a single instance of the ligand per box, '\
            'optimized Boresh-style restraints, '\
            'and the non-equilibrium method.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--toppath',
                        type=str,
                        dest='toppath',
                        help='Path to itp and structure files describing'
                            ' the protein and the ligand.',
                        default='../data')
    parser.add_argument('--basepath',
                        type=str,
                        dest='basepath',
                        help='Path where all everything will be done.',
                        default=os.getcwd())
    parser.add_argument('--mdppath',
                        type=str,
                        dest='mdppath',
                        help='Path to mdp files for'
                            ' the protein and ligand simulations.',
                        default='../data/mdp')
    parser.add_argument('--bt',
                        dest='bt',
                        type=str,
                        choices=['triclinic', 'cubic',
                                 'dodecahedron', 'octahedron'],
                        help='Box type.',
                        default='dodecahedron')
    parser.add_argument('-d',
                        dest='d',
                        type=float,
                        help='Distance (nm) between the solute and the box.',
                        default=1.5)
    parser.add_argument('-b',
                        dest='b',
                        type=float,
                        help='Time (ps) at which to start sampling frames '
                        'from equilibrium simulations.',
                        default=2256.0)
    parser.add_argument('--workers',
                        dest='workers',
                        type=int,
                        help='How many tasks can be running simultaneously. 0 -> one task per host/guest/repeat/direction/n_sampling sims combo.',
                        default=0)
    # parser.add_argument('--gmx',
    #                     dest='gmx',
    #                     type=str,
    #                     help='Call to gmx',
    #                     default="gmx")


    mdrun_help=""
    if(SGE):
        #parser.add_argument('--pe',
                            #dest='pe',
                            #type=str,
                            #help='Parellel environment to use for SGE jobs.',
                            #default="openmp_fast")
        parser.set_defaults(rem_sched=False)
        parser.add_argument('--rem_sched',  action='store_true',
                            dest='rem_sched',
                            help='If supplied, luigi will use a central scheduling '
                            'server to manage tasks in the workflow. '
                            'hostname and port will be read from luigi\'s '
                            'config file. '
                            'Otherwize, will use a local scheduler on the '
                            'current machine.')

        mdrun_help=' For best performance on a cluster make this value aliased for the optimal version of mdrun on each node.'

    parser.add_argument('--mdrun',
                            dest='mdrun',
                            type=str,
                            help='Call to mdrun.' + mdrun_help,
                            default="gmx mdrun")

    parser.add_argument('--mdrun_double',
                            dest='mdrun_double',
                            type=str,
                            help='Call to mdrun.' + mdrun_help,
                            default="gmx mdrun")

    parser.add_argument('--mdrun_opts',
                        dest='mdrun_opts',
                        type=str,
                        help='Optional arguments to mdrun. '
                        'Enclose in quotes.',
                        default="")



    args, unknown = parser.parse_known_args()
    check_unknown_cmd(unknown)

    return args
