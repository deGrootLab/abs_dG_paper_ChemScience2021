#!/usr/bin/env python

import glob
import luigi
import os
import sys
import pmx.scripts.analyze_dhdl
from luigi.parameter import ParameterVisibility
from pmx.model import Model
from pmx.xtc import Trajectory
#import pmx.scripts.workflows.SGE_tasks.SGETunedRunner as sge_runner
#from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask, extended_build_qsub_command, _parse_qsub_job_id #tuned for the owl cluster
from pmx.scripts.workflows.SGE_tasks.SGETunedJobTask import SGETunedJobTask
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.morphes import Task_PL_gen_morphes
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.analysis import Task_PL_analysis_aligned
from pmx.scripts.workflows.SGE_tasks.absFE.LinP.decorrelate_algortimically import readii, update_anchors
from pmx.scripts.workflows.postHoc_restraining_python3 import calc_dist, calc_angle, calc_dih
from pmx.scripts.workflows.find_anchors_and_write_ii import wrap_ang, wrap_dih
#from luigi.contrib.sge import logger

import numpy as np
import scipy as sp
import shutil
#import subprocess
from copy import deepcopy

try:
   import cPickle as pickle
except:
   import pickle

# ==============================================================================
#                              Helper Functions
# ==============================================================================
def build_aligned_distib_no_data(ii, struct, traj, expected_frames, indiv_frames=False, tpr=None):
    idx_lig, idx_pro, means, ks = readii(ii)
    """
    Build a 6D distribution of aligned restraints seen after alignment
    and return the mean and covariance matrix in 6D.
    """

    for i in range(6):
        #convert rad^-2 to deg^-2
        if(i>0):
            ks[i]/=(180/np.pi)**2
#         print ("mean={}\t\tks={}\t\t=sigma^2={}".format(means[i], ks[i], R*T/ks[i]))


    dist=[]
    angA=[]
    angB=[]
    dihA=[]
    dihB=[]
    dihC=[]


    m = Model(struct)
    m.a2nm() #set units to nm
    if(not indiv_frames):
        trj_m = Trajectory(traj) #P+L
        iter_m = iter(trj_m)
    else:
        frame_lst=glob.glob(traj)
        frame_id=0


    fr=0


    while True:
        try:
            if(not indiv_frames):
                frame_m = next(iter_m)
                frame_m.update(m)
            else:
                if(frame_id>=len(frame_lst)):
                    raise(StopIteration())
                #fix pbc
                frame_pbc_fn = frame_lst[frame_id][:-4]+"_pbc"+frame_lst[frame_id][-4:]
                if(not os.path.isfile(frame_pbc_fn)):
                    if(not tpr):
                        raise(Exception("Missing tpr for trjconv -pbc mol!"))
                    os.system("echo '0\n' | gmx trjconv -s {tpr} -f {f} -o {o} -pbc mol".format(
                               tpr=tpr, f=frame_lst[frame_id], o=frame_pbc_fn))
                m = Model(frame_pbc_fn)
                frame_id+=1

        except StopIteration:
            break


        anchors=update_anchors(m, idx_lig, idx_pro)

        dist.append(calc_dist(  anchors, 3, 2 )) #nm
        angA.append(calc_angle( anchors, 3, 2, 1)) #deg
        angB.append(calc_angle( anchors, 4, 3, 2)) #deg

        dihC.append(calc_dih(   anchors, 5, 4, 3, 2)) #deg
        dihB.append(calc_dih(   anchors, 4, 3, 2, 1)) #deg
        dihA.append(calc_dih(   anchors, 3, 2, 1, 0)) #deg

        fr+=1
    # print("went through {} frames for repeat {}".format(fr, self.i))
    if(fr != expected_frames):
        raise(Exception("too few frames in {traj}"))



    #make a coarse histogram to find the peak
    bscale=10 # deg/bin
    angA_h = np.histogram(angA, bins=int(180/bscale), range=(0,180), density=True)
    angB_h = np.histogram(angB, bins=int(180/bscale), range=(0,180), density=True)
    dihA_h = np.histogram(dihA, bins=int(360/bscale), range=(-180,180), density=True)
    dihB_h = np.histogram(dihB, bins=int(360/bscale), range=(-180,180), density=True)
    dihC_h = np.histogram(dihC, bins=int(360/bscale), range=(-180,180), density=True)


    dists=[None, angA_h, angB_h, dihA_h, dihB_h, dihC_h]
    data=[dist, angA, angB, dihA, dihB, dihC]
    for i in range(len(data)):
        data[i]=np.array(data[i], dtype=np.float64)

        #wrap values to center them on the peakdists=[dist_h, angA_h, angB_h, dihA_h, dihB_h, dihC_h]
        if(i>0):
            d=dists[i]
            dists[i] = [d[0], d[1][:-1]+(d[1][1]-d[1][0])*0.5]
            m=dists[i][1][np.argmax(dists[i][0])] #coord of max height of distribution
            if(i<3): #wrap angles
                data[i]=wrap_ang(data[i],m)
            else:    #wrap dihedrals
                data[i]=wrap_dih(data[i],m)

    #data[0]*=100.0 #nm * 10^-2 to avoid floating point issues in eigh()

    new_means=np.mean(data, axis=-1)
    nev_cov=np.cov(np.array(data), bias=True)

    #return(new_means, nev_cov, data)
    return(new_means, nev_cov)



# ==============================================================================
def prob_pdf(m, cov):
    rv=sp.stats.multivariate_normal(m, cov)
    return(rv.pdf)

def U(x, pdf, T):
    return(-0.008314463*T*np.log(pdf(x))) #kJ/mol

def MCI(f, minlim_wp, maxlim_wp, N, max_N=int(1e5)):
    #x=np.random.uniform(size=(len(minlim_wp), N))
    #dif=maxlim_wp-minlim_wp

    #x*=dif[:, np.newaxis]
    #x+=minlim_wp[:, np.newaxis]

    #x=np.transpose(x)

    #y=f(x)
    #s=np.sum(y)
    #return((np.prod(dif)/N)*s)

    dif=maxlim_wp-minlim_wp
    N_tot=N
    N_left=N
    s=0
    while N_left>0:
        N=min(max_N, N_left)
        N_left-=N

        x=np.random.uniform(size=(len(minlim_wp), N))#.astype(dtype=np.float128)
        x*=dif[:, np.newaxis]
        x+=minlim_wp[:, np.newaxis]

        x=np.transpose(x)

        y=f(x)
        s+=np.sum(y)
    return((np.prod(dif)/N_tot)*s)
# ==============================================================================
def dG_by_MCI_between_two_distribs(m, cov, m_d, cov_d, T, N=10000000, reps=100, nsigma=4.5):
    def conv_to_rad(m,cov):
        mr = deepcopy(m)
        mr[1:] = np.pi*mr[1:]/180.
        covr = deepcopy(cov)
        covr[1:,:] = np.pi*covr[1:,:]/180.
        covr[:,1:] = np.pi*covr[:,1:]/180.

        minlim=[]
        maxlim=[]
        for i in range(6):
            minlim.append(mr[i]-nsigma*np.sqrt(covr[i,i]))
            maxlim.append(mr[i]+nsigma*np.sqrt(covr[i,i]))

            #hard limits
            if(i==0): #distance
                if(minlim[i]<0.0):
                    minlim[i]=0.0
            elif(i>0 and i<3): #angles
                # limits in gromacs are theta in [0; pi], as theta is evaluated through arccos()
                # see https://manual.gromacs.org/documentation/2021.2/reference-manual/functions/bonded-interactions.html#harmonic-angle-potential
                if(minlim[i]<0.0):
                    minlim[i]=0.0
                if(maxlim[i]>np.pi):
                    maxlim[i]=np.pi

            elif(i>=3): #dihedrals
                # limits in gromacs are (phi-mu) in [-pi; pi]
                if(minlim[i]<mr[i]-np.pi):
                    minlim[i]=mr[i]-np.pi
                if(maxlim[i]>mr[i]+np.pi):
                    maxlim[i]=mr[i]+np.pi

        minlim=np.array(minlim)
        maxlim=np.array(maxlim)
        return(mr, covr, minlim, maxlim)

    m_c, cov_c, minlim_c, maxlim_c = conv_to_rad(m,cov)
    m_d, cov_d, minlim_d, maxlim_d = conv_to_rad(m_d,cov_d)

    pdf_cor = prob_pdf(m_c, cov_c)
    pdf_decor = prob_pdf(m_d, cov_d)

    #U_c_0 = U(m_c, pdf_cor, T)
    #U_d_0 = U(m_d, pdf_decor, T)
    #U_c = lambda x: U(x, pdf_cor, T) - U_c_0
    #U_d = lambda x: U(x, pdf_decor, T) - U_d_0

    R=0.008314463 #kJ/mol
    #f_c = lambda x: x[:,0]*x[:,0]*np.sin(x[:,1])*np.sin(x[:,2])*np.exp(-U_c(x)/(R*T))
    f_c = lambda x: x[:,0]*x[:,0]*np.sin(x[:,1])*np.sin(x[:,2])*(pdf_cor(x)/pdf_cor(m_c))
    #f_d = lambda x: x[:,0]*x[:,0]*np.sin(x[:,1])*np.sin(x[:,2])*np.exp(-U_d(x)/(R*T))

    #since we use rigid rotator approx. in the analytical correction in summary we need to use it here as well.
    f_d = lambda x: m_d[0]*m_d[0]*np.sin(m_d[1])*np.sin(m_d[2])*(pdf_decor(x)/pdf_decor(m_d))

    dG=[]
    for i in range(reps):
        Q_c=MCI(f_c, minlim_c, maxlim_c, N)
        Q_d=MCI(f_d, minlim_d, maxlim_d, N)
        dG.append(-R*T*np.log(Q_c/Q_d))
#     print("dG =", np.mean(dG), "+-", np.std(dG), "kJ/mol")
    return(np.mean(dG), np.std(dG)/np.sqrt(reps))


# ==============================================================================
#                         Derivative Task Classes
# ==============================================================================
class Task_PL_energetic_decorrelation(SGETunedJobTask):

    #Parameters:
    p = luigi.Parameter(description='Protein name')
    l = luigi.Parameter(description='Ligand name')
    i = luigi.IntParameter(description='Repeat number')

    folder_path = luigi.Parameter(significant=False,
        visibility=ParameterVisibility.HIDDEN,
        description='Path to the protein+ligand folder to set up')

    study_settings = luigi.DictParameter(significant=True,
        visibility=ParameterVisibility.HIDDEN,
        description='Dict of study stettings '
        'used to propagate settings to dependencies')

    #request 1 core
    n_cpu = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=1, significant=False)

    job_name_format = luigi.Parameter(
        visibility=ParameterVisibility.HIDDEN,
        significant=False, default="pmx_{task_family}_p{p}_l{l}_{i}",
        description="A string that can be "
        "formatted with class variables to name the job with qsub.")

    T = luigi.FloatParameter(default=298.0, significant=True,
        description="Simulation temperature.")

    MCI_samples = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=int(1e7), significant=False,
                               description="Number samples for MC integration.")

    MCI_repeats = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=100, significant=False,
                               description="Number times to repeat MC integration for statistics.")
    MCI_limits_nsigma = luigi.FloatParameter(default=4.5, significant=True,
                               description="How many standard deviations away from the distribution mean "
                               "to set the integration limits. Should be between 4 and ~7 with MCI_samples=1e7."
                               "Too low value leads to a significant potrion of the distribution "
                               "outside integration limits. Too large leads to a lot of noise and need for more MCI_samples.")
    force_redo_MCI = luigi.BoolParameter(default=False, significant=False,
                               description="Run MCI even if it was run before.")

    n_bins = luigi.IntParameter(visibility=ParameterVisibility.HIDDEN,
                               default=10, significant=False,
                               description="Number of histogram bins in plot.")

    n_bootstrap = luigi.IntParameter(default=0, significant=True,
                               description="Number of times estimators are bootstrapped.")


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        #set variables
        self.base_path = self.study_settings['base_path']
        self.decor_folder=self.folder_path+"/decor_analysis/repeat%d"%self.i
        self.ana_folder=self.folder_path+"/analysis/repeat%d"%self.i

        #make sure numpy doesn't use more threads than we request from sge
        if(not self.run_locally):
            self.qsub_extra_options+=" -v OMP_NUM_THREADS={}".format(self.n_cpu)
            self.qsub_extra_options+=" -v NUMEXPR_NUM_THREADS={}".format(self.n_cpu)
            self.qsub_extra_options+=" -v MKL_NUM_THREADS={}".format(self.n_cpu)


    def work(self):
        if(self.study_settings['n_sampling_sims']!=1):
            raise(Exception("Currently only supports n_sampling_sims=1."))

        os.makedirs(self.decor_folder, exist_ok=True)
        distrib_cache_fn=self.decor_folder+"/distribution_cache.pickle"

        if(os.path.isfile(distrib_cache_fn)):
            m,cov,m_decor,cov_decor = pickle.load( open( distrib_cache_fn, "rb" ) )
        else:
            #step 1: find the covarriance from aligned trajectory
            ii=self.folder_path+f"/ii_C_{self.i}.itp"
            if(not os.path.isfile(ii)):
                ii=self.folder_path+f"/ii_{self.i}.itp" #backwards compatibility for older runs
            aligned_path=self.folder_path+f"/stateC/repeat{self.i}/morphes{0}/"
            struct=aligned_path+"/frame0.gro"
            traj=aligned_path+"/aligned.trr"

            if(not os.path.isfile(struct)):
                raise(Exception(f"First aligned frame file {struct} is missing!"))
            if(not os.path.isfile(traj)):
                raise(Exception(f"Aligned trajectory file {traj} is missing!"))
            nframes = len(glob.glob1(aligned_path,"frame*.gro"))
            m, cov = build_aligned_distib_no_data(ii, struct, traj, nframes)


            #step 2: read the covarriance for uncorrelated restraints
            idx_lig, idx_pro, m_decor, ks_decor = readii(ii)
            R=0.008314463 #kJ/mol K
            var_decor=R*self.T/ks_decor #FC=kT/(var) -> var=kT/FC # still in rad^2
            var_decor[1:]*=(180./np.pi)**2 #convert to deg^2
            cov_decor=np.diag(var_decor)

            #save distribution cache
            pickle.dump( (m,cov,m_decor,cov_decor), open( distrib_cache_fn, "wb" ) )


        corr_cache_fn=self.decor_folder+"/correction_cache.pickle"
        dG_dif=None
        err_dif=None
        run_MCI=True
        if(os.path.isfile(corr_cache_fn) and not self.force_redo_MCI):
            dG_dif, err_dif, saved_MCI_samples, saved_MCI_repeats, saved_MCI_limits_nsigma = pickle.load( open( corr_cache_fn, "rb" ) )
            #if MCI setting were the same, no need to rerun MCI
            if(saved_MCI_samples==self.MCI_samples and saved_MCI_repeats==self.MCI_repeats and saved_MCI_limits_nsigma==self.MCI_limits_nsigma):
                run_MCI=False
        if(run_MCI):
            #step 3: find free energy difference between correlated and decorrelated ensembles
            dG_dif, err_dif = dG_by_MCI_between_two_distribs(m, cov,
                                    m_decor, cov_decor,
                                    self.T,
                                    N=self.MCI_samples, reps=self.MCI_repeats,
                                    nsigma=self.MCI_limits_nsigma)
            pickle.dump( (dG_dif, err_dif,self.MCI_samples,self.MCI_repeats,self.MCI_limits_nsigma), open( corr_cache_fn, "wb" ) )

        if(dG_dif is None or err_dif is None):
            raise(Exception("This should never happen."))

        print(f"Energetic correction for correlation of restrained DOFs: {dG_dif}+-{err_dif}")
        if(err_dif>0.2):
            print("Warning: MCI uincertainty (SEM) for the correction is fairly high. Consider increasing MCI_samples or MCI_repeats and reruning.")

        #step 4: make new integrated W files with shifted reverse W
        shutil.copyfile(self.ana_folder+"/integA.dat", self.decor_folder+"/integA.dat")

        with open(self.ana_folder+"/integB.dat",'r') as fin:
            origB = fin.readlines()

        with open(self.decor_folder+"/integB.dat",'w') as fout:
            for l in origB:
                s=l.split()
                if(len(s)>0):
                    W = float(s[1]) - dG_dif
                    newl = "{} {}\n".format(s[0], W)
                    fout.write(newl)

        #step 5: run analysis script on new W distributions
        os.chdir(self.decor_folder)
        orig_argv  =sys.argv
        orig_stdout=sys.stdout
        orig_stderr=sys.stderr
        sys.argv = [['analyze_dhdl.py'],
                    ['-iA', self.decor_folder+"/integA.dat"],
                    ['-iB', self.decor_folder+"/integB.dat"],
                    ['--nbins', str(int(self.n_bins))],
                    ['-b', str(int(self.n_bootstrap))]]
        sys.argv = [item for sublist in sys.argv for item in sublist] #flatten argv

        with open("analysis.log", "w") as f:
            sys.stdout = f
            sys.stderr = f
            pmx.scripts.analyze_dhdl.entry_point()

        sys.argv  =orig_argv #reset argv
        sys.stdout=orig_stdout #reset output
        sys.stderr=orig_stderr

        os.chdir(self.base_path)#reset cwd

    def output(self):
        files=['wplot.png', 'analysis.log', 'integA.dat', 'integB.dat', 'results.txt']
        return([luigi.LocalTarget(os.path.join(self.decor_folder, f)) for f in files])

    def complete(self):
        """
        Check if Bootstraping was done (if required).
        """
        c=super().complete()
        if(c and self.n_bootstrap>0):
            found_boots=0
            with open(os.path.join(self.decor_folder, 'analysis.log'), "r") as anaf:
                content = anaf.readlines()
                for l in content:
                    if("Bootstrap (Std Err): iteration" in l):
                        s=l.split('/')
                        found_boots=int(s[-1])
                        break;
            if(found_boots!=self.n_bootstrap):
                c=False

        return(c)

    def requires(self):
        tasks=[]
        tasks.append(Task_PL_analysis_aligned(
              p=self.p, l=self.l, i=self.i,
              study_settings=self.study_settings,
              folder_path=self.folder_path,
              parallel_env=self.parallel_env) )
        
        #tasks.append( Task_PL_gen_morphes(p=self.p, l=self.l,
                        #i=self.i, m=0, sTI="C",
                        #study_settings=self.study_settings,
                        #folder_path=self.folder_path,
                        #parallel_env=self.parallel_env,
                        #restr_scheme="Aligned") )

        return(tasks)

