#! /bin/bash
source /usr/local/gromacs/GMXRC2020   # will use gmx 2020 for setup, but for simulations will use gmx selected in $MDR and $MDRD
export GMXLIB=/home/vgapsys/gromacs-454/share/gromacs/top # replace with your own
module load sge

cwd=$(pwd)

# wrapers to use specific gromacs versions during simulations that can be different from the one used for setup.
MDRD=$cwd/../mdrun_optimal_double_2019.4_fixed.sh
MDR=$cwd/../mdrun_optimal_2019.4_fixed.sh

# to run with a remote luigi central scheduler (it has a gui for monitoring progress):
#python ./Workflow_aligned_Complete.py --mdrun $MDR  --mdrun_double $MDRD --mdrun_opts=" -ntmpi 1 -notunepme" --toppath ./data --mdppath ./data/mdp --rem_sched

# to run locally without a remote scheduler
python ./Workflow_aligned_Complete.py --mdrun $MDR  --mdrun_double $MDRD --mdrun_opts=" -ntmpi 1 -notunepme " --toppath ./data --mdppath ./data/mdp --workers 8
