# abs_dG_paper_data
Structures and data for https://doi.org/10.1039/D1SC03472C 

# Software
Throughout this work a modified version of gromacs 2019.4 was used, fixing a bug (gmx issue 3403) with how the free energy kernel handled exclusions beyond the electrostatic cutoff. The code is available in the fix_for_gmx_2019.4 folder.

Most of the absolute free energy calculations were carried out with an automated workflow available separately here: https://github.com/deGrootLab/pmx/tree/abs_dG_workflow/src/pmx/scripts/workflows
A manual and a working example of the workflows use is available in the workflow_example folder.

# Starting Structures and simulation settings
These are found in the other folders herein, named by the corresponding protein. These folders also include text files with the final results for the calcualtions.
