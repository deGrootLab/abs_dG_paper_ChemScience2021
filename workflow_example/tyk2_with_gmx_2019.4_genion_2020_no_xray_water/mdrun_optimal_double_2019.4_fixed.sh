#!/bin/bash

source /home/ykhalak/soft/gromacs-2019.4_debug_fixed/InstalledHere_owl/bin/GMXRC

# Comment out this block if you do not have hardware optimised compilations of mdrun
######################################################################################
CPUFLAGS=$( cat /proc/cpuinfo | grep flags     | head -1 )
VENDORID=$( cat /proc/cpuinfo | grep vendor_id | head -1 )

if     [ $( echo $VENDORID | grep -c "AuthenticAMD" ) -eq "1" ]; then
    ACCEL="_AVX2_128"
elif   [ $( echo $CPUFLAGS | grep -c "avx2"         ) -eq "1" ]; then
    ACCEL="_AVX2_256"
elif   [ $( echo $CPUFLAGS | grep -c "avx"          ) -eq "1" ]; then
    ACCEL="_AVX_256"
elif   [ $( echo $CPUFLAGS | grep -c "sse4_1"       ) -eq "1" ]; then
    ACCEL="_SSE4_1"
elif   [ $( echo $CPUFLAGS | grep -c "sse2"         ) -eq "1" ]; then
    ACCEL="_SSE2"
else
    ACCEL=""
fi
######################################################################################

MDRUN=mdrun_threads_d$ACCEL

$MDRUN $@
