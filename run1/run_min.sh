#!/bin/bash
#source /wynton/home/grabe/shared/gromacs/gromacs-2020.6_CUDA10_SSE4/bin/GMXRC
source /wynton/home/rotation/jessicaflowers/software/gromacs/2025.4_CUDA/bin/GMXRC

initial_gro="../input/step3_input.gro"
top="../input/topol.top"
ndx="../input/index.ndx"
mdp="../mdp/seg_01_MIN.mdp"
out="seg_01_MIN"

gmx grompp -f $mdp -c $initial_gro -r $initial_gro -p $top -n $ndx -o $out.tpr
gmx mdrun -v -deffnm $out -ntmpi 1

#gmx mdrun -v -deffnm $out
