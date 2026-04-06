# Modified RIP (Rotamerically Induced Perturbation) Procedure

The pipeline has two steps:
1. **rip_perturb.py** — generates perturbed .gro files (one per residue)
2. **submit.sh** — submits MD production runs to the SGE cluster


## Directory structure

    modified_RIP/
    ├── pdbstruct/              # shared Python library (chi topology, atom classification)
    ├── 1hck/                   # protein 1 (CDK2 kinase)
    │   ├── input/              # topol.top, index.ndx, step3_input.gro, toppar/
    │   ├── mdp/                # seg_11_NPT_PROD.mdp (production MDP)
    │   └── run2/               # working directory
    │       ├── rip_perturb.py
    │       ├── run_dynamics.sh
    │       ├── submit.sh
    │       ├── seg_11_NPT_PROD_01.gro   # equilibrated starting structure
    └── 3g33/                   # protein 2 (CDK4 kinase)
        ├── input/
        ├── mdp/
        └── run2/
            └── (same layout as 1hck/run2/)


## Step 1: Generate perturbed structures with rip_perturb.py

Run from the run2/ directory. The script reads the equilibrated .gro file, applies a single chi-angle velocity pulse to the specified residue(s), and writes a perturbed .gro into a pulse_res_<num>/ subdirectory.

### Usage

    python rip_perturb.py [gro_base] [topology] [temperature] [residue_numbers]

### Arguments

    gro_base          Base name of the starting .gro file (without .gro extension).
                      Default: seg_11_NPT_PROD_01

    topology          Path to the GROMACS topology file.
                      Default: ../input/topol.top

    temperature       (Optional) Perturbation temperature in Kelvin. Controls how much
                      kinetic energy the chi-angle pulses carry. Should be
                      higher than the simulation temperature (310.15 K) to
                      produce a meaningful perturbation.
                      Default: 350.0

    residue_numbers   (Optional) Which residue(s) to perturb. Can be a single
                      number or a comma-separated list. If omitted, ALL
                      rotatable residues are perturbed (one directory each).

### Examples

Perturb a single residue (residue 203 at 350 K):

    python rip_perturb.py seg_11_NPT_PROD_01 ../input/topol.top 350 203

Perturb two specific residues:

    python rip_perturb.py seg_11_NPT_PROD_01 ../input/topol.top 350 203,131

Perturb ALL rotatable residues using defaults:

    python rip_perturb.py

This will create a pulse_res_<num>/ directory for every residue that has at least one chi angle. Residues that already have an existing pulse_res_<num>/ directory are skipped automatically.

### What it creates

For each perturbed residue, the script creates:

    pulse_res_<num>/
    ├── seg_11_NPT_PROD_01_PERT.gro     # perturbed starting structure
    └── topol.top                       # copy of the topology


## Step 2: Run MD with submit.sh

submit.sh is an SGE array job script that finds pulse_res_<num>/ directories
and runs run_dynamics.sh inside each one.

### Usage

From the run2/ directory:

    qsub submit.sh

### How it works

- Each SGE array task (1-150) searches for the next pulse_res_<num>/ directory starting from its task ID
- Inside that directory, it runs run_dynamics.sh which:
  1. Runs gmx grompp to prepare the .tpr
  2. Runs gmx mdrun for production MD
  3. Loops PROD_REPEAT times (default: 10), each time continuing from the previous iteration's output
  4. If a task gets killed at the walltime limit (2 hours), the next array task picks up from the checkpoint (.cpt) file

### Key parameters in run_dynamics.sh

    PROD_REPEAT=10     # number of sequential MD chunks
    NTHREAD=8          # CPU threads
    NTOMP=8            # OpenMP threads

### Key parameters in mdp/seg_11_NPT_PROD.mdp

    dt       = 0.004          # 4 fs timestep (with HMR factor 3)
    nsteps   = 500000         # 2 ns per chunk (x10 repeats = 20 ns total)
    ref-t    = 310.15         # simulation temperature (K)
    tcoupl   = v-rescale      # thermostat (turining this off is necessary for RIP bus also crashed the simulation for some reason ... need to debug)
    pcoupl   = c-rescale      # barostat (turining this off is necessary for RIP bus also crashed the simulation for some reason ... need to debug)
    nstxout-compressed = 25000  # write frame every 100 ps

### Notes

Steps to create the /input and /mdp folders if you wish to start from scratch:

Use solution builder from CHARMM-GUI for the initial set up. It should generate a bunch of folders and files with the general format:

    charmm-gui/
    ├── tons of input files (ie .pdb, .cif, .in, .crd, etc.)
    ├── gromacs/              
        ├── README
        ├── step3_input.gro
        ├── other important files
        ├── THIS IS WHERE YOU SHOULD COPY dynamics-parameter-generation.py
    ├── toppar/                   
    │   ├── so many files, not important to discuss

When this is complete, you will want to take the dynamics-parameter-generation.py from this repo, copy it into the \gromacs folder, and run it. This will generate the /input and /mdp folders that are required to run the MD simulations.  