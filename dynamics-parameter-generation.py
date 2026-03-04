###################################################################################################################
# Adapted from Jonathan Borowsky
# This program edits charmmgui output files and generates new files to execute a particular equilibration protocol 
# In the end, the input/ and mdp/ directories are ones that can be used to run the simulations.

# 1. MINIMIZATION MDP
# 
# 2. EQUILIBRATION MDP
# Containing the following information:
#   length (ps)
#   thermostat, if any (enter None otherwise)
#   barostat, if any (enter None otherwise)
#       if neither a thermostat nor barostat is specified, run minimization
#       if only a thermostat is specified, the system is in an NVT ensemble\
#       if both a thermostat and barostat are specified, the system is in an NPT ensemble
#   the electric field vector 
#   a list of atom group names
#   a matching list of the restraint strengths for each atom group (restraints are assumed to be isotropic, as supporting anisotropic restraints would involve more topology .itp file editing)
# 3. PRODUCTION MDP
# Same as equilibration MDP, but with different time step (dt=0.004 ps) and mass repartition factor (hmr_factor=3)
###################################################################################################################3
import os
import sys
import datetime

#----------------------------------------------------------------------------------
#                               MDP FILE PARAMETERS
#----------------------------------------------------------------------------------

# This section contains lists of the parameters common to all or most minimization equilibration runs,
# as well as methods to fill in the variable parameters and return lists in the same format


#see https://manual.gromacs.org/current/user-guide/mdp-options.html

# Path to the Gromacs folder
GROMACS_FOLDER = "../gromacs"

# Physical parameters for thermostat and barostat
TEMPERATURE = 310.15 #K
PRESSURE = 1.01325 #bar (1 atmosphere)
WATER_COMPRESSIBILITY = 4.5*10**-5 #isothermal compressibility of water, bar^-1, at 300K and 1 atmosphere
                                   #We should either be using the compressibility at 310.15K or running at 300 

# Time step, which is important for determining when to split segments and so goes out here
dt = 0.002 #time step in femtoseconds, should be 0.002 for dynamics

#------------------ Below are sections of the mdp file -------------------
minimization_integration_parameters = [
    ["integration"           ],
    ["integrator",           "steep"],              
    ["emstep",               0.001],
    ["nsteps",               10000],
    ["emtol",                1000]                    #kJ/mol/nm; Stop minimization early if the maximum force is less than emtol
]

def dynamics_integration_parameters(dt, nsteps, tinit, hmr_factor):
    
    #nsteps = int(round(length/dt))
    return [
        ["integration"           ],
        ["integrator",           "md"],               #
        ["dt",                   dt],                 #time step in femtoseconds, should be 0.002 for dynamics
        ["nsteps",               nsteps],             #number of steps
        ["tinit",                tinit],              
        ["mass-repartition-factor",           hmr_factor]           
        ]

#use these only for dynamics runs
dynamics_constraint_parameters = [
    ["bonded interactions"   ],
    ["constraints",          "h-bonds"],
    ["constraint-algorithm", "lincs"],
    ["lincs-order",          4],                      #a value of 8 should be used for minimization in double precision, but values greater than 4 should not be used in single precision
    ["lincs-iter",           1]                       #a value of 2 should be used for minimization in double precision, but values greater than 1 should not be used in single precision
]

#Note: do not use a hard cutoff with a monte carlo barostat
#   see Yessica's paper on problematic barostat-cutoff combinations for more information: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8776593/
def nonbonded_force_calculation_parameters(electric_field_vector):
    nb_params = [
        ["neighbor searching"    ],
        ["cutoff-scheme",        "Verlet"],               #
        ["nstlist",              20],                     #
        ["pbc",                  "xyz"],                  #periodic boundary conditions along all three axes
        
        ["electrostatics"        ],
        ["coulombtype",          "PME"],                  #use PME electrostatics
        ["pme-order",            4],                      #order of the PME interpolation
        ["fourierspacing",       0.12],                   #PME grid spacing
        ["rcoulomb",             1.2],                    #minimum electrostatic cutoff radius

        ["Van der Waals"         ],
        ["vdwtype",              "cut-off"],              #assume van der waals forces are 0 beyond the cutoff distance
        ["vdw-modifier",         "force-switch"],         #make the force go smoothly to 0 between rvdw and rvdw-switch
        ["rvdw-switch",          1],                      #see above
        ["rvdw",                 1.2]                     #van der waals cutoff radius

    ]
    
    if electric_field_vector != [0,0,0,0]:
        nb_params += [
            ["external electric field"                ],
            ["electric-field-z", electric_field_vector]]
    
    return nb_params

minimization_output_parameters = [
    ["output control; frequencies in steps"],
    ["nstxout",              50],                     #positions to uncompressed log file
    ["nstenergy",            10]                      #energies to energy file

]

#comments refer to what is being written and to where
def dynamics_output_parameters(log_frequency):
    return [
        ["output control; frequencies in steps"],
    #    ["nstxout",              log_frequency],          #positions to uncompressed log file
    #    ["nstvout",              log_frequency],          #velocities to uncompressed log file
        ["nstlog",               log_frequency],          #energies to uncompressed log file
        ["nstenergy",            log_frequency],          #energies to energy file
        ["nstxout-compressed",   log_frequency]           #positions to compressed log file
    ]

#for equilibration with restraints
def restraint_parameters(restraint_groups, restraint_strengths):

   # restraint_string = "-DPOSRES " + " ".join([f"-DPOSRES_FC_{grp}={fc}" for (grp, fc) in zip(restraint_groups, restraint_strengths)])

    restraint_string = "-DPOSRES " + " ".join(
    ["-DPOSRES_FC_{}={}".format(grp, fc) for grp, fc in zip(restraint_groups, restraint_strengths)])


    return [
        ["restraints"            ],
        ["define",               restraint_string]                      #restraints
        ]

#for equilibration after restraints are lifted
#comm removal is not done when restraints are present because restraints already prevent 
# accumulation of net momentum by the system
#the gromacs file checker will print a note about this regardless of whether your mdp file has comm removal
# see https://gromacs.bioexcel.eu/t/com-motion-and-position-restraints-in-an-npt-ensemble/1778/4 
unrestrained_comm_parameters = [
    ["center of mass motion removal"],
    ["comm-mode",            "linear"]                #what type of center of mass motion to remove; 
                                                      #periodic systems do not require angular COM motion removal because the PBC already prevents this
#    ["nstcomm",              100]                     #frequency of center of mass motion removal
]

#see https://manual.gromacs.org/documentation/nightly/reference-manual/algorithms/molecular-dynamics.html group temperature coupling section
#see 37 M.P. Eastwood, K.A. Stafford, R.A. Lippert, M.Ø. Jensen, P. Maragakis, C. Predescu, R.O. Dror, and D.E. Shaw, “Equipartition and the calculation of temperature in biomolecular simulations,” J. Chem. Theory Comput., ASAP DOI: 10.1021/ct9002916 (2010).
#see https://manual.gromacs.org/current/user-guide/terminology.html 
#current relevance of multiple TC groups is unclear: https://gromacs.bioexcel.eu/t/multiple-tc-grps-does-make-any-sense/6129
#It might also be good to have multiple groups in early equilibration stages but not in production: https://simtk.org/tracker/index.php?func=detail&aid=886&group_id=161&atid=436 
#see also https://arxiv.org/pdf/cond-mat/0511629 for theory
def thermostat_parameters(thermostat):
    return [
        ["thermostat"            ],
        ["tcoupl",               thermostat],                     #which thermostat to use
        ["tc_grps",              "System"],                         #groups of atoms to which to apply thermostat
        ["tau-t",                1],            #frequency of temperature coupling
        ["ref-t",                TEMPERATURE]   #target temperature
]

def barostat_parameters(barostat, nstpcouple=-1):
    barostat_parameters = [
        ["barostat"             ],
        ["pcoupl",              barostat],                   #which barostat to use
        ["pcoupltype",          "isotropic"],            #one barostat for xy and one for z due to membrane in xy
        ["tau-p",               5],                          #pressure coupling frequency
        ["ref-p",               [PRESSURE]],       #target pressure
        ["compressibility",     [WATER_COMPRESSIBILITY]],
        ["refcoord-scaling",    "com"]                       #keep the center of mass of the reference foodinates centered but do not rescale them; important when unrestrained waters may change the box size and shape but protein is still restrained
    ]
    if nstpcouple != -1:
        barostat_parameters += [["nstpcouple",           nstpcouple]]                       #frequency of pressure coupling
    return barostat_parameters

#for the first dynamics run after minimization
dynamics_init_parameters = [
    ["input handling"         ],
    ["gen-vel",              "yes"],                  #whether to generate new velocities from a maxwell boltzmann distribution. 
                                                      #Only do this if not loading from a .trr file at the desired temperature
    ["gen-temp",             TEMPERATURE],            #the temperature of the MB distribution from which veloci
    ["gen-seed",             -1],                     #seed for random velocity generation; -1 is random
    ["continuation",         "no"]                    #apply constraints to the starting configuration
]

#for subsequent dynamics runs which continue other dynamics runs
dynamics_continuation_parameters = [
    ["input handling"           ],
    ["gen-vel",              "no"],                   #whether to generate new velocities from a maxwell boltzmann distribution. 
    ["continuation",         "yes"]                   #do not apply constraints to the starting configuration because they were already applied in the previous trajectory segment
]

#----------------------------------------------------------------------------------
#                    TABLE OF EQUILIBRATION SEGMENT PARAMETERS
#----------------------------------------------------------------------------------

# This section contains a table of the equilibration parameters \
# which vary between different equilibration rounds and different systems


# in each row:  [[0]segment length in ps (except for minimization, where length is specified above and does not have units of time), 
#               [1]thermostat, 
#               [2]barostat, 
#               [3]electric field parameters. One can make a constant electric field using [E, 0, 0, 0], where E is in V/nm; 
#                      see https://manual.gromacs.org/current/user-guide/mdp-options.html#electric-fields, 
#                      and https://manual.gromacs.org/documentation/current/reference-manual/special/electric-fields.html
#               [4]list of gromacs atom group names,
#               [5]list of harmonic restraint spring coefficients for each atom group, in kJ/(mol*nm^2)]


prot_fn = ["PROC"] # protein subunits
solvent_groups = ["SOD", "CLA", "TIP3"] # water, ions, lipid
small_molecules = [] # ligand

equilibration_segment_parameters = [
    # Minimization
    [0, None, None, [0,0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [4000, 4000, [4000 for _ in range(len(small_molecules))], [400 for _ in range(len(solvent_groups))]]],
    # NVT
    [100, "v-rescale", None, [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [4000, 4000, [4000 for _ in range(len(small_molecules))], [400 for _ in range(len(solvent_groups))]]],
    # NPT 1
    [100, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [4000, 4000, [4000 for _ in range(len(small_molecules))], [400 for _ in range(len(solvent_groups))]]],
    # NPT 2
    [100, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [2000, 2000, [2000 for _ in range(len(small_molecules))], [80 for _ in range(len(solvent_groups))]]],
    # NPT 3
   # [100, "v-rescale", "c-rescale", [0.0, 0,0,0], 
   #  ["BB", "SC"] + small_molecules + solvent_groups, 
   #  [2000, 2000, [2000 for _ in range(len(small_molecules))], [16 for _ in range(len(solvent_groups))]]],
    # NPT 4: unrestrain solvent
    [5000, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [2000, 2000, [2000 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    # NPT 5: release side chain restraints
    [500, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [2000, 1000, [2000 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    # NPT 6: release backbone, side chain, and ligand restraints
    [500, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [1000, 400, [1000 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    #NPT 7: release backbone, side chain, and ligand restraints
   # [500, "v-rescale", "c-rescale", [0.0, 0,0,0], 
   #  ["BB", "SC"] + small_molecules + solvent_groups, 
   #  [400, 200, [400 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    #NPT 8: release backbone, side chain, and ligand restraints
    [500, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [200, 20, [200 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    #NPT 9: unrestrained sidechain
    [5000, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [20, 0, [20 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    #NPT 10: unrestrained everything
    [10000, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [0, 0, [0 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]],
    #NPT 9
    # This is where you add electric field.
    # PROD: production run. This will be run X times.
    [100000, "v-rescale", "c-rescale", [0.0, 0,0,0], 
     ["BB", "SC"] + small_molecules + solvent_groups, 
     [0, 0, [0 for _ in range(len(small_molecules))], [0 for _ in range(len(solvent_groups))]]]]

def flatten_restraints(restraint_list):
    flat_list = []
    for item in restraint_list:
        if isinstance(item, list):
            flat_list.extend(item)
        else:
            flat_list.append(item)
    return flat_list

#----------------------------------------------------------------------------------
#                           ASSEMBLE AND WRITE MDP FILES
#----------------------------------------------------------------------------------

#copy inputs out of charmmgui output automatically
#unless an input directory already exists
def extract_inputs():

    #make input folder
    if not os.path.exists("input"):
        os.system("mkdir input")

        os.system(f"cp {GROMACS_FOLDER}/index.ndx input")
        os.system(f"cp {GROMACS_FOLDER}/step3_input.gro input")
        os.system(f"cp {GROMACS_FOLDER}/topol.top input")
        os.system(f"cp -r {GROMACS_FOLDER}/toppar input")

    #for manually prepared inputs
    else:
        print(f"using existing inputs: {os.listdir('input')}. Input directory should contain:\n    index.ndx  step5_input.gro  topol.top  toppar/")


# Remove existing position restraint information from the topology.itp files in toppar/ and add new position restraints.
# The new restraints harmonically and isotropically restrain heavy atoms with a force constant defined in the .mpp file
# Charmmgui has already provided the desired restraints for proteins (and these are more complicated due to the backbone-sidechain distinction), 
# so the PROA.itp file is not edited

def edit_topology_itp_files(small_molecules, solvent_groups, prot_name):
    top_itp_files = os.listdir("input/toppar")

    #check that the list of small molecules above matches the ones in toppar/
    known_molecules = set(prot_name) | set(small_molecules) | set(solvent_groups)
    known_molecules.add("forcefield")

    for fn in top_itp_files:
        if not fn.endswith(".itp") or fn.endswith("-original.itp"):
            continue
        mol_name = fn[:-4]
        print(mol_name)
        if mol_name not in known_molecules:
            print(f"warning: no restraints given for {fn}")
   
    #edit the topology .itp file for each small molecule
    for mol in small_molecules + solvent_groups:
        #rename original .itp file as a backup/reference copy
        if not os.path.exists(f"input/toppar/{mol}-original.itp"):
            os.system(f"mv input/toppar/{mol}.itp input/toppar/{mol}-original.itp")
        
        #open both the reference file and the new .itp file
        with open(f"input/toppar/{mol}-original.itp", 'r') as file_in:
            with open(f"input/toppar/{mol}.itp", 'w') as file_out:
            
                #variables to track the reader's state
                
                #When the first line not belonging to the existing header is detected, 
                # additional header information is written and write_header is 
                # set to False to ensure no duplicate header is written
                write_header = True

                #gets set to True when the [atoms] section header is detected, 
                # then False again when the next section header is detected.
                # while it is true atom information is read into heavy_atom_inds
                reading_atoms = False
                heavy_atom_inds = []

                #whether to copy the current line to the input file. 
                # Gets set to False when existing position restraints are detected 
                # and to True again after the end of those position restraints,
                # causing those restraints not to be copied
                copylines = True
                
                # read through each line in the input file. 
                # most lines are written unchanged to the output file. 
                # lines to be deleted are not written to the output file.
                # additional lines are added to the output file in places defined by the input file contents
                
                for line in file_in:

                    #write header following the existing charmmgui header
                    if line[0:2] != ";;" and write_header:
                        file_out.write(
                        f"\n;; This gromacs .itp file was edited\n\
;; as {os.path.abspath(f'{mol}.itp')}\n\
;; by {os.path.realpath(__file__)}\n\
;; on {datetime.datetime.now()}\n\n"
                        )
                        write_header = False

                    #identify end of atoms section
                    if line[0] == "[" and line.strip() != "[ atoms ]":
                        reading_atoms = False

                    #read [atoms] section to identify heavy atoms
                    if reading_atoms:
                        if line.strip() != "":
                            linedata = list(filter(None, line.strip().split(" ")))
                            if float(linedata[7]) != 1.008:
                                heavy_atom_inds.append(int(linedata[0]))

                    #identify start of [atoms] section
                    if line.strip() == "; nr	type	resnr	residu	atom	cgnr	charge	mass": #"[ atoms ]":
                        reading_atoms = True

                    #identify start of existing position restraints and disable writing to output file
                    if line.strip() == "#ifdef POSRES":
                        copylines = False

                    #copy input file contents to output file
                    if copylines:
                        file_out.write(line)

                    #identify end of existing position restraints and reenable writing to output file                    
                    if line.strip() == "#endif":
                        copylines = True

                #write new position restraints
                file_out.write("\n")
                file_out.write("#ifdef POSRES\n")
                file_out.write("[ position_restraints ]\n")

                for hai in heavy_atom_inds:
                    #the first column is the atom index, IDK what the second column is, and the other three are harmonic restraint force constants
                    file_out.write(f"{str(hai).rjust(5)}     1    POSRES_FC_{mol}    POSRES_FC_{mol}    POSRES_FC_{mol}\n")
                file_out.write("#endif\n")

                file_out.close()
            file_in.close()


#write a single mdp file using a list of mdp file parameters for one equilibration segment
def write_mdp_file(filename, parameter_list):

    f = open("mdp/" + filename, "w") #use x in the future

    #write header
    #TODO: make this print the server and user name
    f.write(
f"; This gromacs .mdp file was generated\n\
; as {os.path.abspath(filename)}\n\
; by {os.path.realpath(__file__)}\n\
; on {datetime.datetime.now()}\n"
)
    
    #write parameters and comments
    for p in parameter_list:
        
        #write comments
        if len(p) == 1:
            row = "\n; " + p[0]
        
        #write parameters which are lists of values 
        elif type(p[1]) == list:
            row = p[0].ljust(25) + "= " + " ".join([str(i) for i in p[1]])
        
        #write parameters which are single values         
        else:
            row = p[0].ljust(25) + "= " + str(p[1])
        
        f.write(row + "\n")
    

# Assemble complete lists of mdp file parameters for each equilibration segment and write the corresponding mdp files
def assemble_mdp_file_inputs(segment, i, dt, n_steps, t_init, hmr_factor, log_frequency, run_type):
    
    i_num = str(i+1).zfill(2)

    #control inclusion of restraints and restraint-dependent parameters, 
    # which are orthogonal to the thermostat and barostat,
    # although in practice only NPT simulations are run without restraints for my purposes
    flat_restraints = flatten_restraints(segment[5])
    if all(r == 0 for r in flat_restraints):
        restraint_dependent_parameters = unrestrained_comm_parameters
    else:
        restraint_dependent_parameters = restraint_parameters(segment[4], flat_restraints)

    # Minimization (non-minimization NVE is not currently supported)
    if segment[1] is None and segment[2] is None:
        write_mdp_file(f"seg_{i_num}_{run_type}.mdp", 
                  minimization_integration_parameters 
                + nonbonded_force_calculation_parameters(segment[3])                 
                + restraint_dependent_parameters
                + minimization_output_parameters)
    
    # NVT ensemble
    # Warning: dynamics_init_parameters is used here, implicitly assuming a single NVT segment immediately after equilibration
    elif segment[2] is None:
        write_mdp_file(f"seg_{i_num}_NVT_{run_type}.mdp", 
                  dynamics_init_parameters
                + dynamics_integration_parameters(dt, n_steps, t_init, hmr_factor) 
                + dynamics_constraint_parameters
                + nonbonded_force_calculation_parameters(segment[3]) 
                + restraint_dependent_parameters
                + thermostat_parameters(segment[1])
                + dynamics_output_parameters(log_frequency))
    
    # NP ensemble is not supported; I've never heard of anyone using one
    elif segment[1] is None:
        print("error: NP ensemble not supported")
        sys.exit(0)

    # NPT ensemble
    else:
        if i == 2: # NPT3
            nstpcouple = 5 # more accurate pressure coupling for first step of NPT 
        else:
            nstpcouple = -1
        write_mdp_file(f"seg_{i_num}_NPT_{run_type}.mdp", 
                  dynamics_continuation_parameters
                + dynamics_integration_parameters(dt, n_steps, t_init, hmr_factor) 
                + dynamics_constraint_parameters
                + nonbonded_force_calculation_parameters(segment[3]) 
                + restraint_dependent_parameters
                + thermostat_parameters(segment[1])
                + barostat_parameters(segment[2], nstpcouple)
                + dynamics_output_parameters(log_frequency))
    
    return t_init + dt*n_steps


def write_all_mdp_files(equilibration_segment_parameters, num_prod_runs):

    if not os.path.exists("mdp"):
        os.system("mkdir mdp")

    total_dynamics_length = 0 #in picoseconds

    for i, segment in enumerate(equilibration_segment_parameters):
        if (segment[1] is None and segment[2] is None): # Minimization segment
            dt = 0
            log_frequency  = 0
            hmr_factor = 0
            total_dynamics_length = assemble_mdp_file_inputs(segment, i, dt, 0, 0, hmr_factor, log_frequency, run_type="MIN")
        else:
            if i < len(equilibration_segment_parameters) - 1: # Equilibration segments
                dt = 0.002
                log_frequency = 5000 # 5000 * 0.002 = 10 ps / frame
                hmr_factor = 1
                total_steps = int(round(segment[0]/dt))
                total_dynamics_length = assemble_mdp_file_inputs(segment, i, dt, total_steps, total_dynamics_length, hmr_factor, log_frequency, run_type="EQ")    
            else: # Production segment
                dt = 0.004
                log_frequency = 25000 # 25000 * 0.004 = 100 ps / frame
                hmr_factor = 3
                total_steps = int(round(segment[0]/dt))
                total_dynamics_length = assemble_mdp_file_inputs(segment, i, dt, total_steps, total_dynamics_length, hmr_factor, log_frequency, run_type="PROD") * num_prod_runs
    return total_dynamics_length

# Generate input files
extract_inputs()
edit_topology_itp_files(small_molecules, solvent_groups, prot_fn)
total_dynamics_length = write_all_mdp_files(equilibration_segment_parameters, 10)
print(f"Total dynamics length: {total_dynamics_length} ps")
