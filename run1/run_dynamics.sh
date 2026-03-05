#!/bin/bash
# Run production steps

NTHREAD=8
NTOMP=8
PROD_REPEAT=10
source /wynton/home/rotation/jessicaflowers/software/gromacs/2025.4_CUDA/bin/GMXRC


# Locate topology and mdp directories (search upwards)
TOP=""
NDX=""
MDP_DIR=""
for up in . .. ../.. ../../..; do
    if [ -f "${up}/input/topol.top" ]; then
        TOP="${up}/input/topol.top"
        NDX="${up}/input/index.ndx"
        break
    fi
done
[ -z "${TOP}" ] && TOP="../input/topol.top" && NDX="../input/index.ndx"
for up in . .. ../.. ../../..; do
    if [ -d "${up}/mdp" ]; then
        MDP_DIR="${up}/mdp"
        break
    fi
done
[ -z "${MDP_DIR}" ] && MDP_DIR="../mdp"

# Default starting structure base name (without .gro)
MIN_OUTPUT="seg_11_NPT_PROD_01_PERT"

# Pick production MDP: argument or first seg_*PROD*.mdp in mdp dir
ARG_MDP="$1"
MDP_PATH=""
if [ -n "${ARG_MDP}" ]; then
    if [ -f "${ARG_MDP}" ]; then
        MDP_PATH="${ARG_MDP}"
    elif [ -f "${MDP_DIR}/${ARG_MDP}" ]; then
        MDP_PATH="${MDP_DIR}/${ARG_MDP}"
    elif [ -f "${MDP_DIR}/${ARG_MDP}.mdp" ]; then
        MDP_PATH="${MDP_DIR}/${ARG_MDP}.mdp"
    fi
fi
if [ -z "${MDP_PATH}" ]; then
    MDP_PATH=$(ls "${MDP_DIR}"/seg_*PROD*.mdp 2>/dev/null | head -n1 || true)
fi
if [ -z "${MDP_PATH}" ]; then
    echo "Error: no production MDP found in ${MDP_DIR} and none provided as arg."
    exit 1
fi

# Optional overrides
[ -n "$2" ] && MIN_OUTPUT="$2"
[ -n "$3" ] && TOP="$3"

if [ ! -f "${MIN_OUTPUT}.gro" ]; then
    echo "Error: starting structure ${MIN_OUTPUT}.gro not found in working directory."
    exit 1
fi

ISTEP_BASE=$(basename "${MDP_PATH}" .mdp)
ISTEP="${ISTEP_BASE}_01"

echo "Running production MDP: ${MDP_PATH}"
echo "Start structure: ${MIN_OUTPUT}.gro -> output deffnm: ${ISTEP}"
echo "Topology: ${TOP}"

# If a checkpoint for this deffnm exists, use it, otherwise start from MIN_OUTPUT
CPT_ARG=""
[ -f "${ISTEP}.cpt" ] && CPT_ARG="-t ${ISTEP}.cpt"

gmx grompp -f "${MDP_PATH}" -o "${ISTEP}".tpr ${CPT_ARG} -c "${MIN_OUTPUT}".gro -r "${MIN_OUTPUT}".gro -p "${TOP}" -n "${NDX}" -maxwarn 10
gmx mdrun -v -deffnm "${ISTEP}" -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off

echo "Production run finished (or failed). Check ${ISTEP}.log and ${ISTEP}.trr/.cpt for outputs."

exit 0
if [ -z "${PROD_MDPS}" ]; then

    echo "No production MDP files (seg_*PROD*.mdp) found in ${MDP_DIR}; nothing to run."
