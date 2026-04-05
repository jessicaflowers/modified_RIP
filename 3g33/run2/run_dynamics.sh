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
# Append _PERT to base so outputs (e.g. seg_11_NPT_PROD_PERT_01) are
# distinct from the perturbed input (seg_11_NPT_PROD_01_PERT.gro)
ISTEP_BASE="${ISTEP_BASE}_PERT"

echo "Running production MDP: ${MDP_PATH}"
echo "Topology: ${TOP}"

for (( r=1; r<=PROD_REPEAT; r++ )); do
    iter=$(printf "%02d" "$r")
    ISTEP="${ISTEP_BASE}_${iter}"

    # Skip iterations that already finished
    if [ -f "${ISTEP}.log" ] && grep -q "Finished mdrun" "${ISTEP}.log" 2>/dev/null; then
        echo "Iteration ${iter} already finished. Skipping."
        continue
    fi

    # If a checkpoint exists from a killed run, resume directly with mdrun
    if [ -f "${ISTEP}.cpt" ] && [ -f "${ISTEP}.tpr" ]; then
        echo "Resuming ${ISTEP} from checkpoint."
        gmx mdrun -v -deffnm "${ISTEP}" -cpi "${ISTEP}".cpt -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off
        continue
    fi

    # Determine the previous step's structure
    if [ "$r" -eq 1 ]; then
        PSTEP="${MIN_OUTPUT}"
    else
        prev_iter=$(printf "%02d" $((r - 1)))
        PSTEP="${ISTEP_BASE}_${prev_iter}"
    fi

    echo "----Running ${ISTEP} (Previous step: ${PSTEP})----"

    CPT_ARG=""
    [ -f "${PSTEP}.cpt" ] && CPT_ARG="-t ${PSTEP}.cpt"

    gmx grompp -f "${MDP_PATH}" -o "${ISTEP}".tpr ${CPT_ARG} -c "${PSTEP}".gro -r "${PSTEP}".gro -p "${TOP}" -n "${NDX}" -maxwarn 10
    gmx mdrun -v -deffnm "${ISTEP}" -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off
done

echo "Production run finished. Check logs for outputs."
