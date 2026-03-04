#!/bin/bash
# Run equilibration and production steps
# The first step is NVT after minimization: seg_02_NVT_EQ.mdp
source /wynton/home/rotation/jessicaflowers/software/gromacs/2025.4_CUDA/bin/GMXRC
# 1. Environment setup

# 2. Configuration
# Original static paths (kept commented for reference)
# TOP="../input/topol.top"
# NDX="../input/index.ndx"
# MIN_OUTPUT="seg_01_MIN"
# ALL_MDPS=$(ls ../mdp/seg_*.mdp | sort -V )

# Threading
NTHREAD=8
NTOMP=8
PROD_REPEAT=10

# Determine input/top/mdp locations relative to current directory so the
# script can be executed either from the parent run1 directory or from inside
# a per-residue pulse directory (e.g. pulse_res_3). We search a few ancestor
# levels for an input/ and mdp/ directory.
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
if [ -z "${TOP}" ]; then
    echo "Warning: could not find topol.top in ancestor ../input or ../../input; using ../input/topol.top as fallback"
    TOP="../input/topol.top"
    NDX="../input/index.ndx"
fi

for up in . .. ../.. ../../..; do
    if [ -d "${up}/mdp" ]; then
        MDP_DIR="${up}/mdp"
        break
    fi
done
if [ -z "${MDP_DIR}" ]; then
    echo "Warning: could not find mdp/ directory in expected ancestors; using ../mdp as fallback"
    MDP_DIR="../mdp"
fi

# choose the minimization/starting .gro available in the working directory
# if [ -f "seg_01_PERT.gro" ]; then
#     MIN_OUTPUT="seg_01_PERT"
# elif [ -f "seg_01_MIN.gro" ]; then
#     MIN_OUTPUT="seg_01_MIN"
# else
#     # fallback to original name; this will probably fail unless the file exists
#     MIN_OUTPUT="seg_01_MIN"
# fi
MIN_OUTPUT="seg_01_PERT" # the stuff commented out above is bad. it should only use the pert file

ALL_MDPS=$(ls "${MDP_DIR}"/seg_*.mdp | sort -V )
TOTAL_STEPS=$(echo "${ALL_MDPS}" | wc -l)

echo "Found $TOTAL_STEPS MDPS in ${MDP_DIR}"

# 3. Resume
# Determine the last gro/log in this working directory (ignore MIN)
LAST_GRO=$(ls seg_*.gro 2>/dev/null | grep -v "MIN" | sort -V | tail -n1 || true)
LAST_LOG=$(ls seg_*.log 2>/dev/null | grep -v "MIN" | sort -V | tail -n1 || true)
echo "last gro: ${LAST_GRO}"
echo "last log: ${LAST_LOG}"

# If there is no previous run output, start from the first equilibration (step 2)
if [ -z "$LAST_GRO" ]; then
    echo "Starting from first equilibration step."
    START_CNT=2
    CONTINUE_CPT=false
else
    # If there is a log file from a previous step, use it to determine whether
    # the last step finished. If there is no log file (e.g. only seg_01_PERT.gro
    # exists), treat this as a fresh start and begin at step 2.
    if [ -z "$LAST_LOG" ]; then
        echo "Found GRO (${LAST_GRO}) but no log files; starting from first equilibration step."
        START_CNT=2
        CONTINUE_CPT=false
    else
        # extract the numeric ID from the last log (seg_05_NPT_EQ.log -> 05)
        LAST_ID=$(echo "$LAST_LOG" | awk -F '[._]' '{print $2}')
        LAST_BASE=$(basename "$LAST_LOG" .log)
        echo "last log: ${LAST_LOG}"
        if grep -q "Finished mdrun" "$LAST_LOG" 2>/dev/null; then
            START_CNT=$((10#$LAST_ID + 1))
            CONTINUE_CPT=false
            echo "Last step $LAST_ID finished. Moving to ${START_CNT}."
        else
            START_CNT=$((10#$LAST_ID))
            CONTINUE_CPT=true
            echo "Last step $LAST_ID unfinished. Continuing..."
        fi
    fi
fi

# 4. Run gmx
for (( cnt=$START_CNT; cnt<=$TOTAL_STEPS; cnt++)); do
    curr_id=$(printf "%02d" "$cnt")

    # grab MDP and base name for this step
    # previously: MDP_PATH=$(ls ../mdp/seg_${curr_id}*.mdp | head -n1)
    MDP_PATH=$(ls "${MDP_DIR}"/seg_${curr_id}*.mdp 2>/dev/null | head -n1)
    if [ -z "${MDP_PATH}" ]; then
        echo "Warning: no MDP file found for step ${curr_id} in ${MDP_DIR} (pattern seg_${curr_id}*.mdp). Skipping."
        continue
    fi
    ISTEP_BASE=$(basename "$MDP_PATH" .mdp)

    if echo "$ISTEP_BASE" | grep -q "PROD"; then
        # --- PRODUCTION SUBLOOP ---    
        for (( r=1; r<=$PROD_REPEAT; r++)); do
            iter=$(printf "%02d" "$r")
            ISTEP="${ISTEP_BASE}_${iter}"

            if [ -f "${ISTEP}.gro" ] && grep -q "Finished mdrun" "${ISTEP}.log"; then
                echo "Iteration $iter finished. Moving to next iteration."
                continue
            fi

            # Determine previous step
            if [ "$r" -eq 1 ]; then # first iteration of production run
                PREV_ID=$(printf "%02d" $((10#$cnt - 1)))
                PSTEP=$(ls seg_${PREV_ID}*.gro 2>/dev/null | sort -V | tail -n1 | sed 's/\.gro//')
            else
                PREV_ITER=$(printf "%02d" $((10#$r - 1)))
                PSTEP=${ISTEP_BASE}_${PREV_ITER}
            fi

            echo "----Running $ISTEP (Previous step: $PSTEP)----"
            if [ "$CONTINUE_CPT" = true ] && [ -f "${ISTEP}.cpt" ]; then
                gmx mdrun -v -deffnm "${ISTEP}" -cpi "${ISTEP}".cpt -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off
            else
                CPT_ARG=""
                [ -f "${PSTEP}.cpt" ] && CPT_ARG="-t ${PSTEP}.cpt"
                gmx grompp -f "${MDP_PATH}" -o "${ISTEP}".tpr ${CPT_ARG} -c "${PSTEP}".gro -r "${PSTEP}".gro -p "${TOP}" -n "${NDX}" -maxwarn 10
                gmx mdrun -v -deffnm "${ISTEP}" -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off
            fi
        done
        # --- END PRODUCTION SUBLOOP ---

    else
        # --- EQUILIBRATION  ---
        ISTEP="${ISTEP_BASE}"
        if [ "$CONTINUE_CPT" = true ]; then
            echo "Continuing from checkpoint for step $curr_id"
            gmx mdrun -v -deffnm "${ISTEP}" -cpi "${ISTEP}".cpt -nt "${NTHREAD}" -ntomp "${NTOMP}" -ntmpi 1 -nb gpu -bonded gpu -pin off
            CONTINUE_CPT=false
        else
            echo "Starting new run for step $curr_id"
            if [ "$cnt" -eq 2 ]; then
                PSTEP="${MIN_OUTPUT}"
            else
                PREV_ID=$(printf "%02d" $((10#$cnt - 1)))
                PSTEP=$(ls seg_${PREV_ID}*.gro 2>/dev/null | sort -V | tail -n1 | sed 's/\.gro//')
            fi

            echo "----Running $ISTEP (Previous step: $PSTEP)----"

            CPT_ARG=""
            [ -f "${PSTEP}.cpt" ] && CPT_ARG="-t ${PSTEP}.cpt"

            gmx grompp \
                -f "${MDP_PATH}" \
                -o "${ISTEP}".tpr \
                ${CPT_ARG} \
                -c "${PSTEP}".gro \
                -r "${PSTEP}".gro \
                -p "${TOP}" \
                -n "${NDX}" \
                -maxwarn 10
            gmx mdrun \
                -v -deffnm "${ISTEP}" \
                -nt "${NTHREAD}" \
                -ntomp "${NTOMP}" \
                -ntmpi 1 \
                -nb gpu \
                -bonded gpu \
                -pin off
        fi
    fi
done

