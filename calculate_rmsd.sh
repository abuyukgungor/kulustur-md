#!/bin/bash
# MD Analysis - RMSD Calculation against Average Structure
# Usage: ./calculate_rmsd.sh <topology> <traj_list> <out_file>

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
MASK=${4:-"@CA"}
REF=${5:-"average"}

# Check arguments
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file>"
    exit 1
fi

# Extract output directory from the out_file path
OUT_DIR=$(dirname "$OUT_FILE")
OUT_BASE=$(basename "$OUT_FILE" .dat)
AVG_FILE="$OUT_DIR/avg.pdb"

# Generate cpptraj trajin lines from trajectory list file
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

# Check if average structure exists in the output directory. If not, generate it.
if [ ! -f "$AVG_FILE" ]; then
    echo "--> Average structure ($AVG_FILE) not found. Generating average structure first..."
    cat <<EOF > "$OUT_DIR/step1_avg.in"
parm $TOPOLOGY
$TRAJIN_LINES
rms first $MASK
average $AVG_FILE
EOF
    cpptraj -i "$OUT_DIR/step1_avg.in" > "$OUT_DIR/cpptraj_avg.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: Average structure calculation failed! Details below:"
        cat "$OUT_DIR/cpptraj_avg.log"
        rm -f "$OUT_DIR/step1_avg.in"
        exit 1
    fi
    rm -f "$OUT_DIR/step1_avg.in"
fi

# Calculate RMSD using the specified reference
if [ "$REF" = "first" ]; then
    echo "--> Calculating RMSD against the first frame of the trajectory..."
    cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
rms first $MASK out $OUT_FILE
EOF
elif [ "$REF" = "average" ]; then
    echo "--> Calculating RMSD against the average structure ($AVG_FILE)..."
    cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
reference $AVG_FILE
$TRAJIN_LINES
rms reference $MASK out $OUT_FILE
EOF
else
    # It must be a path to a custom reference PDB file
    if [ ! -f "$REF" ]; then
        echo "ERROR: Custom reference PDB file '$REF' was not found!"
        exit 1
    fi
    echo "--> Calculating RMSD against custom reference structure ($REF)..."
    cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
reference $REF
$TRAJIN_LINES
rms reference $MASK out $OUT_FILE
EOF
fi

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: RMSD calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> RMSD calculation completed: $OUT_FILE"
