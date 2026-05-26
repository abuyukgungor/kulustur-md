#!/bin/bash
# MD Analysis - Secondary Structure (DSSP) Calculation
# Usage: ./calculate_secstruct.sh <topology> <traj_list> <out_file> <residues>

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
RESIDUES=$4

# Input validation
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ] || [ -z "$RESIDUES" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file> <residues_selection>"
    exit 1
fi

OUT_DIR=$(dirname "$OUT_FILE")
OUT_BASE=$(basename "$OUT_FILE" .dat)
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

# Ensure residues starts with colon for cpptraj mask selection
if [[ "$RESIDUES" != :* ]]; then
    MASK=":$RESIDUES"
else
    MASK="$RESIDUES"
fi

echo "--> Calculating secondary structure for residues $MASK..."
cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
secstruct ds1 $MASK out $OUT_DIR/dssp_time.gnu sumout $OUT_DIR/dssp_sum.dat
EOF

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Secondary structure calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> Secondary structure calculation completed. Output files generated in $OUT_DIR"
