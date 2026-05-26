#!/bin/bash
# MD Analysis - Dihedral Calculation
# Usage: ./calculate_dihedral.sh <topology> <traj_list> <out_file> <mask> <mask> <mask> <mask>

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
MASK1=$4
MASK2=$5
MASK3=$6
MASK4=$7

# Input validation
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ] || [ -z "$MASK1" ] || [ -z "$MASK2" ] || [ -z "$MASK3" ] || [ -z "$MASK4" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file> <mask> <mask> <mask> <mask>"
    exit 1
fi

OUT_DIR=$(dirname "$OUT_FILE")
OUT_BASE=$(basename "$OUT_FILE" .dat)
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

echo "--> Calculating dihedral angle between $MASK1, $MASK2, $MASK3, $MASK4..."
cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
dihedral dihd $MASK1 $MASK2 $MASK3 $MASK4 out $OUT_FILE
EOF

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Dihedral calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> Dihedral calculation completed: $OUT_FILE"
