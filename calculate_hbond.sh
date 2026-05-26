#!/bin/bash
# MD Analysis - Hydrogen Bond Calculation
# Usage: ./calculate_hbond.sh <topology> <traj_list> <out_file> [mask]

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
MASK=${4:-":1-60"}
DISTANCE=${5:-"3.0"}
ANGLE=${6:-"135.0"}

# Input validation
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file> [mask] [distance] [angle]"
    exit 1
fi

OUT_DIR=$(dirname "$OUT_FILE")
OUT_BASE=$(basename "$OUT_FILE" .dat)
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

echo "--> Calculating internal hydrogen bonds for mask $MASK (dist: ${DISTANCE} A, angle: ${ANGLE} deg)..."
cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
hbond hb_int $MASK out $OUT_FILE dist $DISTANCE angle $ANGLE avgout ${OUT_DIR}/${OUT_BASE}_avg.dat
EOF

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Hydrogen bond calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> Hydrogen bond calculation completed: $OUT_FILE"
