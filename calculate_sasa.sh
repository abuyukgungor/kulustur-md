#!/bin/bash
# MD Analysis - Solvent Accessible Surface Area (SASA) Calculation
# Usage: ./calculate_sasa.sh <topology> <traj_list> <out_file> [mask]

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
MASK=${4:-"!@H*"}

# Check arguments
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file> [mask]"
    exit 1
fi

# Extract output directory from the out_file path
OUT_DIR=$(dirname "$OUT_FILE")
OUT_BASE=$(basename "$OUT_FILE" .dat)

# Generate cpptraj trajin lines from trajectory list file
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

# Calculate Solvent Accessible Surface Area using LCPO method
echo "--> Calculating Solvent Accessible Surface Area (SASA) for $MASK atoms..."
cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
surf sasa_out $MASK out $OUT_FILE
EOF

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: SASA calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> SASA calculation completed: $OUT_FILE"
