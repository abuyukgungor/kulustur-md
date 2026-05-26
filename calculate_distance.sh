#!/bin/bash
# MD Analysis - Distance Calculation
# Usage: ./calculate_distance.sh <topology> <traj_list> <out_file> [mask1] [mask2]

TOPOLOGY=$1
TRAJ_LIST=$2
OUT_FILE=$3
MASK1=${4:-":1@CA"}
MASK2=${5:-":60@CA"}

# Check arguments
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$OUT_FILE" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <output_file> [mask1] [mask2]"
    exit 1
fi

# Extract output directory from the out_file path
OUT_DIR=$(dirname "$OUT_FILE")

# Generate cpptraj trajin lines from trajectory list file
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

# Extract filename base of OUT_FILE to use for .in and .log filenames
OUT_BASE=$(basename "$OUT_FILE" .dat)

# Calculate distance
echo "--> Calculating distance between $MASK1 and $MASK2..."
cat <<EOF > "$OUT_DIR/${OUT_BASE}.in"
parm $TOPOLOGY
$TRAJIN_LINES
distance dist_val $MASK1 $MASK2 out $OUT_FILE
EOF

cpptraj -i "$OUT_DIR/${OUT_BASE}.in" > "$OUT_DIR/cpptraj_${OUT_BASE}.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Distance calculation failed! Details below:"
    cat "$OUT_DIR/cpptraj_${OUT_BASE}.log"
    rm -f "$OUT_DIR/${OUT_BASE}.in"
    exit 1
fi
rm -f "$OUT_DIR/${OUT_BASE}.in"

echo "--> Distance calculation completed: $OUT_FILE"
