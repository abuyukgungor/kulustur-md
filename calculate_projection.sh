#!/bin/bash
# MD Analysis - Principal Component Projection (PCA Projection)
# Usage: ./calculate_projection.sh <topology> <traj_list> <pc_num> <out_file>

TOPOLOGY=$1
TRAJ_LIST=$2
PC_NUM=$3
OUT_FILE=$4
MASK=${5:-"@CA"}
MODES_VISUALIZE=${6:-"4"}

# Check arguments
if [ -z "$TOPOLOGY" ] || [ -z "$TRAJ_LIST" ] || [ -z "$PC_NUM" ] || [ -z "$OUT_FILE" ]; then
    echo "Error: Missing arguments."
    echo "Usage: $0 <topology_file> <trajectory_list_file> <pc_number> <output_file>"
    exit 1
fi

# Extract output directory from the out_file path
OUT_DIR=$(dirname "$OUT_FILE")
VECT_FILE="$OUT_DIR/vect.out"
AVG_FILE="$OUT_DIR/avg.pdb"
COVAR_FILE="$OUT_DIR/covar.out"

# Generate cpptraj trajin lines from trajectory list file
TRAJIN_LINES=$(awk 'NF {print "trajin "$1}' "$TRAJ_LIST")

# Step 1: Average structure calculation and Covariance Matrix / Eigenvectors diagonalization
# (Run only if they have not been generated yet)
if [ ! -f "$VECT_FILE" ] || [ ! -f "$AVG_FILE" ]; then
    echo "--> Step 1a: Calculating average structure ($AVG_FILE)..."
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

    echo "--> Step 1b: Calculating covariance matrix and eigenvectors ($VECT_FILE)..."
    DISPL_LINES=""
    for ((i=1; i<=MODES_VISUALIZE; i++)); do
        DISPL_LINES="${DISPL_LINES}analyze modes displ stack vector beg $i end $i out $OUT_DIR/mode${i}.out
"
    done
    cat <<EOF > "$OUT_DIR/step1_cov.in"
parm $TOPOLOGY
reference $AVG_FILE
$TRAJIN_LINES
rms reference $MASK
matrix mwcovar name mwcov $MASK out $COVAR_FILE
run
analyze matrix mwcov name vector out $VECT_FILE vecs 10
$DISPL_LINES
EOF
    cpptraj -i "$OUT_DIR/step1_cov.in" > "$OUT_DIR/cpptraj_cov.log" 2>&1
    if [ $? -ne 0 ]; then
        echo "ERROR: Covariance matrix calculation failed! Details below:"
        cat "$OUT_DIR/cpptraj_cov.log"
        rm -f "$OUT_DIR/step1_cov.in"
        exit 1
    fi
    rm -f "$OUT_DIR/step1_cov.in"
fi

# Validate Step 1 outputs
if [ ! -f "$VECT_FILE" ] || [ ! -f "$AVG_FILE" ]; then
    echo "Error: PCA preparation files ($VECT_FILE or $AVG_FILE) were not created."
    exit 1
fi

# Step 2: Projection onto the specified PC index
echo "--> Step 2: Projecting trajectory coordinates onto PC$PC_NUM..."
cat <<EOF > "$OUT_DIR/step2.in"
parm $TOPOLOGY
reference $AVG_FILE
readdata $VECT_FILE name map
$TRAJIN_LINES
rms reference $MASK
projection proj modes map beg 1 end $PC_NUM $MASK out $OUT_DIR/temp_proj.dat
EOF

cpptraj -i "$OUT_DIR/step2.in" > "$OUT_DIR/cpptraj_proj.log" 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Coordinate projection failed! Details below:"
    cat "$OUT_DIR/cpptraj_proj.log"
    rm -f "$OUT_DIR/step2.in"
    exit 1
fi
rm -f "$OUT_DIR/step2.in"

# Filter the desired PC mode from the cpptraj output (Column 1: Frame, Column PC_NUM+1: Value)
awk -v col=$((PC_NUM + 1)) -v pc="$PC_NUM" '/^#/ {if ($1 == "#Frame") {print "#Frame", "PC"pc} else {print "#"}; next} {print $1, $col}' "$OUT_DIR/temp_proj.dat" > "$OUT_FILE"
rm -f "$OUT_DIR/temp_proj.dat"

echo "--> PC$PC_NUM projection calculation completed: $OUT_FILE"
