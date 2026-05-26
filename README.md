# MD Analysis & Free Energy Landscape (FEL) Automation Suite

This suite is an all-in-one automation tool designed to analyze Molecular Dynamics (MD) simulation trajectories, calculate various structural parameters (RMSD, Rg, SASA, Hydrogen Bonds, DSSP, Dihedral, and Distance) using `cpptraj`, and automatically plot 1D, 2D, or 3D Free Energy Landscapes (FEL).

---

## Features

- **Multi-Dimensional FEL Plotting:**
  - **1D:** matplotlib-based line plot ($\Delta G$ vs Parameter).
  - **2D:** matplotlib-based contour plot (with custom energy level contour lines).
  - **3D:** Plotly-based interactive, rotatable, and zoomable 3D HTML scatter plot.
- **Discrete Binning for Integer Datasets:** Automatically centers bin boundaries to half-integers for variables like Hydrogen Bond counts or secondary structure counts to prevent gaps or NaN grid artifacts in contour maps.
- **Secondary Structure (DSSP) Analysis:** Tracks secondary structure timeline and generates a 2D timeline residue heatmap (`dssp_heatmap.png`) if pandas/seaborn are available.
- **Centralized Parameter File (`parameters.in`):** Manage all analysis masks, temperature settings, axis labels, and colormaps from a single file.
- **Performance Caching:** Queries `cpptraj` for the C-terminal residue index only once and caches it for all active mask auto-resolutions, speeding up initialization.
- **Clean Output Architecture:** Keeps intermediate `.in` input templates and `.log` logs in the output directory under dynamic names (`OUT_BASE`) and cleans them up automatically at the start of each run.

---

## nstallation

### 1. System Prerequisites
The backend uses **`cpptraj`** (provided by AmberTools) for calculations. Ensure `cpptraj` is installed and added to your system's `PATH`.

### 2. Python Dependencies
Install the required Python packages using:
```bash
pip install -r requirements.txt
```

---

## Configuration (`parameters.in`)

All defaults are loaded from the **`parameters.in`** file in the project directory. Key settings include:

| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `topology` | Path to the default topology file (prmtop/top) | `h2h3mono.top` |
| `trajectory` | Input file containing trajectory paths | `traj.in` |
| `output_dir` | Default directory to save output files and plots | `analysis` |
| `temperature` | Temperature in Kelvin used for FEL calculation | `310.0` |
| `energy_threshold_3d` | Maximum energy cutoff shown in 3D plots | `auto` (Dynamic estimation) |
| `pca_mask` | Atom mask used for PCA coordinate calculation | `@CA` |
| `pca_modes_visualize` | Number of PCA modes to calculate displacement trajectories for | `4` |
| `hbond_distance` | Cutoff distance for hydrogen bond acceptance (Å) | `3.0` |
| `hbond_angle` | Cutoff angle for hydrogen bond acceptance (degrees) | `135.0` |
| `secstruct_res` | Residue selection range for DSSP | `1-auto` (Auto C-terminal resolution) |
| `secstruct_type` | Secondary structure types to count (comma-separated) | `total` (All structures except coil) |

---

## Usage

### 1. Interactive Execution
Run the following command to start the analysis suite. You will be prompted to enter the output directory, dimensions, temperature, and target variables. Press **Enter** to accept the defaults from `parameters.in`:
```bash
python3 main.py
```

### 2. HPC Batch Job Execution with Slurm
For cluster queuing systems, you can submit jobs using the **Here Document (`<<EOF`)** approach to feed interactive inputs automatically in batch mode. Here is a sample Slurm submission script:

```bash
#!/bin/bash
#SBATCH -J md_fel_analysis
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -o analysis_%j.out
#SBATCH -e analysis_%j.err

cd $SLURM_SUBMIT_DIR

# Load required modules
module load python/3.9

# Run the program feeding inputs automatically
python3 main.py <<EOF
analysis                          # Output directory
2                                 # Dimensions (1, 2, or 3)
310.0                             # Temperature (Kelvin)
h2h3mono.top                      # Topology path
traj.in                           # Trajectory list path
RMSD                              # Variable 1
RG                                # Variable 2
EOF
```

---

## Project Structure

- **`main.py`:** Main orchestrator Python script that parses user inputs, loads configurations, and triggers analyses.
- **`plot_fel.py`:** Visualization utility that constructs histograms and computes free energies to generate 1D/2D/3D plots.
- **`calculate_*.sh`:** Specialized bash wrappers that generate input files and invoke `cpptraj` to calculate metrics (RMSD, Rg, SASA, Distance, Dihedral, HBOND, and DSSP).
- **`parameters.in`:** Configuration file.
- **`requirements.txt`:** List of Python library dependencies.
- **`.gitignore`:** Git pattern filter to prevent byte-code, logs, data files, and figures from cluttering the repository.

