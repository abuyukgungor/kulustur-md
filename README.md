# MD Analysis & Free Energy Landscape (FEL) Automation Suite

This suite is an all-in-one automation tool designed to analyze Molecular Dynamics (MD) simulation trajectories, calculate various structural parameters (RMSD, Rg, SASA, Hydrogen Bonds, DSSP, Dihedral, and Distance) using `cpptraj`, and automatically plot 1D, 2D, or 3D Free Energy Landscapes (FEL).

---

## Features

- **Multi-Dimensional FEL Plotting:**
  - **1D:** matplotlib-based line plot ($\Delta G$ vs Parameter).
  - **2D:** matplotlib-based contour plot (with custom energy level contour lines).
  - **3D:** Plotly-based interactive, rotatable, and zoomable 3D HTML scatter plot.
- **Non-interactive Batch Mode:** Supports execution via command-line arguments and flags, facilitating seamless integration with HPC Slurm schedulers.
- **Discrete Binning for Integer Datasets:** Automatically centers bin boundaries to half-integers for variables like Hydrogen Bond counts or secondary structure counts to prevent gaps or NaN grid artifacts in contour maps.
- **Secondary Structure (DSSP) Analysis:** Tracks secondary structure timeline and generates a 2D timeline residue heatmap (`dssp_heatmap.png` inside your output directory) if pandas is available.
- **PCA Mode Displacement:** Generates pseudo-trajectory files (`mode{i}.out` inside the output directory) along the major eigenvectors to visualize directional motions in PyMOL/VMD (controlled by `pca_modes_visualize`).
- **Centralized Parameter File (`parameters.in`):** Manage all analysis masks, temperature settings, axis labels, and colormaps from a single file.
- **Performance Caching:** Queries `cpptraj` for the C-terminal residue index only once and caches it for all active mask auto-resolutions, speeding up initialization.
- **Clean Output Architecture:** Keeps intermediate `.in` input templates and `.log` logs in the output directory under dynamic names (`OUT_BASE`) and cleans them up automatically at the start of each run.

---

## Installation

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
| `dimension` | Default FEL dimension (1, 2, or 3) | `2` |
| `temperature` | Temperature in Kelvin used for FEL calculation | `310.0` |
| `kB` | Boltzmann constant in kcal/mol·K | `0.001987` |
| `energy_threshold_3d` | Maximum energy cutoff shown in 3D plots | `auto` (Dynamic estimation) |
| `plot_format_3d` | Output format for 3D FEL plots (`html`, `png`, or `both`) | `html` |
| `pca_mask` | Atom mask used for PCA coordinate calculation | `@CA` |
| `pca_modes_visualize` | Number of PCA modes to calculate displacement trajectories for | `4` |
| `rmsd_reference` | Alignment reference structure (`average`, `first`, or path to PDB file) | `average` |
| `hbond_distance` | Cutoff distance for hydrogen bond acceptance (Å) | `3.0` |
| `hbond_angle` | Cutoff angle for hydrogen bond acceptance (degrees) | `135.0` |
| `secstruct_res` | Residue selection range for DSSP | `1-auto` (Auto C-terminal resolution) |
| `secstruct_type` | Secondary structure types to count (comma-separated) | `total` (All structures except coil) |
| `color_1d` | Hex color code for 1D line plotting | `#1f77b4` |
| `colormap_2d` | Color map for 2D plotting | `viridis` |
| `colormap_3d` | Color map for 3D Plotly plotting | `jet` |

---

## Usage

The suite supports both **Interactive Mode** (prompting for variables) and **Non-interactive Batch Mode** (passing options via CLI flags).

### 1. Command-Line Arguments & Flags

The orchestrator `main.py` accepts the following arguments:

| Flag | Long Flag | Description | Default Value |
| :--- | :--- | :--- | :--- |
| `-b` | `--batch` | Enable non-interactive batch mode (bypasses all prompts). | `False` |
| `-p` | `--params` | Path to the configuration file (e.g. `parameters.in`). | `parameters.in` |
| `-d` | `--dim` | Dimension of the FEL (1, 2, or 3). | Loaded from config or prompted |
| `-t` | `--top` | Topology file path (`.top`/`.prmtop`). | Loaded from config or prompted |
| `-x` | `--traj` | Trajectory list file path. | Loaded from config or prompted |
| `-o` | `--out` | Output directory name. | Loaded from config or prompted |
| | `--temp` | Simulation temperature in Kelvin. | Loaded from config or prompted |
| | `--var1` | First reaction coordinate/variable (e.g. `RMSD`, `RG`, etc.). | Prompted |
| | `--var2` | Second reaction coordinate (required for 2D/3D). | Prompted |
| | `--var3` | Third reaction coordinate (required for 3D). | Prompted |

### 2. Interactive Execution
Run the following command to start the analysis suite. You will be prompted to enter parameters. Any values specified via CLI flags will override the configuration file defaults and become the preselected default values for the prompts:
```bash
# Standard interactive mode (using default parameters.in):
python main.py

# Interactive mode overriding the config file and dimension:
python main.py -p local_parameters.in --dim 3
```

> [!NOTE]
> If no parameter file is specified using `-p` or `--params`, the program strictly falls back to the global `parameters.in` located in the script's directory. It does not automatically load a `parameters.in` in your current working directory.

### 3. Non-interactive Batch Execution (HPC/Slurm Recommended)
Specify the `-b`/`--batch` flag to bypass all prompts. If any required parameters are missing from both the CLI flags and the configuration file, the script will exit with an error.

```bash
# Run a 2D FEL analysis completely non-interactively
python main.py --batch --dim 2 --top h2h3mono.top --traj traj.in --var1 RMSD --var2 RG --out my_analysis
```

Here is a sample Slurm submission script utilizing the batch mode:

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

# Run program in batch mode (no user interaction required)
python main.py --batch --dim 2 --top h2h3mono.top --traj traj.in --var1 RMSD --var2 RG --out analysis
```
---


## Project Structure

- **`main.py`:** Main orchestrator Python script that parses user inputs, loads configurations, and triggers analyses.
- **`plot_fel.py`:** Visualization utility that constructs histograms and computes free energies to generate 1D/2D/3D plots.
- **`calculate_*.sh`:** Specialized bash wrappers that generate input files and invoke `cpptraj` to calculate metrics (RMSD, Rg, SASA, Distance, Dihedral, HBOND, and DSSP).
- **`parameters.in`:** Configuration file.
- **`requirements.txt`:** List of Python library dependencies.
