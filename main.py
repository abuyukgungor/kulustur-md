import os
import sys
import subprocess

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def check_file(filepath, description):
    if not os.path.exists(filepath):
        print(f"Error: {description} '{filepath}' was not found!")
        sys.exit(1)

def get_last_residue_index(topology_path):
    try:
        # Run parminfo and resinfo to get both summary and detail
        result = subprocess.run(
            ["cpptraj", "-p", topology_path],
            input="parminfo\nresinfo\nexit\n",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        # Method 1: Look for "Total residues:" or "Residues:" in the output
        for line in result.stdout.splitlines():
            line_clean = line.strip().lower()
            # Handle "total residues: 60" or "residues: 60"
            for marker in ["total residues:", "residues:"]:
                if marker in line_clean:
                    parts = line_clean.split(marker)
                    if len(parts) > 1:
                        val_str = parts[1].strip().split()[0]
                        # Remove any non-digit chars (e.g. punctuation)
                        val_digits = ''.join(c for c in val_str if c.isdigit())
                        if val_digits:
                            return int(val_digits)
                            
        # Method 2: Fallback to parsing the resinfo table (first column is residue index)
        last_res = 1
        found_table_val = False
        for line in result.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith('#') or 'Residue' in line or '---' in line:
                continue
            parts = line.split()
            if parts and parts[0].isdigit():
                val = int(parts[0])
                if val > last_res:
                    last_res = val
                    found_table_val = True
                    
        if found_table_val:
            return last_res
            
        return None
    except Exception:
        return None


def parse_dssp_time_gnu(gnu_file, sec_type):
    from collections import defaultdict
    counts = defaultdict(int)
    max_frame = 1
    
    name_to_code = {
        'extended': 1,
        'para': 1,
        'bridge': 2,
        'anti': 2,
        '3-10': 3,
        'alpha': 4,
        'pi': 5,
        'turn': 6,
        'bend': 7
    }
    
    sec_type_clean = str(sec_type).lower().strip()
    allowed_codes = []
    
    if sec_type_clean == 'total':
        allowed_codes = [1, 2, 3, 4, 5, 6, 7]
    else:
        # Support comma-separated strings/digits (e.g. "extended,bridge" or "1,2")
        parts = sec_type_clean.split(',')
        for p in parts:
            p_strip = p.strip()
            if p_strip in name_to_code:
                allowed_codes.append(name_to_code[p_strip])
            elif p_strip.isdigit():
                allowed_codes.append(int(p_strip))
                
    if not allowed_codes:
        allowed_codes = [1, 2, 3, 4, 5, 6, 7]
        
    with open(gnu_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    frame = int(float(parts[0]))
                    if frame > max_frame:
                        max_frame = frame
                        
                    if len(parts) == 3:
                        # 3-column pm3d format: Frame  Residue  Code
                        code = int(float(parts[2]))
                        if code in allowed_codes:
                            counts[frame] += 1
                    else:
                        # Matrix format: Frame  Code1  Code2  ...
                        for p in parts[1:]:
                            code = int(float(p))
                            if code in allowed_codes:
                                counts[frame] += 1
                except ValueError:
                    continue
                    
    results = []
    for f in range(1, max_frame + 1):
        results.append((f, counts[f]))
    return results

def plot_dssp_heatmap(output_dir):
    try:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
        
        gnu_file = os.path.join(output_dir, "dssp_time.gnu")
        output_image = os.path.join(output_dir, "dssp_heatmap.png")
        
        if not os.path.exists(gnu_file):
            return
            
        print("--> Generating 2D DSSP Heatmap...")
        
        # Parse the gnu file manually to extract only numerical data lines and handle formats
        data_rows = []
        with open(gnu_file, 'r') as f:
            for line in f:
                line_strip = line.strip()
                if not line_strip or line_strip.startswith('#'):
                    continue
                parts = line_strip.split()
                try:
                    # Convert to float first, to filter out text commands (e.g. "set", "splot")
                    row_vals = [float(p) for p in parts]
                    if len(row_vals) > 0:
                        data_rows.append(row_vals)
                except ValueError:
                    # Skip gnuplot script commands or text labels
                    continue
                    
        if not data_rows:
            print("Warning: No numerical data found in DSSP gnu file.")
            return
            
        # Determine the format based on the number of columns in the first data row
        num_cols = len(data_rows[0])
        if num_cols == 3:
            # 3-column long format: Frame  Residue  Code
            df = pd.DataFrame(data_rows, columns=['Frame', 'Residue', 'Code'])
            df['Frame'] = df['Frame'].astype(int)
            df['Residue'] = df['Residue'].astype(int)
            df['Code'] = df['Code'].astype(int)
            # Pivot into index=Residue, columns=Frame, values=Code matrix
            df_matrix = df.pivot(index='Residue', columns='Frame', values='Code')
        else:
            # Wide matrix format: Frame  Res1  Res2  ...
            df = pd.DataFrame(data_rows)
            df.iloc[:, 0] = df.iloc[:, 0].astype(int)
            df.set_index(df.columns[0], inplace=True)
            df_matrix = df.T
            # Ensure the index represents residue numbers starting from 1
            df_matrix.index = range(1, len(df_matrix) + 1)
            
        cmap = ListedColormap(['white', 'cyan', 'blue', 'pink', 'red', 'magenta', 'green', 'yellow'])
        
        plt.figure(figsize=(15, 6))
        
        extent = [df_matrix.columns[0], df_matrix.columns[-1], df_matrix.index[0], df_matrix.index[-1]]
        im = plt.imshow(df_matrix, cmap=cmap, vmin=-0.5, vmax=7.5, aspect='auto', interpolation='nearest', extent=extent, origin='lower')
        
        # Invert y-axis so residue 1 is at the top and residue N is at the bottom
        plt.gca().invert_yaxis()
        
        # Add colorbar
        cbar = plt.colorbar(im, ticks=range(8))
        cbar.set_ticklabels(['Coil', 'Para-Beta', 'Anti-Beta', '3-10 Helix', 'Alpha Helix', 'Pi Helix', 'Turn', 'Bend'])
        cbar.ax.tick_params(labelsize=9)
        
        dir_name = os.path.basename(os.path.abspath(output_dir))
        plt.title(f"DSSP Secondary Structure Timeline ({dir_name})", fontsize=14, fontweight='bold', pad=15)
        plt.xlabel('Frame', fontsize=11, fontweight='bold')
        plt.ylabel('Residue No', fontsize=11, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_image, dpi=300)
        plt.close()
        print(f"--> DSSP Heatmap saved: {output_image}")
    except ImportError:
        print("Note: pandas or seaborn not installed. Skipping DSSP Heatmap generation.")
    except Exception as e:
        print(f"Warning: Failed to generate DSSP Heatmap: {e}")

def run_script(cmd_list):
    print(f"Running command: {' '.join(cmd_list)}")
    try:
        subprocess.run(cmd_list, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed (Exit code: {e.returncode})")
        sys.exit(1)

def process_variable(var_name, topology, traj_list, output_dir, config=None):
    if config is None:
        config = {}
        
    pca_mask = config.get('pca_mask', '@CA')
    pca_modes_visualize = config.get('pca_modes_visualize', 4)
    rmsd_mask = config.get('rmsd_mask', '@CA')
    rg_mask = config.get('rg_mask', '@CA')
    sasa_mask = config.get('sasa_mask', '!@H*')
    hbond_mask = config.get('hbond_mask', ':1-auto')
    hbond_distance = config.get('hbond_distance', 3.0)
    hbond_angle = config.get('hbond_angle', 135.0)
    secstruct_res = config.get('secstruct_res', '1-auto')
    secstruct_type = config.get('secstruct_type', 'total')
    dihedral_mask1 = config.get('dihedral_mask1', ':1-13@CA')
    dihedral_mask2 = config.get('dihedral_mask2', ':14@CA')
    dihedral_mask3 = config.get('dihedral_mask3', ':49@CA')
    dihedral_mask4 = config.get('dihedral_mask4', ':50-auto@CA')
    distance1_mask1 = config.get('distance1_mask1', ':1@CA')
    distance1_mask2 = config.get('distance1_mask2', 'auto')
    distance2_mask1 = config.get('distance2_mask1', ':1@CA')
    distance2_mask2 = config.get('distance2_mask2', 'auto')
    rmsd_reference = config.get('rmsd_reference', 'average')
    # 1. Check if the variable is an existing file (optionally with a column index, e.g. file.dat:2)
    path_part = var_name.strip()
    col_index = 1  # Default to column 1 (the second column, since 0 is usually frame index)
    
    if ":" in path_part:
        parts = path_part.rsplit(":", 1)
        if len(parts) == 2 and parts[1].isdigit():
            path_part = parts[0].strip()
            col_index = int(parts[1])
            
    if os.path.exists(path_part):
        filename = os.path.basename(path_part)
        label = os.path.splitext(filename)[0].upper()
        print(f"--> Found existing data file: {path_part} (using column {col_index})")
        
        # Save a clean copy with just (Frame, Value) in the output directory
        # so plot_fel.py can read it uniformly.
        cleaned_file = os.path.join(output_dir, f"ext_{label}_{col_index}.dat")
        try:
            with open(path_part, 'r') as f_in, open(cleaned_file, 'w') as f_out:
                frame_idx = 1
                for line in f_in:
                    line_strip = line.strip()
                    if not line_strip or line_strip.startswith('#') or line_strip.startswith('@'):
                        f_out.write(line)
                        continue
                    cols = line_strip.strip().split()
                    if len(cols) > col_index:
                        f_out.write(f"{frame_idx}\t{cols[col_index]}\n")
                        frame_idx += 1
            return cleaned_file, f"{label}_COL{col_index}"
        except Exception as e:
            print(f"Error parsing existing file '{path_part}': {e}")
            sys.exit(1)

    # 2. Built-in variables computed via cpptraj
    var_clean = var_name.upper().strip()
    
    # PCA Analysis (e.g., PC1, PC2)
    if var_clean.startswith("PC") and len(var_clean) > 2 and var_clean[2:].isdigit():
        pc_num = var_clean[2:]
        output_file = os.path.join(output_dir, f"pca{pc_num}.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_projection.sh"), topology, traj_list, pc_num, output_file, pca_mask, str(pca_modes_visualize)])
        return output_file, var_clean
        
    # RMSD Analysis
    elif var_clean == "RMSD":
        output_file = os.path.join(output_dir, "rmsd_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_rmsd.sh"), topology, traj_list, output_file, rmsd_mask, rmsd_reference])
        return output_file, var_clean
        
    # Radius of Gyration (Rg) Analysis
    elif var_clean in ["RG", "RADGYR"]:
        output_file = os.path.join(output_dir, "rg_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_rg.sh"), topology, traj_list, output_file, rg_mask])
        return output_file, "RG"
        
    # Solvent Accessible Surface Area (SASA) Analysis
    elif var_clean == "SASA":
        output_file = os.path.join(output_dir, "sasa_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_sasa.sh"), topology, traj_list, output_file, sasa_mask])
        return output_file, "SASA"
        
    # Secondary Structure (DSSP) Analysis
    elif var_clean in ["SECSTRUCT", "DSSP", "SECONDARYSTRUCTURE"]:
        output_file = os.path.join(output_dir, "secstruct_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_secstruct.sh"), topology, traj_list, output_file, secstruct_res])
        
        # Parse dssp_time.gnu in Python to extract counts
        gnu_file = os.path.join(output_dir, "dssp_time.gnu")
        check_file(gnu_file, "DSSP time data")
        parsed_data = parse_dssp_time_gnu(gnu_file, secstruct_type)
        
        # Write clean time series (Frame  Value) to secstruct_output.dat
        with open(output_file, 'w') as f_out:
            f_out.write(f"#Frame  {secstruct_type}_count\n")
            for frame, val in parsed_data:
                f_out.write(f"{frame}\t{val}\n")
                
        return output_file, "SECSTRUCT"
        
    # Hydrogen Bond (HBOND) Analysis
    elif var_clean in ["HBOND", "HBONDS", "HYDROGENBOND", "HYDROGENBONDS"]:
        output_file = os.path.join(output_dir, "hbond_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_hbond.sh"), topology, traj_list, output_file, hbond_mask, str(hbond_distance), str(hbond_angle)])
        return output_file, "HBOND"
        
    # Dihedral Angle Analysis
    elif var_clean in ["DIHEDRAL", "DIHD", "DIH", "DIHEDRALS"]:
        output_file = os.path.join(output_dir, "dihedral_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_dihedral.sh"), topology, traj_list, output_file, dihedral_mask1, dihedral_mask2, dihedral_mask3, dihedral_mask4])
        return output_file, "DIHEDRAL"
        
    # Distance 1 Analysis
    elif var_clean in ["DISTANCE", "DIST", "DISTANCE1", "DIST1"]:
        output_file = os.path.join(output_dir, "distance1_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_distance.sh"), topology, traj_list, output_file, distance1_mask1, distance1_mask2])
        return output_file, "DISTANCE1"
        
    # Distance 2 Analysis
    elif var_clean in ["DISTANCE2", "DIST2"]:
        output_file = os.path.join(output_dir, "distance2_output.dat")
        run_script(["bash", os.path.join(SCRIPT_DIR, "calculate_distance.sh"), topology, traj_list, output_file, distance2_mask1, distance2_mask2])
        return output_file, "DISTANCE2"
        
    else:
        print(f"Error: Unsupported variable name or file not found: '{var_name}'")
        print("Supported built-in variables: PC1, PC2, ... or RMSD, RG (RADGYR), SASA, HBOND, DIHEDRAL, DISTANCE1, DISTANCE2, SECSTRUCT (DSSP)")
        print("Alternatively, you can provide the path to an existing data file (e.g., dihedral.dat or file.dat:2)")
        sys.exit(1)

def safe_input(prompt):
    sys.stdout.write(prompt)
    sys.stdout.flush()
    try:
        raw_bytes = sys.stdin.buffer.readline()
        # Decode as utf-8 but ignore any invalid bytes (e.g. copy-paste artifacts)
        return raw_bytes.decode('utf-8', errors='ignore').strip()
    except Exception:
        # Fallback if stdin buffer is not accessible
        return input().strip()

def load_parameters():
    defaults = {
        'output_dir': 'analysis',
        'dimension': 2,
        'temperature': 310.0,
        'kB': 0.001987,
        'topology': 'h2h3mono.top',
        'trajectory': 'traj.in',
        'energy_threshold_3d': 'auto',
        'bins_1d': 50,
        'bins_2d': 50,
        'bins_3d': 30,
        'pca_mask': '@CA',
        'pca_modes_visualize': 4,
        'rmsd_mask': '@CA',
        'rg_mask': '@CA',
        'sasa_mask': '!@H*',
        'hbond_mask': ':1-auto',
        'hbond_distance': 3.0,
        'hbond_angle': 135.0,
        'secstruct_res': '1-auto',
        'secstruct_type': 'total',
        'dihedral_mask1': ':1-13@CA',
        'dihedral_mask2': ':14@CA',
        'dihedral_mask3': ':49@CA',
        'dihedral_mask4': ':50-auto@CA',
        'distance1_mask1': ':1@CA',
        'distance1_mask2': 'auto',
        'distance2_mask1': ':1@CA',
        'distance2_mask2': 'auto',
        'rmsd_reference': 'average'
    }
    config_file = os.path.join(SCRIPT_DIR, 'parameters.in')
    if os.path.exists(config_file):
        try:
            with open(config_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    if '=' in line:
                        key, val = line.split('=', 1)
                        key = key.strip()
                        val = val.strip()
                        # Convert to correct type with safe fallback
                        if key in ['temperature', 'kB', 'hbond_distance', 'hbond_angle']:
                            try:
                                defaults[key] = float(val)
                            except ValueError:
                                defaults[key] = val
                        elif key in ['dimension', 'bins_1d', 'bins_2d', 'bins_3d', 'pca_modes_visualize']:
                            try:
                                defaults[key] = int(val)
                            except ValueError:
                                defaults[key] = val
                        else:
                            defaults[key] = val
        except Exception as e:
            print(f"Warning: Failed to load parameters from {config_file}: {e}")
    return defaults

def main():
    print("--- MD Analysis & FEL Automation ---")
    
    # Load default configurations
    defaults = load_parameters()
    
    # Initialize flags defensively to prevent NameError if code flow changes
    has_secstruct = False
    
    # Get masks
    pca_mask = defaults.get('pca_mask', '@CA')
    pca_modes_visualize = defaults.get('pca_modes_visualize', 4)
    rmsd_mask = defaults.get('rmsd_mask', '@CA')
    rg_mask = defaults.get('rg_mask', '@CA')
    sasa_mask = defaults.get('sasa_mask', '!@H*')
    hbond_mask = defaults.get('hbond_mask', ':1-auto')
    hbond_distance = defaults.get('hbond_distance', 3.0)
    hbond_angle = defaults.get('hbond_angle', 135.0)
    dihedral_mask1 = defaults.get('dihedral_mask1', ':1-13@CA')
    dihedral_mask2 = defaults.get('dihedral_mask2', ':14@CA')
    dihedral_mask3 = defaults.get('dihedral_mask3', ':49@CA')
    dihedral_mask4 = defaults.get('dihedral_mask4', ':50-auto@CA')
    distance1_mask1 = defaults.get('distance1_mask1', ':1@CA')
    distance1_mask2 = defaults.get('distance1_mask2', 'auto')
    distance2_mask1 = defaults.get('distance2_mask1', ':1@CA')
    distance2_mask2 = defaults.get('distance2_mask2', 'auto')
    rmsd_reference = defaults.get('rmsd_reference', 'average')
    secstruct_res = defaults.get('secstruct_res', '1-auto')
    secstruct_type = defaults.get('secstruct_type', 'total')
    
    # Get the output directory name first
    d_dir = defaults.get('output_dir', 'wow')
    output_dir = safe_input(f"Output directory name [default: {d_dir}]: ").strip()
    if not output_dir:
        output_dir = d_dir
    os.makedirs(output_dir, exist_ok=True)
        
    # Get the dimension selection (1, 2, or 3)
    d_dim = defaults.get('dimension', 2)
    dim_str = safe_input(f"Choose dimension (1, 2, or 3) [default: {d_dim}]: ").strip()
    if not dim_str:
        dim = d_dim
    else:
        try:
            dim = int(dim_str)
            if dim not in [1, 2, 3]:
                raise ValueError
        except ValueError:
            print("Error: Dimension must be 1, 2, or 3!")
            sys.exit(1)
            
    # Get the temperature (default: 310 K)
    d_temp = defaults.get('temperature', 310.0)
    temp_str = safe_input(f"Temperature in Kelvin [default: {d_temp}]: ").strip()
    if not temp_str:
        temp = d_temp
    else:
        try:
            temp = float(temp_str)
            if temp <= 0:
                raise ValueError
        except ValueError:
            print("Error: Temperature must be a positive number!")
            sys.exit(1)

    # For 3D, get the energy threshold silently from parameters.in (default: auto)
    energy_threshold = "auto"
    if dim == 3:
        d_thresh = defaults.get('energy_threshold_3d', 'auto')
        energy_threshold = d_thresh
        
        # If it's a numeric override in parameters.in, validate it
        if energy_threshold != "auto":
            try:
                energy_threshold = float(energy_threshold)
                if energy_threshold <= 0:
                    raise ValueError
            except ValueError:
                print("Error: energy_threshold_3d in parameters.in must be 'auto' or a positive number!")
                sys.exit(1)
    
    # Get inputs from the user (Interactive Menu)
    d_top = defaults.get('topology', '')
    prompt_top = f" [default: {d_top}]" if d_top else ""
    topology = safe_input(f"Topology file path{prompt_top}: ").strip()
    if not topology:
        topology = d_top
    check_file(topology, "Topology file")
    
    d_traj = defaults.get('trajectory', '')
    prompt_traj = f" [default: {d_traj}]" if d_traj else ""
    traj_list = safe_input(f"Trajectory list file path{prompt_traj}: ").strip()
    if not traj_list:
        traj_list = d_traj
    check_file(traj_list, "Trajectory list file")
    
    print("Options: PC1, RMSD, RG, SASA, HBOND, DIHD, DIST, DSSP")

    if dim == 1:
        var1 = safe_input("Variable 1: ")
    elif dim == 2:
        var1 = safe_input("Variable 1: ")
        var2 = safe_input("Variable 2: ")
    else:
        var1 = safe_input("Variable 1: ")
        var2 = safe_input("Variable 2: ")
        var3 = safe_input("Variable 3: ")

    # Resolve distance1_mask2 / distance2_mask2 / hbond_mask / dihedral masks if calculations are requested
    vars_to_check = [var1]
    if dim >= 2:
        vars_to_check.append(var2)
    if dim == 3:
        vars_to_check.append(var3)
        
    vars_clean_check = [v.upper().strip() for v in vars_to_check]
    
    # Caching variables to query cpptraj at most once for last residue index
    last_res_idx = None
    queried_last_res = False
    
    # Check for Distance 1
    has_dist1 = any(v in ["DISTANCE", "DIST", "DISTANCE1", "DIST1"] for v in vars_clean_check)
    if has_dist1 and distance1_mask2 == 'auto':
        print("--> Auto-detecting last residue index from topology for Distance 1...")
        if not queried_last_res:
            last_res_idx = get_last_residue_index(topology)
            queried_last_res = True
        if last_res_idx:
            distance1_mask2 = f":{last_res_idx}@CA"
            print(f"--> Auto-detected C-terminal residue for Distance 1: {distance1_mask2}")
        else:
            print("Warning: Could not automatically detect the last residue index for Distance 1. Defaulting to :60@CA.")
            distance1_mask2 = ":60@CA"
            
    # Check for Distance 2
    has_dist2 = any(v in ["DISTANCE2", "DIST2"] for v in vars_clean_check)
    if has_dist2 and distance2_mask2 == 'auto':
        print("--> Auto-detecting last residue index from topology for Distance 2...")
        if not queried_last_res:
            last_res_idx = get_last_residue_index(topology)
            queried_last_res = True
        if last_res_idx:
            distance2_mask2 = f":{last_res_idx}@CA"
            print(f"--> Auto-detected C-terminal residue for Distance 2: {distance2_mask2}")
        else:
            print("Warning: Could not automatically detect the last residue index for Distance 2. Defaulting to :60@CA.")
            distance2_mask2 = ":60@CA"
            
    # Check for HBOND
    has_hbond = any(v in ["HBOND", "HBONDS", "HYDROGENBOND", "HYDROGENBONDS"] for v in vars_clean_check)
    if has_hbond and "auto" in hbond_mask:
        print("--> Auto-detecting last residue index for Hydrogen Bond mask...")
        if not queried_last_res:
            last_res_idx = get_last_residue_index(topology)
            queried_last_res = True
        if last_res_idx:
            hbond_mask = hbond_mask.replace("auto", str(last_res_idx))
            print(f"--> Auto-detected Hydrogen Bond mask: {hbond_mask}")
        else:
            print("Warning: Could not automatically detect the last residue index for HBOND. Defaulting to :1-60.")
            hbond_mask = ":1-60"
            
    # Check for SECSTRUCT
    has_secstruct = any(v in ["SECSTRUCT", "DSSP", "SECONDARYSTRUCTURE"] for v in vars_clean_check)
    if has_secstruct and "auto" in secstruct_res:
        print("--> Auto-detecting last residue index for Secondary Structure (DSSP) mask...")
        if not queried_last_res:
            last_res_idx = get_last_residue_index(topology)
            queried_last_res = True
        if last_res_idx:
            secstruct_res = secstruct_res.replace("auto", str(last_res_idx))
            print(f"--> Auto-detected Secondary Structure mask: {secstruct_res}")
        else:
            print("Warning: Could not automatically detect the last residue index for SECSTRUCT. Defaulting to 1-60.")
            secstruct_res = "1-60"
            
    # Check for DIHEDRAL
    has_dihedral = any(v in ["DIHEDRAL", "DIHD", "DIH", "DIHEDRALS"] for v in vars_clean_check)
    if has_dihedral:
        needs_auto = any("auto" in m for m in [dihedral_mask1, dihedral_mask2, dihedral_mask3, dihedral_mask4])
        if needs_auto:
            print("--> Auto-detecting last residue index for Dihedral masks...")
            if not queried_last_res:
                last_res_idx = get_last_residue_index(topology)
                queried_last_res = True
            if last_res_idx:
                if "auto" in dihedral_mask1: dihedral_mask1 = dihedral_mask1.replace("auto", str(last_res_idx))
                if "auto" in dihedral_mask2: dihedral_mask2 = dihedral_mask2.replace("auto", str(last_res_idx))
                if "auto" in dihedral_mask3: dihedral_mask3 = dihedral_mask3.replace("auto", str(last_res_idx))
                if "auto" in dihedral_mask4: dihedral_mask4 = dihedral_mask4.replace("auto", str(last_res_idx))
                print(f"--> Resolved Dihedral masks: {dihedral_mask1}, {dihedral_mask2}, {dihedral_mask3}, {dihedral_mask4}")
            else:
                print("Warning: Could not automatically detect the last residue index for DIHEDRAL. Using default/fallback 60.")
                if "auto" in dihedral_mask1: dihedral_mask1 = dihedral_mask1.replace("auto", "60")
                if "auto" in dihedral_mask2: dihedral_mask2 = dihedral_mask2.replace("auto", "60")
                if "auto" in dihedral_mask3: dihedral_mask3 = dihedral_mask3.replace("auto", "60")
                if "auto" in dihedral_mask4: dihedral_mask4 = dihedral_mask4.replace("auto", "60")
    
    # Clean up old intermediate files and logs from previous runs in the output directory
    print(f"\n--> Cleaning old intermediate analysis files and logs in '{output_dir}'...")
    files_to_clean = [
        "vect.out", "avg.pdb", "covar.out", 
        "cpptraj_avg.log", "cpptraj_cov.log", "cpptraj_proj.log", 
        "cpptraj_rmsd.log", "cpptraj_rg.log", "cpptraj_sasa.log", "cpptraj_hbond_output.log", "cpptraj_dihedral_output.log",
        "cpptraj_distance1.log", "cpptraj_distance2.log", "cpptraj_secstruct_output.log",
        "cpptraj_rmsd_output.log", "cpptraj_rg_output.log", "cpptraj_sasa_output.log",
        "cpptraj_distance1_output.log", "cpptraj_distance2_output.log",
        "step1_avg.in", "step1_cov.in", "step2.in", "rmsd.in", "rg.in", "sasa.in", "hbond_output.in", "dihedral_output.in",
        "distance1.in", "distance2.in", "secstruct_output.in",
        "rmsd_output.in", "rg_output.in", "sasa_output.in", "distance1_output.in", "distance2_output.in",
        "temp_proj.dat", "dssp_time.gnu", "dssp_sum.dat", "hbond_output_avg.dat"
    ]
    for i in range(1, 21):
        files_to_clean.append(f"mode{i}.out")
    for f in files_to_clean:
        path = os.path.join(output_dir, f)
        if os.path.exists(path):
            try:
                os.remove(path)
            except Exception:
                pass
                
    # Pack configuration for process_variable
    analysis_config = {
        'pca_mask': pca_mask,
        'pca_modes_visualize': pca_modes_visualize,
        'rmsd_mask': rmsd_mask,
        'rg_mask': rg_mask,
        'sasa_mask': sasa_mask,
        'hbond_mask': hbond_mask,
        'hbond_distance': hbond_distance,
        'hbond_angle': hbond_angle,
        'secstruct_res': secstruct_res,
        'secstruct_type': secstruct_type,
        'dihedral_mask1': dihedral_mask1,
        'dihedral_mask2': dihedral_mask2,
        'dihedral_mask3': dihedral_mask3,
        'dihedral_mask4': dihedral_mask4,
        'distance1_mask1': distance1_mask1,
        'distance1_mask2': distance1_mask2,
        'distance2_mask1': distance2_mask1,
        'distance2_mask2': distance2_mask2,
        'rmsd_reference': rmsd_reference
    }

    # Step 1: Run calculations
    print("\n>>> STEP 1: Running Analysis Scripts <<<")
    if dim == 1:
        data_file1, label1 = process_variable(var1, topology, traj_list, output_dir, analysis_config)
        check_file(data_file1, f"Analysis output for {var1}")
    elif dim == 2:
        data_file1, label1 = process_variable(var1, topology, traj_list, output_dir, analysis_config)
        data_file2, label2 = process_variable(var2, topology, traj_list, output_dir, analysis_config)
        check_file(data_file1, f"Analysis output for {var1}")
        check_file(data_file2, f"Analysis output for {var2}")
    else:
        data_file1, label1 = process_variable(var1, topology, traj_list, output_dir, analysis_config)
        data_file2, label2 = process_variable(var2, topology, traj_list, output_dir, analysis_config)
        data_file3, label3 = process_variable(var3, topology, traj_list, output_dir, analysis_config)
        check_file(data_file1, f"Analysis output for {var1}")
        check_file(data_file2, f"Analysis output for {var2}")
        check_file(data_file3, f"Analysis output for {var3}")
    
    # Step 2: Plot the Free Energy Landscape (FEL)
    print("\n>>> STEP 2: Plotting Free Energy Landscape (FEL) <<<")
    if dim == 1:
        output_image = os.path.join(output_dir, f"{var1.lower().strip()}_fel.png")
        run_script([sys.executable, os.path.join(SCRIPT_DIR, "plot_fel.py"), data_file1, label1, output_image, "--temp", str(temp)])
    elif dim == 2:
        output_image = os.path.join(output_dir, f"{var1.lower().strip()}_{var2.lower().strip()}_fel.png")
        run_script([sys.executable, os.path.join(SCRIPT_DIR, "plot_fel.py"), data_file1, data_file2, label1, label2, output_image, "--temp", str(temp)])
    else:
        output_image = os.path.join(output_dir, f"{var1.lower().strip()}_{var2.lower().strip()}_{var3.lower().strip()}_fel.html")
        run_script([sys.executable, os.path.join(SCRIPT_DIR, "plot_fel.py"), data_file1, data_file2, data_file3, label1, label2, label3, output_image, "--temp", str(temp), "--threshold", str(energy_threshold)])
    
    # Generate 2D Heatmap if secondary structure analysis was performed
    if has_secstruct:
        plot_dssp_heatmap(output_dir)
        
    print(f"\n[Success] Process completed! Your plot has been created: {output_image}")

if __name__ == "__main__":
    main()
