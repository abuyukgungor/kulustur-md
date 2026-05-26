import sys
import os
import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def load_parameters():
    defaults = {
        'kB': 0.001987,
        'bins_1d': 50,
        'bins_2d': 50,
        'bins_3d': 30,
        'color_1d': '#1f77b4',
        'colormap_2d': 'viridis',
        'colormap_3d': 'Jet',
        'plot_title': 'Free Energy Landscape ({dir_name})',
        'plot_format_3d': 'html'
    }
    
    config_file = None
    if "--config" in sys.argv:
        try:
            idx = sys.argv.index("--config")
            if idx + 1 < len(sys.argv):
                config_file = sys.argv[idx + 1]
        except Exception:
            pass
            
    if config_file is None:
        config_file = os.path.join(SCRIPT_DIR, 'parameters.in')
        
    config_file = os.path.abspath(config_file)
    
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
                        if key in ['kB', 'temperature']:
                            try:
                                defaults[key] = float(val)
                            except ValueError:
                                defaults[key] = val
                        elif key in ['bins_1d', 'bins_2d', 'bins_3d']:
                            try:
                                defaults[key] = int(val)
                            except ValueError:
                                defaults[key] = val
                        else:
                            defaults[key] = val
        except Exception:
            pass
    return defaults

def get_custom_label(var_label, defaults):
    key = f"label_{var_label.lower().strip()}"
    return defaults.get(key, var_label)

def resolve_colormap(name, default='viridis'):
    import matplotlib as mpl
    name_clean = str(name).strip()
    if hasattr(mpl, 'colormaps'):
        available_cmaps = list(mpl.colormaps.keys())
    else:
        available_cmaps = plt.colormaps()
    cmap_lower = name_clean.lower()
    for c in available_cmaps:
        if c.lower() == cmap_lower:
            return c
    return default

DEFAULTS = load_parameters()
kB = DEFAULTS['kB']

def load_data(filepath):
    # Reads cpptraj data file, skipping comments (# or @) and extracts
    # the value from the second column (index 1).
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append(float(parts[1]))
                except ValueError:
                    continue
    return np.array(data)

def plot_1d(file1, label1, output_image, T=310.0):
    print(f"--> Loading 1D data: {file1}...")
    x = load_data(file1)
    
    if len(x) == 0:
        raise ValueError("No readable numeric data found in 1D data file.")
        
    print("--> Calculating 1D histogram and Free Energy Landscape...")
    nbins = DEFAULTS.get('bins_1d', 50)
    
    # Auto-detect if data consists entirely of discrete integers
    is_discrete = np.all(np.abs(x - np.round(x)) < 1e-9)
    if is_discrete:
        # Align bins to half-integers to center each integer in a bin and prevent empty bin gaps
        bins_x = np.arange(np.min(x) - 0.5, np.max(x) + 1.5, 1.0)
        # Fallback if the range of integers is very large
        if len(bins_x) > 100:
            bins_x = nbins
    else:
        bins_x = nbins
        
    hist, xedges = np.histogram(x, bins=bins_x, density=True)
    
    hist = np.where(hist == 0, np.nan, hist)
    
    G = -(kB * T) * np.log(hist)
    G = G - np.nanmin(G)
    
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    
    # 1D Plotting
    dir_name = os.path.basename(os.path.dirname(os.path.abspath(output_image)))
    lbl1 = get_custom_label(label1, DEFAULTS)
    color1d = DEFAULTS.get('color_1d', '#1f77b4')
    plt.figure(figsize=(8, 6))
    plt.plot(xcenters, G, color=color1d, linewidth=2.5)
    
    plt.xlabel(lbl1, fontsize=12, fontweight='bold', labelpad=10)
    plt.ylabel('ΔG (kcal/mol)', fontsize=12, fontweight='bold', labelpad=10)
    
    # Resolve custom plot title
    title_tmpl = DEFAULTS.get('plot_title', 'Free Energy Landscape ({dir_name})')
    try:
        resolved_title = title_tmpl.format(dir_name=dir_name)
    except Exception:
        resolved_title = title_tmpl
    plt.title(resolved_title, fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=300)
    plt.close()
    print(f"--> Plot saved successfully: {output_image}")

def plot_2d(file1, file2, label1, label2, output_image, T=310.0):
    print(f"--> Loading 2D data: {file1} and {file2}...")
    data1 = load_data(file1)
    data2 = load_data(file2)
    
    if len(data1) == 0 or len(data2) == 0:
        raise ValueError("No readable numeric data found in 2D data files.")
        
    if len(data1) != len(data2):
        diff = abs(len(data1) - len(data2))
        if diff <= 10:
            min_len = min(len(data1), len(data2))
            print(f"--> [Warning] Data lengths differ by {diff} frames (Var1: {len(data1)}, Var2: {len(data2)}).")
            print(f"    Truncating both datasets to {min_len} frames to maintain synchronization.")
            data1 = data1[:min_len]
            data2 = data2[:min_len]
        else:
            raise ValueError(f"Data lengths differ significantly! Variable 1 has {len(data1)} frames, but Variable 2 has {len(data2)} frames. Difference is {diff} frames (max allowed for auto-truncation is 10). Frame synchronization is not guaranteed.")
        
    print("--> Calculating 2D histogram and Free Energy Landscape...")
    nbins = DEFAULTS.get('bins_2d', 50)
    
    # Auto-detect if data1 consists entirely of discrete integers
    is_discrete_1 = np.all(np.abs(data1 - np.round(data1)) < 1e-9)
    if is_discrete_1:
        bins_1 = np.arange(np.min(data1) - 0.5, np.max(data1) + 1.5, 1.0)
        if len(bins_1) > 100:
            bins_1 = nbins
    else:
        bins_1 = nbins
        
    # Auto-detect if data2 consists entirely of discrete integers
    is_discrete_2 = np.all(np.abs(data2 - np.round(data2)) < 1e-9)
    if is_discrete_2:
        bins_2 = np.arange(np.min(data2) - 0.5, np.max(data2) + 1.5, 1.0)
        if len(bins_2) > 100:
            bins_2 = nbins
    else:
        bins_2 = nbins
        
    hist, xedges, yedges = np.histogram2d(data1, data2, bins=[bins_1, bins_2], density=True)
    
    hist = np.where(hist == 0, np.nan, hist)
    
    G = -(kB * T) * np.log(hist)
    G = G - np.nanmin(G)
    
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(xcenters, ycenters)
    
    # 2D Plotting
    plt.figure(figsize=(9, 7.5))
    levels = np.linspace(0, np.nanmax(G), 50)
    cmap2d = resolve_colormap(DEFAULTS.get('colormap_2d', 'viridis'), 'viridis')
    cf = plt.contourf(X, Y, G.T, levels=levels, cmap=cmap2d)
    
    # Add energy contour lines (from 0.5 to 10 kcal/mol, step size 1.0)
    max_contour = min(10.0, np.nanmax(G))
    contour_levels = np.arange(0.5, max_contour, 1.0)
    if len(contour_levels) > 0:
        contours = plt.contour(X, Y, G.T, levels=contour_levels, colors='black', linewidths=0.5, alpha=0.3)
        plt.clabel(contours, inline=True, fontsize=8, fmt='%.1f kcal')
        
    cbar = plt.colorbar(cf, format='%.2f')
    cbar.set_label('ΔG (kcal/mol)', rotation=270, labelpad=20, fontsize=12, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    dir_name = os.path.basename(os.path.dirname(os.path.abspath(output_image)))
    lbl1 = get_custom_label(label1, DEFAULTS)
    lbl2 = get_custom_label(label2, DEFAULTS)
    plt.xlabel(lbl1, fontsize=12, fontweight='bold', labelpad=10)
    plt.ylabel(lbl2, fontsize=12, fontweight='bold', labelpad=10)
    
    # Resolve custom plot title
    title_tmpl = DEFAULTS.get('plot_title', 'Free Energy Landscape ({dir_name})')
    try:
        resolved_title = title_tmpl.format(dir_name=dir_name)
    except Exception:
        resolved_title = title_tmpl
    plt.title(resolved_title, fontsize=14, fontweight='bold', pad=15)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_image, dpi=300)
    plt.close()
    print(f"--> Plot saved successfully: {output_image}")

def plot_3d(file1, file2, file3, label1, label2, label3, output_path, T=310.0, threshold='auto', fmt='html'):
    print(f"--> Loading 3D data: {file1}, {file2}, {file3}...")
    x = load_data(file1)
    y = load_data(file2)
    z = load_data(file3)
    
    if len(x) == 0 or len(y) == 0 or len(z) == 0:
        raise ValueError("No readable numeric data found in 3D data files.")
        
    if len(x) != len(y) or len(x) != len(z):
        lens = [len(x), len(y), len(z)]
        min_len = min(lens)
        max_len = max(lens)
        diff = max_len - min_len
        if diff <= 10:
            print(f"--> [Warning] Data lengths differ by up to {diff} frames (Var1: {len(x)}, Var2: {len(y)}, Var3: {len(z)}).")
            print(f"    Truncating all datasets to {min_len} frames to maintain synchronization.")
            x = x[:min_len]
            y = y[:min_len]
            z = z[:min_len]
        else:
            raise ValueError(f"Data lengths differ significantly! Variable 1: {len(x)} frames, Variable 2: {len(y)} frames, Variable 3: {len(z)} frames. Difference is {diff} frames (max allowed for auto-truncation is 10). Frame synchronization is not guaranteed.")
        
    print("--> Calculating 3D histogram and Free Energy Landscape...")
    bins_xyz = DEFAULTS.get('bins_3d', 30)
    
    # Auto-detect if x consists entirely of discrete integers
    is_discrete_x = np.all(np.abs(x - np.round(x)) < 1e-9)
    if is_discrete_x:
        bins_x = np.arange(np.min(x) - 0.5, np.max(x) + 1.5, 1.0)
        if len(bins_x) > 100:
            bins_x = bins_xyz
    else:
        bins_x = bins_xyz
        
    # Auto-detect if y consists entirely of discrete integers
    is_discrete_y = np.all(np.abs(y - np.round(y)) < 1e-9)
    if is_discrete_y:
        bins_y = np.arange(np.min(y) - 0.5, np.max(y) + 1.5, 1.0)
        if len(bins_y) > 100:
            bins_y = bins_xyz
    else:
        bins_y = bins_xyz
        
    # Auto-detect if z consists entirely of discrete integers
    is_discrete_z = np.all(np.abs(z - np.round(z)) < 1e-9)
    if is_discrete_z:
        bins_z = np.arange(np.min(z) - 0.5, np.max(z) + 1.5, 1.0)
        if len(bins_z) > 100:
            bins_z = bins_xyz
    else:
        bins_z = bins_xyz
        
    hist, edges = np.histogramdd((x, y, z), bins=(bins_x, bins_y, bins_z), density=True)
    
    # Avoid ln(0) error by setting 0 values to NaN (matching 1D and 2D methods)
    hist = np.where(hist == 0, np.nan, hist)
    
    G = -(kB * T) * np.log(hist)
    G = G - np.nanmin(G)
    
    # Calculate threshold if set to auto
    if threshold == 'auto':
        try:
            detected_thresh = min(6.0, max(2.0, np.nanpercentile(G, 35)))
            print(f"--> Auto-detected 3D energy threshold: {detected_thresh:.2f} kcal/mol")
            threshold_val = detected_thresh
        except Exception:
            threshold_val = 3.0
    else:
        try:
            threshold_val = float(threshold)
        except ValueError:
            threshold_val = 3.0
            
    X, Y, Z = np.meshgrid(edges[0][:-1], edges[1][:-1], edges[2][:-1], indexing='ij')
    
    # Filtering lower energy states (e.g. G < threshold_val kcal/mol)
    mask = G < threshold_val
    X_plot = X[mask]
    Y_plot = Y[mask]
    Z_plot = Z[mask]
    G_plot = G[mask]
    
    base_path, ext = os.path.splitext(output_path)
    
    # 1. HTML generation using Plotly
    if fmt in ['html', 'both']:
        try:
            import plotly.graph_objects as go
        except ImportError:
            raise ImportError("plotly library is required for 3D HTML plotting! Please install it using: pip install plotly")
            
        print("--> Generating interactive 3D HTML plot (Plotly)...")
        fig = go.Figure(data=[go.Scatter3d(
            x=X_plot,
            y=Y_plot,
            z=Z_plot,
            mode='markers',
            marker=dict(
                size=5,
                color=G_plot,
                colorscale=DEFAULTS.get('colormap_3d', 'Jet'),
                opacity=0.8,
                colorbar=dict(title='ΔG (kcal/mol)')
            )
        )])
        
        dir_name = os.path.basename(os.path.dirname(os.path.abspath(output_path)))
        title_tmpl = DEFAULTS.get('plot_title', 'Free Energy Landscape ({dir_name})')
        try:
            resolved_title = title_tmpl.format(dir_name=dir_name)
        except Exception:
            resolved_title = title_tmpl
            
        fig.update_layout(
            title=dict(
                text=resolved_title,
                x=0.5,
                y=0.95,
                font=dict(size=16, color='black', family='Arial')
            ),
            scene=dict(
                xaxis_title=get_custom_label(label1, DEFAULTS),
                yaxis_title=get_custom_label(label2, DEFAULTS),
                zaxis_title=get_custom_label(label3, DEFAULTS)
            ),
            margin=dict(l=0, r=0, b=40, t=60)
        )
        
        html_out = base_path + ".html"
        fig.write_html(html_out)
        print(f"--> Interactive 3D plot saved successfully: {html_out}")
        
    # 2. PNG generation using Matplotlib 3D scatter
    if fmt in ['png', 'both']:
        print("--> Generating static 3D PNG plot (Matplotlib)...")
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        cmap_name = resolve_colormap(DEFAULTS.get('colormap_3d', 'jet'), 'jet')
        sc = ax.scatter(X_plot, Y_plot, Z_plot, c=G_plot, cmap=cmap_name, alpha=0.7, marker='o')
        
        cbar = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.1)
        cbar.set_label('ΔG (kcal/mol)', fontsize=11, fontweight='bold')
        
        ax.set_xlabel(get_custom_label(label1, DEFAULTS), fontsize=11, fontweight='bold', labelpad=10)
        ax.set_ylabel(get_custom_label(label2, DEFAULTS), fontsize=11, fontweight='bold', labelpad=10)
        ax.set_zlabel(get_custom_label(label3, DEFAULTS), fontsize=11, fontweight='bold', labelpad=10)
        
        dir_name = os.path.basename(os.path.dirname(os.path.abspath(output_path)))
        title_tmpl = DEFAULTS.get('plot_title', 'Free Energy Landscape ({dir_name})')
        try:
            resolved_title = title_tmpl.format(dir_name=dir_name)
        except Exception:
            resolved_title = title_tmpl
            
        plt.title(resolved_title, fontsize=13, fontweight='bold', pad=15)
        
        png_out = base_path + ".png"
        plt.savefig(png_out, dpi=300)
        plt.close()
        print(f"--> Static 3D plot saved successfully: {png_out}")

def main():
    args = sys.argv[1:]
    
    # Extract config if provided via --config to prevent it interfering with positional arguments
    if "--config" in args:
        try:
            idx = args.index("--config")
            if idx + 1 < len(args):
                del args[idx + 1]
            del args[idx]
        except Exception:
            pass
    
    # Extract temperature if provided via --temp
    temp = 310.0
    if "--temp" in args:
        try:
            idx = args.index("--temp")
            if idx + 1 < len(args):
                temp = float(args[idx + 1])
                del args[idx:idx+2]
        except Exception:
            pass
            
    # Extract threshold if provided via --threshold
    threshold = "auto"
    if "--threshold" in args:
        try:
            idx = args.index("--threshold")
            if idx + 1 < len(args):
                threshold = args[idx + 1]
                del args[idx:idx+2]
        except Exception:
            pass
    # Extract format if provided via --format
    fmt = "html"
    if "--format" in args:
        try:
            idx = args.index("--format")
            if idx + 1 < len(args):
                fmt = args[idx + 1].lower().strip()
                del args[idx:idx+2]
        except Exception:
            pass
            
    try:
        if len(args) == 3:
            # 1D: file1, label1, output
            plot_1d(args[0], args[1], args[2], T=temp)
        elif len(args) == 5:
            # 2D: file1, file2, label1, label2, output
            plot_2d(args[0], args[1], args[2], args[3], args[4], T=temp)
        elif len(args) == 7:
            # 3D: file1, file2, file3, label1, label2, label3, output
            plot_3d(args[0], args[1], args[2], args[3], args[4], args[5], args[6], T=temp, threshold=threshold, fmt=fmt)
        else:
            print("Error: Invalid arguments.")
            print("Usage:")
            print("  1D: python plot_fel.py <file1> <label1> <out>")
            print("  2D: python plot_fel.py <file1> <file2> <label1> <label2> <out>")
            print("  3D: python plot_fel.py <file1> <file2> <file3> <label1> <label2> <label3> <out> [--format html|png|both]")
            sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
