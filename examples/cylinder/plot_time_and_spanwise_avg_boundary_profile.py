# A simple script to calculate time-weighted and spanwise-integrated spatial averages from Neko's boundary_profile simcomp data.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


data_file = 'boundary_profile_cylinder_wall_data.csv'
mesh_file = 'boundary_profile_cylinder_wall_mesh.csv'
target_variable = 'tau_mag'

# Time Averaging
# Set t_end to None to average until the last snapshot.
t_start = 5.0
t_end = None


# ==========================================
def main():
    print("Loading mesh data...")
    df_mesh = pd.read_csv(mesh_file)

    print(f"Loading boundary_profile data to extract '{target_variable}'...")
    data_lines = []
    with open(data_file, 'r') as f:
        columns = f.readline().strip().split(',')

        n_points = int(f.readline().strip())
        print(f"Detected {n_points} spatial points per time step.")

        current_time = 0.0
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith('t='):
                current_time = float(line.split('=')[1].strip())
            else:
                row = line.split(',')
                row.append(current_time)
                data_lines.append(row)

    columns_with_time = columns + ['time']
    df_data = pd.DataFrame(data_lines, columns=columns_with_time).astype(float)

    if target_variable not in df_data.columns:
        raise ValueError(f"Error: '{target_variable}' not found. Available: {columns_with_time}")

    total_snapshots = df_data['time'].unique()
    print(f"Read {len(total_snapshots)} total time snapshots (t={total_snapshots.min()} to t={total_snapshots.max()}).")

    if t_start is not None:
        df_data = df_data[df_data['time'] >= t_start]

    if t_end is not None:
        df_data = df_data[df_data['time'] <= t_end]

    if df_data.empty:
        raise ValueError(f"Error: No data found within the specified time window (t_start={t_start}, t_end={t_end}).")

    filtered_snapshots = df_data['time'].unique()
    print(f"---> Using {len(filtered_snapshots)} snapshots for time averaging (t={filtered_snapshots.min()} to t={filtered_snapshots.max()}).")

    # Calculate Time-Weighted Average
    print("Calculating time-weighted average...")
    df_time_avg = df_data.groupby('gll_id').apply(
        lambda group: time_weighted_average(group, target_variable)
    ).reset_index(name=target_variable)

    df_merged = pd.merge(df_mesh, df_time_avg, on='gll_id')
    df_merged['theta'] = np.degrees(np.arctan2(df_merged['y_ref'], df_merged['x_ref']))
    df_merged['theta_r'] = df_merged['theta'].round(6)

    # Spanwise-Average (Average over Z)
    print("Calculating spatial-weighted spanwise average...")
    df_final = df_merged.groupby('theta_r').apply(
        lambda group: spanwise_weighted_average(group, target_variable)
    ).reset_index(name=target_variable)

    df_final = df_final.sort_values(by='theta_r')

    # Plot
    print("Plotting...")
    plt.figure(figsize=(10, 6))
    plt.plot(df_final['theta_r'], df_final[target_variable], '-', color='blue', linewidth=2, label=f"{target_variable}_mean")

    window_title = f" (t={filtered_snapshots.min()} to {filtered_snapshots.max()})"
    plt.title(f'Time- and Spanwise-Averaged {target_variable} vs. Theta\n{window_title}', fontsize=14)
    plt.xlabel(r'Angle $\theta$ (degrees)', fontsize=12)
    plt.ylabel(f'Mean {target_variable}', fontsize=12)
    plt.xlim(-180, 180)
    plt.xticks(np.arange(-180, 181, 45))
    plt.grid(True, linestyle='--')
    plt.legend()
    plt.tight_layout()
    plt.show()

# ==========================================
def time_weighted_average(group, target_var):
    """Calculates the time integral of a variable to handle variable time steps."""
    group = group.sort_values('time')

    t = group['time'].values
    y = group[target_var].values

    if len(t) > 1:
        total_time = t[-1] - t[0]
        if total_time > 0:
            return np.trapz(y, x=t) / total_time

    return np.mean(y)

def spanwise_weighted_average(group, target_var):
    """Handles non-uniform GLL points by spatial integration."""
    # Sort along the spanwise direction
    group = group.sort_values('z_ref')

    # Merge duplicate Z-coordinates
    group = group.groupby('z_ref', as_index=False)[target_var].mean()

    z = group['z_ref'].values
    y = group[target_var].values

    # Integrate over space
    if len(z) > 1:
        span_length = z[-1] - z[0]
        if span_length > 0:
            return np.trapz(y, x=z) / span_length

    return np.mean(y)

# ==========================================
if __name__ == "__main__":
    main()
