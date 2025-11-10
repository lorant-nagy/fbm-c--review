"""
FBM Analysis Script
Analyzes fractional Brownian motion data generated from C code and compares with Python-generated control data.
All plots are saved to files with a run identifier.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from statsmodels.tsa.stattools import adfuller
import seaborn as sns
from fbm import FBM

# =============================================================================
# USER CONFIGURATION
# =============================================================================
RUN_IDENTIFIER = "H0.8_T128_n128_run2"
BASE_DIR = "/Users/lorantnagy/Documents/CODE/MLDTRADE/fbm-c--review/data/fbm_H0.8_T128.0_n128"
OUTPUT_DIR = "./python/analysis_results_CImp_patched"
# =============================================================================

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)


def read_data(BASE_DIR=None):
    """Read FBM data from CSV files in the specified directory."""
    file_count = 0
    for file in os.listdir(BASE_DIR):
        if file.endswith('.csv'):
            file_count += 1

    # Read all files and store DataFrames in a list
    dataframes = []
    time_column = None
    
    for i in range(1, file_count + 1):
        df = pd.read_csv(f'{BASE_DIR}/fbm_{i}.csv', sep=',')
        if time_column is None:
            time_column = df['t'].copy()  # Store time column from first file
        # Create a DataFrame with just the path data and proper column name
        path_df = pd.DataFrame({f'path_{i}': df['BH']})
        dataframes.append(path_df)
    
    # Concatenate all path DataFrames at once
    if dataframes:
        data = pd.concat(dataframes, axis=1)
        # Add time column at the beginning
        data.insert(0, 'time', time_column)
    else:
        data = pd.DataFrame()

    # extract H T n from the base dir
    parts = BASE_DIR.split('_')
    H = float(parts[1][1:])
    T = float(parts[2][1:])
    n = int(parts[3][1:])

    return data, (H, T, n)


def _calc_diff(data):
    """Calculate increments (differences) of paths."""
    path_columns = [col for col in data.columns if col.startswith('path_')]
    data_increments = data[path_columns].diff()
    data_increments = data_increments.dropna()
    return data_increments.values


def save_figure(filename):
    """Save figure with run identifier in filename."""
    full_path = os.path.join(OUTPUT_DIR, f"{RUN_IDENTIFIER}_{filename}")
    plt.savefig(full_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {full_path}")
    plt.close()


def main():
    print("="*80)
    print("FBM ANALYSIS")
    print("="*80)
    print(f"Run Identifier: {RUN_IDENTIFIER}")
    print(f"Data Directory: {BASE_DIR}")
    print(f"Output Directory: {OUTPUT_DIR}")
    print("="*80)
    
    # Load data
    print("\n[1/9] Loading data...")
    data, (H, T, n) = read_data(BASE_DIR=BASE_DIR)
    number_of_paths = data.shape[1] - 1
    print(f"  Parameters: H={H}, T={T}, n={n}")
    print(f"  Number of paths: {number_of_paths}")
    
    # Generate control data using Python FBM
    print("\n[2/9] Generating control data with Python FBM...")
    f = FBM(n=n, hurst=H, length=T, method='daviesharte')
    t = np.linspace(0, T, n+1)
    fbm_paths = np.zeros((n+1, number_of_paths))
    for i in range(number_of_paths):
        fbm_paths[:, i] = f.fbm()

    path_dataframes = []
    for i in range(number_of_paths):
        path_df = pd.DataFrame({f'path_{i+1}': fbm_paths[:, i]})
        path_dataframes.append(path_df)

    data_control = pd.concat(path_dataframes, axis=1)
    data_control.insert(0, 'time', data['time'])
    
    # Calculate envelope
    print("\n[3/9] Calculating theoretical envelope...")
    time_stamps = data.time.values
    n_fill = 4
    time_stamps_wth0_plus_e = time_stamps[n_fill:].copy()
    envelop = time_stamps_wth0_plus_e**H * np.sqrt(np.log(np.log(time_stamps_wth0_plus_e + 1e-10)) + 1e-10)
    envelop = np.concatenate((np.zeros(n_fill), envelop))
    
    # Plot trajectories
    print("\n[4/9] Plotting trajectories...")
    N_sample = 15
    
    plt.figure(figsize=(4, 3))
    for i in range(1, N_sample):
        plt.plot(data['time'], data[f'path_{i}'], linewidth=0.5, alpha=0.8)
    plt.plot(time_stamps, envelop, 'k--', label='envelope +')
    plt.plot(time_stamps, -envelop, 'k--', label='envelope -')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.title('Data generated with C')
    plt.legend()
    save_figure("trajectories_C.png")

    plt.figure(figsize=(4, 3))
    for i in range(1, N_sample):
        plt.plot(data_control['time'], data_control[f'path_{i}'], linewidth=0.5, alpha=0.8)
    plt.plot(time_stamps, envelop, 'k--', label='envelope +')
    plt.plot(time_stamps, -envelop, 'k--', label='envelope -')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.title('Data generated with Python')
    plt.legend()
    save_figure("trajectories_Python.png")
    
    # Plot standard deviation
    print("\n[5/9] Plotting standard deviation comparison...")
    std = data.iloc[:, 1:].std(axis=1)
    t_ad_H = time_stamps**(H)
    std_control = data_control.iloc[:, 1:].std(axis=1)

    plt.figure()
    plt.plot(time_stamps, t_ad_H, label="t^H")
    plt.plot(time_stamps, std, label="std (C)")
    plt.plot(time_stamps, std_control, label="std (Python)")
    plt.xlabel('Time')
    plt.ylabel('Standard Deviation')
    plt.title('Standard Deviation Comparison')
    plt.legend()
    plt.grid(True)
    save_figure("std_comparison.png")
    
    # Calculate increments
    print("\n[6/9] Calculating increments and correlations...")
    data_diff = _calc_diff(data)
    control_diff = _calc_diff(data_control)
    
    # Stationarity tests
    print("\n[7/9] Running stationarity tests...")
    path_1_diff_C = data_diff[:, 0]
    path_1_diff_control = control_diff[:, 0]

    print("\n  Stationarity test for C generated data:")
    adf_result = adfuller(path_1_diff_C)
    print(f"    ADF Statistic: {adf_result[0]:.6f}")
    print(f"    p-value: {adf_result[1]:.6f}")
    print(f"    Critical Values: {adf_result[4]}")
    if adf_result[1] < 0.05:
        print("    Result: Increments are stationary (reject null hypothesis)")
    else:
        print("    Result: Increments may be non-stationary (fail to reject null hypothesis)")

    print("\n  Stationarity test for Python generated data:")
    adf_result = adfuller(path_1_diff_control)
    print(f"    ADF Statistic: {adf_result[0]:.6f}")
    print(f"    p-value: {adf_result[1]:.6f}")
    print(f"    Critical Values: {adf_result[4]}")
    if adf_result[1] < 0.05:
        print("    Result: Increments are stationary (reject null hypothesis)")
    else:
        print("    Result: Increments may be non-stationary (fail to reject null hypothesis)")
    
    # Autocorrelation analysis
    print("\n[8/9] Computing autocorrelation functions...")
    
    # C generated data
    corrcoeffs = {}
    left = []
    right = []
    for l in range(0, 10):
        for i in range(data_diff.shape[1]):
            path = data_diff[:, i]
            for j in range(len(path) - l):
                left.append(path[j])
                right.append(path[j+l])
        corrcoeffs[l] = np.corrcoef(left, right)[0, 1]
        left = []
        right = []

    plt.figure()
    plt.plot(list(corrcoeffs.keys()), list(corrcoeffs.values()), marker='o', label='empirical decay')
    plt.plot(list(corrcoeffs.keys()), 
             np.insert(H*(2*H-1)*(np.array(list(corrcoeffs.keys()))[1:])**(2*H-2), 0, 1.0), 
             marker='x', label='theoretical decay (∼ l^{2H-2})')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation Coefficient')
    plt.title('Autocorrelation of Increments (C generated data)')
    plt.grid()
    plt.legend()
    save_figure("autocorrelation_C.png")

    # Python generated data
    corrcoeffs = {}
    left = []
    right = []
    for l in range(0, 10):
        for i in range(control_diff.shape[1]):
            path = control_diff[:, i]
            for j in range(len(path) - l):
                left.append(path[j])
                right.append(path[j+l])
        corrcoeffs[l] = np.corrcoef(left, right)[0, 1]
        left = []
        right = []

    plt.figure()
    plt.plot(list(corrcoeffs.keys()), list(corrcoeffs.values()), marker='o', label='empirical decay')
    plt.plot(list(corrcoeffs.keys()), 
             np.insert(H*(2*H-1)*(np.array(list(corrcoeffs.keys()))[1:])**(2*H-2), 0, 1.0), 
             marker='x', label='theoretical decay (∼ l^{2H-2})')
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation Coefficient')
    plt.title('Autocorrelation of Increments (Python generated data)')
    plt.grid()
    plt.legend()
    save_figure("autocorrelation_Python.png")
    
    # Histograms
    print("\n[9/9] Creating histograms...")
    
    # C generated only
    plt.figure(figsize=(6, 4))
    plt.hist(data_diff.flatten(), bins=50, alpha=0.7, label='C generated data')
    plt.xlabel('Increment Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Increments (C generated data)')
    plt.legend()
    save_figure("histogram_C.png")

    # Python generated only
    plt.figure(figsize=(6, 4))
    plt.hist(control_diff.flatten(), bins=50, alpha=0.7, label='Python generated data')
    plt.xlabel('Increment Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Increments (Python generated data)')
    plt.legend()
    save_figure("histogram_Python.png")

    # Comparison histogram
    plt.figure(figsize=(6, 4))
    plt.hist(data_diff.flatten(), bins=50, alpha=0.5, label='C generated data')
    plt.hist(control_diff.flatten(), bins=50, alpha=0.5, label='Python generated data')
    plt.xlabel('Increment Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Increments (Comparison)')
    plt.legend()
    save_figure("histogram_comparison.png")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print(f"All plots saved to: {OUTPUT_DIR}")
    print("="*80)


if __name__ == "__main__":
    main()