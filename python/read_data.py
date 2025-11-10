import pandas as pd
import os

# this is how we create the dir sprintf(dir, "data/fbm_H%0.1f_T%0.1f_n%d", H, T, n);
# sample /Users/lorantnagy/Documents/CODE/MLDTRADE/fbm-c--review/data/fbm_H0.7_T128.0_n128

def read_data(BASE_DIR = None):
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
