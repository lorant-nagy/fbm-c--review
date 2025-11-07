import pandas as pd
import os

BASE_DIR = "/Users/lorantnagy/Documents/CODE/MLDTRADE/fbm-c--review/data"

def read_data(BASE_DIR = BASE_DIR):
    file_count = 0
    for file in os.listdir(BASE_DIR):
        if file.endswith('.csv'):
            file_count += 1

    data = pd.DataFrame()
    for i in range(1, file_count + 1):
        df = pd.read_csv(f'../../data/fbm_{i}.csv', sep=',')
        if 'time' not in data.columns:
            data['time'] = df['t']
        data[f'path_{i}'] = df['BH']

    return data
