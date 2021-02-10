import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk')

def load_data(file_dir, set_name):
    dict = {}
    key = ""
    for i in range(1, 101):
        file_path = ""
        if i < 10:
            key = f"{set_name}00{i}"
            file_path = f"{file_dir}/{key}.txt"
        elif i == 100:
            key = f"{set_name}100"
            file_path = f"{file_dir}/{key}.txt"
        else:
            key = f"{set_name}0{i}"
            file_path = f"{file_dir}/{key}.txt"

        with open(file_path) as eeg_file:
            tmp = eeg_file.readlines()
            data = [float(k) for k in tmp]
            dict[f'{set_name}{i}'] = data

    return dict

def plot_segment(data_set, set_name, start_seg=1, num_seg=5, sample_no=4097):
    fig, axes = plt.subplots(num_seg, figsize=(20, 20))
    axes_index = 0
    for i in range(start_seg, start_seg + num_seg):
        series = pd.Series(data_set[f'{set_name}{i}'])
        axes[axes_index].plot(series)
        axes_index += 1

