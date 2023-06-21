import numpy as np
from matplotlib import pyplot as plt
from ssparser import parse_dataset_freq

from constants import *

"""
target_residues = []
for i, val in enumerate(beta[:-2]):
    if val / total_chain_count > threshold:
        target_residues.append(i+1)
"""

def plot(plot_dict: dict):

    x = np.arange(ABETA_LONG_RESIDUE_COUNT)  
    width = 0.25  # the width of the bars
    multiplier = 3

    fig, ax = plt.subplots(layout='constrained')
    plt.xticks(np.arange(1,43, dtype=int), ABETA_LONG_SEQUENCE, rotation=-60)

    # add line at 50%
    plt.axhline(y=0.5, color='r', linestyle='-')

    for index, (attribute, measurement) in enumerate(plot_dict.items()):
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        multiplier += 1

    ax.set_ylabel("Frequency of Secondary Structure")
    ax.set_xlabel("Residue ID")
    ax.set_title("Frequency of Secondary Structure per Residue of ABeta-Monomer")
    ax.legend(loc='upper right')

    plt.show()

plot(parse_dataset_freq(dataset_dir=DATASET_DIR))