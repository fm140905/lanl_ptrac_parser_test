# Description: This script reads the output file and plots the histogram of deposited energies
import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
    data = np.loadtxt("photon_energy_deposition_events.txt", skiprows=1)
    deposited_energies = data[:, 6]
    weights = data[:, 8]
    # draw histogram of deposited energies
    counts, bins = np.histogram(deposited_energies, bins=100, weights=weights)
    # step-style histogram
    fig, ax = plt.subplots()
    ax.step(bins[:-1], counts, where='post')
    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel("Counts")
    plt.show()
