import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk

def display_data(data, model):
    #Generate Plots
    # Plot collected data
    try:
        data_to_plot = [
            [0, 0, data["time"], data["acid_flow"], "r", "Acid Feed", "Flow Rate (mL/s)"], 
            [0, 0, data["time"], data["base_flow"], "b", "Base Feed", None],

            [0, 1, data["time"], data["pH"], "k", "Experimental pH", None], 
            [0, 1, model["t"], model["pH"], "m.", "Predicted pH", None],

            [0, 2, model["t"], model["Error"], "o", "Model Error", "Experimental pH - Model pH"],

            [1, 0, model["t"], model["[H+]"], "r.", "Predicted [H+]", "Concentration (mol/L)"],
            [1, 0, data["time"], data["[H+]"], "k", "Experimental [H+]", None],

            [1, 1, model["t"], model["[OH-]"], "b.", "Predicted [OH-]", None], 
            [1, 1, data["time"], data["[OH-]"], "k", "Experimental [OH-]", None],

            [1, 2, model["t"], model["[A2-]"], "m.", "Predicted [A2-]", None], 
            [1, 2, model["t"], model["[B+]"], "g.", "Predicted [B+]", None]
        ]
    except Exception as e:
        print(f"Error parsing data to plot: {e}")
        return

    fig, ax = plt.subplots(2, 3, figsize = (12, 6))
    for i, j, x, y, marker, label, ylabel in data_to_plot:
        ax[i, j].plot(x, y, marker, markersize = 3, label = label)

        ax[i, j].set_xlabel("Time (s)")
        ax[i, j].legend()
        if ylabel is not None:
            ax[i, j].set_ylabel(ylabel, labelpad = 10)

    plt.tight_layout()
    plt.show()
    return
