import numpy as np
import math
import itertools
from decimal import Decimal
import pandas as pd
import tkinter as tk
from tkinter import filedialog

import calculations
import import_csv
import gui

def main():
    # Initialize Tkinter but hide the main window
    root = tk.Tk()
    root.withdraw()
    # Prompt user to select a file
    file_path = filedialog.askopenfilename(title="Select CSV File", filetypes=[("CSV Files", "*.csv")])
    root.destroy()
    #file_name = "2025.02.24-acid-only.csv"
    
    # Prompt for Initial Guess for [A2-] and [B+] in vial
    # "Note: Converging on a solution may require iterating initial guesses"
    #acid only flow: 
    #
    # Guess for C_a_vial and C_b_vial
    guess = [0.00200439, 0.004008]
    # Concentration of A2- and B+ in the acid and base feeds, respectively
    C_a_feed = 0.00275
    C_b_feed = 0.02

    data_dictionary = import_csv.import_data_from_csv(file_path)
    model = calculations.solve_command(guess, C_a_feed, C_b_feed, data_dictionary)
    print(f"Model for {file_path}\n")  
    gui.display_data(data_dictionary, model)
  
if __name__ == "__main__":
    main()

