import csv
import pandas as pd
### --- Access Experimental Data from CSV File ---
def import_data_from_csv(file_name):
    data_dictionary = {
        'time': [],
        'acid_flow': [],
        'base_flow': [],
        'pH': [], 
        '[H+]': [], 
        '[OH-]': []
    }
        
    # Open the CSV file and read the data
    with open(file_name, mode='r') as file:
        reader = csv.reader(file) 
        # Skip the header 
        next(reader, None)
        first_row = next(reader)
        time_initial = float(first_row[1]) if first_row[1] else 0.0

        for row in reader:
            if not row or all(cell == '' for cell in row):  
                break  # Exit the loop if an empty or blank line is encountered

            time = float(row[1]) if row[1] else float(0)
            # Default flow rates and pH to zero if empty
            acid_flow = float(row[4]) if row[4] else float(0)  
            base_flow = float(row[5]) if row[5] else float(0)  
            pH = float(row[6]) if row[6] else float(0)  
    
            # Append the data as a dictionary
            data_dictionary['time'].append(time - time_initial)
            data_dictionary['acid_flow'].append(acid_flow / 1000)
            data_dictionary['base_flow'].append(base_flow / 1000)
            data_dictionary['pH'].append(pH)
            data_dictionary['[H+]'].append(10**(-pH))
            data_dictionary['[OH-]'].append(10**(-(14 - pH)))
    return data_dictionary
