# Nonlinear pH Model of Sulferic Acid and NaOH Mixing

This project was completed as part of the coursework for the 
Chemical Engineering undergradudate course, Process Controls Lab (CME-453L)

The experimental setup consisted of diluted Sulferic Acid (H2SO4) and
Sodium Hydroxide (NaOH) pumped from container vats by diaphragm pumps
into an agitated round-bottom flask with an overflow outlet, and a pH meter.
See "visuals" for experimental setup visual. 

## Usage 

1. Install requirements:
    pip install -r requirements.txt

2. Run `main.py`

3. You will be prompted to select a CSV file. Sample data is provided in sample_data.

4. The initial guess for the initial concentrations of SO4^2- (A2-) and Na+ (B+) in 
    the vial is built into main.py. This may need to be edited based on the csv file, 
    as the minimize function can be finicky.  

5. Graphical output of predicted vs actual pH, [H+], [OH-], etc. Example output in "visuals"

## Requirements

- matplotlib
- tkinter
- numpy
- csv
- scipy
- scikit-learn

## File Structure

- `main.py` : Main runner
- `import_csv_data.py` : Imports acid/base flow rate and pH data 
- `calculations.py ` : Nonlinear model calculations
- `gui.py` : Matplotlib graph display

## Author

Caroline Howley – c.howley3@yahoo.com
