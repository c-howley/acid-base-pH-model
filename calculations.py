from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from sklearn.metrics import r2_score
import numpy as np

# requirements:
#scipy
#scikit-learn

### --- Constants ---
V = 1.7 #Liters
AcidpH = 1.8
BasepH = 11.5
Kw_eq = (10**-14)
data_dictionary = {}
model_dict = {
    "t": [], 
    "pH": [],
    "[H+]": [],
    "[OH-]": [],
    "[H+]i": [],
    "[OH-]i": [],
    "[A2-]": [],
    "[HA-]": [],
    "[B+]": [],
    "pOH": [], 
    "Exp. pH (Interpd to sol.t)": [], 
    "Error": []
}

### --- Nonlinear pH Model ---
def model(t, z, C_a_feed, C_b_feed):
    C_a, C_b = z
    #Differential Equations
    #Format:
        #Line 1: Inputs
        #Line 2: Formation/Reaction
        #Line 3: Outputs
    dC_adt = ((wA(t)* C_a_feed)  
              + 0
              - (wE(t) * C_a))/V
    
    dC_bdt = ((wB(t) * C_b_feed)  
              + 0
             - (wE(t) * C_b))/V
    return dC_adt, dC_bdt

### --- Find data at a particular time (with interpolation) ---
def get_data_by_time(time_value):
    times = data_dictionary['time']
    acid_flows = data_dictionary['acid_flow']
    base_flows = data_dictionary['base_flow']
    pHs = data_dictionary['pH']

    # Check if time_value exactly exists
    if time_value in times:
        idx = times.index(time_value)
        return acid_flows[idx], base_flows[idx], pHs[idx]

    # If not, find closest indices around time_value
    times_arr = np.array(times)
    if time_value < times_arr[0] or time_value > times_arr[-1]:
        raise ValueError(f"Time {time_value} out of data range.")

    # Find indices bounding time_value
    idx_upper = np.searchsorted(times_arr, time_value, side='right')
    idx_lower = idx_upper - 1

    # Linear interpolate values
    t0, t1 = times_arr[idx_lower], times_arr[idx_upper]
    alpha = (time_value - t0) / (t1 - t0)
    acid_flow_interp = acid_flows[idx_lower] + alpha * (acid_flows[idx_upper] - acid_flows[idx_lower])
    base_flow_interp = base_flows[idx_lower] + alpha * (base_flows[idx_upper] - base_flows[idx_lower])
    pH_interp = pHs[idx_lower] + alpha * (pHs[idx_upper] - pHs[idx_lower])

    return acid_flow_interp, base_flow_interp, pH_interp

### --- Forcing Functions ---
def wA(t):
    return float(get_data_by_time(t)[0])

def wB(t):
    return float(get_data_by_time(t)[1])

def wE(t):
    return wB(t) + wA(t)

# --- Calculate equilibrium concentrations of H+, O-, and H20 given ions ---
def water_formation_equilibrium(h, oh, Kw_eq):
        # Define quadratic eqn constants
        # based on Kw_eq = ([H+]-x)([OH-]-x)
        a = 1
        b = -(h + oh)
        c = h * oh - Kw_eq
    
        # Solve for x using discriminant
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            raise ValueError("No real roots exist.")

        root1 = (-b + np.sqrt(discriminant)) / (2 * a)
        root2 = (-b - np.sqrt(discriminant)) / (2 * a)
    
        # Condition: root must be positive and less than min(h, oh)
            # i.e: reaction always tends towards water formation &
            # the amt of water formed cannot be greater than the H or OH present
        threshold = min(h, oh)
        valid_roots = [r for r in (root1, root2) if r >= 0 and r < threshold]
    
        if not valid_roots:
            return np.nan
    
        # Return the smallest valid root (closest to zero)
        return max(valid_roots)

### --- Function to Calculate pH from A2- and B+ ---
def pH_from_A_B(a, b, source):
    pH_model = []

    #RETURN - add this to the model dictionary instead, and pass that
    deprotonation_factor = 2
    protonation_factor = 1

    #Create list of initial Ch and Coh concentrations
    for i in range(len(a)):
        h = a[i] * deprotonation_factor
        oh = b[i] * protonation_factor
        ha = a[i]* (2 - deprotonation_factor)

        # Only append to dictionary if called by main solve function
        if source == "main":
            model_dict["[HA-]"].append(ha)
            model_dict["[OH-]i"].append(oh)
            model_dict["[H+]i"].append(h)
            model_dict["[A2-]"].append(a[i])
            model_dict["[B+]"].append(b[i])

        # Solve quadratic eqn for change in [H+] and [OH-]
        x = water_formation_equilibrium(h, oh, Kw_eq)

        # Append NaN to model if quadratic eqn cannot be solved
        if np.isnan(x):
            if source == "main":
                model_dict["[H+]"].append(np.Nan)
                model_dict["[OH-]"].append(np.Nan)
                model_dict["pH"].append(np.Nan)
                model_dict["pOH"].append(np.Nan)
        # Else, append solutions
        else:
            # If x is larger than [H+]i, throw error
            assert (h-x) > 0, f"Error at i = {i}: x = {x} must be smaller \
            than [H]i = {h}. (h-x) = {h-x} for a2- = {a[i]} and b+ = {b[i]}."
            # If x is larger than [OH-]i, throw error
            assert (oh-x) > 0, f"Error at i = {i}: x = {x} must be smaller \
            than [OH]i = {oh}. (oh-x) = {oh-x} for a2- = {a[i]} and b+ = {b[i]}."
                
            positive_charges = h + (b[i])
            negative_charges = oh + (2 * a[i]) + ha
            vial_charge_concentration = positive_charges - negative_charges
            assert abs(vial_charge_concentration) < 10**-7, \
                f"Error: Solver calculated a nonzero vial charge, \
                    {vial_charge_concentration:.5e}, at i = {i}"
            
            pH = -np.log10(h - x)
            pOH = -np.log10(oh - x)
            pH_model.append(pH)

            if source == "main": 
                model_dict["[H+]"].append(h - x)
                model_dict["[OH-]"].append(oh - x)
                model_dict["pH"].append(pH)
                model_dict["pOH"].append(pOH)
    return pH_model

### --- Minimize error function ---
def minimize_error(guess):
    intervals = 200
    t = np.linspace(data_dictionary['time'][0], data_dictionary['time'][-2], intervals)

    def rmse_calculation(IC, t, pH_data):
        sol = solve_ivp(
            model,
            t_span=[t[0], t[3]],
            y0=IC,
            max_step=10,
            args=(C_a_feed, C_b_feed),
            t_eval=t[0:3]
        )
        # how to add to model within function but not include minimize stuff:
        # pass a dictionary name to add to
        # only need pH_model
        pH_model = pH_from_A_B(sol.y[0], sol.y[1], "minimize")
        if np.isnan(pH_model).any():
            return 1e2
        
        pH_data_norm = np.array([get_data_by_time(time)[2] for time in sol.t])
        sum_error_sqd = np.sum((np.array(pH_model) - pH_data_norm) ** 2)
        return sum_error_sqd

    test_error = rmse_calculation(guess, t, data_dictionary["pH"])

    if test_error > 10:
        result = minimize(
            rmse_calculation,
            x0=guess,
            args=(t, data_dictionary["pH"]),
            bounds=[(1e-14, 0.1), (1e-14, 0.1)],
            options={'maxiter': 100, 'disp': True, 'gtol': 1e-3, 'eps': 1e-8}
        )
        IC = result.x
        assert result.success == True, f"Error: Minimize unsuccessful. {result}"
    else:
        IC = guess
    return IC

### --- Error calculation ---
def error(pH_data, pH_model):
    error = list(np.array(pH_data) - np.array(pH_model))
    rmse = np.sqrt(np.mean(np.square(np.array(pH_data) - np.array(pH_model))))
    r2 = r2_score(pH_data, pH_model)
    return error, r2, rmse

### --- Main solve command ---
def solve_command(guess, Caf, Cbf, data_dict):
    # Declare data and initial guesses as globals
    global data_dictionary
    global C_a_feed
    global C_b_feed
    data_dictionary, C_a_feed, C_b_feed = data_dict, Caf, Cbf

    # Solve for initial [A2-] & [B+] in vial using guess
    vial_IC = minimize_error(guess)

    # Solve linear model
    intervals = 200
    t = np.linspace(data_dictionary['time'][0], data_dictionary['time'][-2], intervals)
    sol = solve_ivp(
        model,
        t_span=[t[0], t[-2]],
        y0=vial_IC,
        max_step=5,
        args=(C_a_feed, C_b_feed),
    )
    model_dict["t"] = sol.t

    # Find pH from A2- and B+ (values appended to model_dict in function)
    pH_from_A_B(sol.y[0], sol.y[1], "main")
    
    # Interpolate pH data to find values at solveivp timesteps
    pH_data_solveivp_timestep = np.array([
        get_data_by_time(time)[2]
        for time in sol.t
    ])
    model_dict["Exp. pH (Interpd to sol.t)"] = pH_data_solveivp_timestep 

    # Solve for error stats
    model_error, r2, rmse = error(pH_data_solveivp_timestep, model_dict["pH"])
    model_dict["Error"] = model_error

    return model_dict
