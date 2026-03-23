import pandas as pd

"""
DON'T CHANGE THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING
ALSO DON'T CHANGE FULL_DATA.TXT UNLESS YOU KNOW WHAT YOU ARE DOING
Function to read the data from the CSV file and convert it into a pandas DataFrame.
Columns:
- config: configuration identifier
- de: elevator deflection
- dr: rudder deflection
- run: run number
- AoA: angle of attack
- AoS: angle of sideslip
- V: velocity
- Re: Reynolds number
- rpsM2: RPS * M^2
- J: advance ratio
- CL: lift coefficient
- CD: drag coefficient
- CM25c: pitching moment coefficient at 25% chord
- CYaw: yawing moment coefficient
- eps_total: total downwash angle
- eps_sb: sidebody downwash angle
- eps_wb0: wing-body downwash angle
- CL_w: wing lift coefficient
- CL_alpha_mapped: mapped CL alpha
- Delta_Alpha: alpha correction
- Delta_CD: drag correction
- Delta_Cm: moment correction
- AoA_corrected: corrected angle of attack
- CD_corrected: corrected drag coefficient
- CM25c_corrected: corrected pitching moment coefficient

To use this function, add this at the top of your code:
    from read_data import process_data
Then call:
    data = process_data()
And access columns like:
    data['CL']
"""

def process_data(csv_file='fully_corrected_data.csv'):
    data = pd.read_csv(csv_file)

    # Convert all columns to numeric where possible, keeping strings otherwise
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='ignore')

    return data