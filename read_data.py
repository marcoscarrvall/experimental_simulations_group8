import pandas as pd
"""
DON'T CHANGE THIS FILE UNLESS YOU KNOW WHAT YOU ARE DOING
ALSO DON'T CHANGE FULL_DATA.TXT UNLESS YOU KNOW WHAT YOU ARE DOING
Function to read the data from the text file and convert it into a pandas DataFrame.
Columns: 
- Run_nr: the number of the run
- Alpha: the angle of attack
- CL: the lift coefficient
- CD: the drag coefficient
- Cyaw: the yawing moment coefficient
- Cm_p_qc: the pitching moment coefficient
- Ct: the thrust coefficient
- Cn: the rolling moment coefficient
- Cside: the side force coefficient
- Cm_roll: the rolling moment coefficient
- Cm_pitch: the pitching moment coefficient
- Cm_yaw: the yawing moment coefficient
- V: the velocity
- Re: the Reynolds number
- M: the Mach number
- delta_e: the elevator deflection
- J: the advance ratio

To input function add this at the top of your code:
from read_data import process_data
To find for example CL, we run process_data()
and then we can access the CL column by data['CL']
"""

def process_data(text_file = 'Data/full_data.txt'):
    with open(text_file, 'r') as f:
        text_file = f.read()

    columns = {}
    
    lines = text_file.strip().split('\n')   # Read line by line
    headers = lines[0].split()  # Extract headers from the first line

    for h in headers:
        columns[h] = []     # Set header names as keys and initialize empty lists for values
        
    for line in lines[2:-5]:
        parts = line.split()
        if not parts:
            continue
            
        
        for i in range(0, len(headers)):
            val_str = parts[i]
            try:
                val = float(val_str)  # Try to convert to float
            except ValueError:
                val = val_str  # Otherwise keep as string
            columns[headers[i]].append(val)
    data = pd.DataFrame(columns)
    return data

