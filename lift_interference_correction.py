import pandas as pd
import numpy as np

# 1. Define filenames
main_file = 'blockage_corrected_data.csv'   # Your main data file
# Updated path based on your folder structure!
tail_off_file = r'TAILOFF\TAILOFF/unc_tailOff_beta0_balance.txt' 

# 2. Load datasets
df_main = pd.read_csv(main_file)
# We use skiprows=[1] to skip the second row that contains the units!
df_tail_off = pd.read_csv(tail_off_file, sep='\s+', skiprows=[1])


# --- THE MATCHING PROCESS ---

# 3. Create rounded matching columns in BOTH datasets
df_main['AoA_match'] = df_main['AoA'].round().astype(int)
df_tail_off['AoA_match'] = df_tail_off['Alpha'].round().astype(int)

# Velocity rounding works perfectly (e.g., 19.38 -> 20, 38.76 -> 40)
df_main['V_match'] = (df_main['V'] / 10).round().astype(int) * 10
df_tail_off['V_match'] = (df_tail_off['V'] / 10).round().astype(int) * 10

# 4. Isolate the CL column from the tail-off data
df_tail_off_subset = df_tail_off[['AoA_match', 'V_match', 'CL']].copy()
df_tail_off_subset.rename(columns={'CL': 'CL_w'}, inplace=True)

# 5. Merge the datasets (Tail-off CL)
df = pd.merge(df_main, df_tail_off_subset, on=['AoA_match', 'V_match'], how='left')


# --- NEW: APPLYING THE C_L_alpha LOOKUP TABLE (V & J only, de=0) ---

# Create a dataframe using your 6 exact V and J combinations for de=0
cla_data = {
    'V_match': [40, 40, 40, 20, 20, 20],
    'J_match': [1.6, 1.9, 2.2, 1.6, 1.9, 2.2],
    'CL_alpha_mapped': [
        0.107854,  # <-- exact value for V=40, J=1.6
        0.106074,  # <-- exact value for V=40, J=1.9
        0.105445,  # <-- exact value for V=40, J=2.2
        0.106699,  # <-- exact value for V=20, J=1.6
        0.104644,  # <-- exact value for V=20, J=1.9
        0.103870  # <-- exact value for V=20, J=2.2
    ] 
}
df_cla = pd.DataFrame(cla_data)

# --- THE LENIENT ROUNDING FIX FOR J ---
# Find the absolute closest nominal J value (1.6, 1.9, or 2.2) for each row
nominal_Js = [1.6, 1.9, 2.2]
df['J_match'] = df['J'].apply(lambda x: min(nominal_Js, key=lambda target: abs(target - x)) if pd.notna(x) else np.nan)

# Merge the correct CLa values into the main dataframe
df = pd.merge(df, df_cla, on=['V_match', 'J_match'], how='left')

# NOW you can safely clean up all the matching columns at once!
df.drop(columns=['AoA_match', 'V_match', 'J_match'], inplace=True)


# --- THE CORRECTION MATH ---

# 6. Define your constants
S = 0.2172   
St = 0.0858
ARt = 3.87             
C = 2.07                # <-- UPDATE THIS to your tunnel cross-sectional area
delta = 0.1033            
tau_2 = 0.120      
tau_2t = 0.75       
Cmac =  0.165 
lt = 0.535

dCm_dalpha_t = - ((St*lt)/(S*Cmac))*((0.1*ARt*0.8)/(ARt+2))  # Using CL_alpha = 0.1*AR/(AR+2) for the tail

area_ratio = S / C

# Assign the dynamically mapped column to your Cla variable
Cla = df['CL_alpha_mapped']

CMuw = (1/8)*tau_2*(0.5*Cmac)*delta*area_ratio*df['CL_w']*Cla 
CMt = dCm_dalpha_t*delta*area_ratio*df['CL_w']*(1+tau_2t*lt)

# 7. Calculate the correction terms using the newly matched 'CL_w'
df['Delta_Alpha'] = delta * area_ratio * df['CL_w'] * (1+tau_2*(0.5*Cmac))*57.3 
df['Delta_CD'] = delta * area_ratio * (df['CL_w'] ** 2)
df['Delta_Cm'] = (CMuw+CMt) 

# 8. Apply the corrections by OVERWRITING the original columns
df['AoA'] = df['AoA'] + df['Delta_Alpha']
df['CD'] = df['CD'] + df['Delta_CD']
df['CM25c'] = df['CM25c'] + df['Delta_Cm']  

# Check for missing matches (Tail-off CL_w)
missing_matches = df['CL_w'].isna().sum()
if missing_matches > 0:
    print(f"WARNING: {missing_matches} rows could not find a matching CL_w in the tail-off data.")

# Check for missing slopes (C_L_alpha)
missing_slopes = df['CL_alpha_mapped'].isna().sum()
if missing_slopes > 0:
    print(f"WARNING: {missing_slopes} rows could not find a matching CL_alpha in the lookup table. Check your J=0 or wind-off runs.")

# --- FINAL FORMATTING ---

# 9. Filter to keep exactly the requested columns
final_columns = [
    'config', 'de', 'dr', 'run', 'AoA', 'AoS', 'V', 'Re', 
    'rpsM2', 'J', 'CL', 'CD', 'CM25c', 'CYaw', 
    'eps_total', 'eps_sb', 'eps_wb0'
]

# Create the final dataframe with only these columns
df_final = df[final_columns]

# 10. Save back to the main file (or a new file)
output_filename = 'fully_corrected_data.csv'
df_final.to_csv(output_filename, index=False)

print(f"Data successfully matched and corrections applied! Saved to {output_filename}")