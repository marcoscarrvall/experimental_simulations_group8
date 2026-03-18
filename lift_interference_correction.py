import pandas as pd

# 1. Define filenames
main_file = 'blockage_corrected_data.csv'   # Your main data file

# Updated path based on your folder structure!
tail_off_file = 'TAILOFF\TAILOFF/unc_tailOff_beta0_balance.txt' 

# 2. Load datasets
df_main = pd.read_csv(main_file)

# We use skiprows=[1] to skip the second row that contains the units!
df_tail_off = pd.read_csv(tail_off_file, sep='\s+', skiprows=[1])


# --- THE MATCHING PROCESS ---

# 3. Create rounded matching columns in BOTH datasets
df_main['AoA_match'] = df_main['AoA'].round().astype(int)
df_tail_off['AoA_match'] = df_tail_off['Alpha'].round().astype(int)

df_main['V_match'] = (df_main['V'] / 10).round().astype(int) * 10
df_tail_off['V_match'] = (df_tail_off['V'] / 10).round().astype(int) * 10

# 4. Isolate the CL column from the tail-off data
df_tail_off_subset = df_tail_off[['AoA_match', 'V_match', 'CL']].copy()
df_tail_off_subset.rename(columns={'CL': 'CL_w'}, inplace=True)

# 5. Merge the datasets
df = pd.merge(df_main, df_tail_off_subset, on=['AoA_match', 'V_match'], how='left')
df.drop(columns=['AoA_match', 'V_match'], inplace=True)


# --- NEW: APPLYING THE C_L_alpha LOOKUP TABLE ---
# Create a dataframe using the exact values from your Excel image
# --- NEW: APPLYING THE C_L_alpha LOOKUP TABLE ---
# Added the 'de0dr0_extra' configuration to catch those missing 10 rows!
cla_data = {
    'config': ['de0dr0', 'de0dr0', 'de0dr0', 'de0dr0', 'de0dr0', 
               'de5dr0', 'de5dr0', 'de5dr0', 
               'dem5dr0', 'dem5dr0', 'dem5dr0',
               'de0dr0_extra', 'de0dr0_extra'], # <-- Added missing config
    'J_match': [1.5, 1.6, 1.8, 1.9, 2.2, 
                1.6, 1.9, 2.2, 
                1.6, 1.9, 2.2,
                1.6, 2.2],                      # <-- J values for the extra runs
    'CL_alpha_mapped': [0.11714, 0.10850, 0.10664, 0.10427, 0.10432, 
                        0.11304, 0.10542, 0.10260, 
                        0.10580, 0.10505, 0.10443,
                        0.10850, 0.10432]       # <-- Copied slopes from de0dr0
}
df_cla = pd.DataFrame(cla_data)

# Round J to 1 decimal place in the main data so it matches the table
df['J_match'] = df['J'].round(1)

# Merge the correct CLa values into the main dataframe based on config and J
df = pd.merge(df, df_cla, on=['config', 'J_match'], how='left')
df.drop(columns=['J_match'], inplace=True) # clean up


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
df['Delta_Cm'] = (CMuw+CMt) * 57.3

# 8. Apply the corrections to create new columns
df['AoA_corrected'] = df['AoA'] + df['Delta_Alpha']
df['CD_corrected'] = df['CD'] + df['Delta_CD']
df['CM25c_corrected'] = df['CM25c'] - df['Delta_Cm']  # Note the subtraction!

# Check for missing matches (Tail-off CL_w)
missing_matches = df['CL_w'].isna().sum()
if missing_matches > 0:
    print(f"WARNING: {missing_matches} rows could not find a matching CL_w in the tail-off data.")

# Check for missing slopes (C_L_alpha)
missing_slopes = df['CL_alpha_mapped'].isna().sum()
if missing_slopes > 0:
    print(f"WARNING: {missing_slopes} rows could not find a matching CL_alpha in the lookup table.")

# 9. Save back to the main file (or a new file)
output_filename = 'fully_corrected_data.csv'
df.to_csv(output_filename, index=False)

print(f"Data successfully matched and corrections applied! Saved to {output_filename}")