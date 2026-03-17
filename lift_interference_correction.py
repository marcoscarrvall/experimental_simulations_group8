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


# --- THE CORRECTION MATH ---

# 6. Define your constants (UPDATE C TO YOUR TUNNEL AREA!)
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

print(dCm_dalpha_t)
area_ratio = S / C
Cla = 1 #...........................................................................

CMuw = (1/8)*tau_2*(0.5*Cmac)*delta*area_ratio*df['CL_w']*Cla # Cmac???????????????????????
CMt = dCm_dalpha_t*delta*area_ratio*df['CL_w']*(1+tau_2t*lt)
# 7. Calculate the correction terms using the newly matched 'CL_w'
df['Delta_Alpha'] = delta * area_ratio * df['CL_w'] * (1+tau_2*(0.5*Cmac))*57.3 # NOT SURE IF IT SHOULD BE Cmac
df['Delta_CD'] = delta * area_ratio * (df['CL_w'] ** 2)
df['Delta_Cm'] = (CMuw+CMt) * 57.3

# 8. Apply the corrections to create new columns
df['AoA_corrected'] = df['AoA'] + df['Delta_Alpha']
df['CD_corrected'] = df['CD'] + df['Delta_CD']
df['CM25c_corrected'] = df['CM25c'] - df['Delta_Cm']  # Note the subtraction!

# Check for any missing matches and warn the user
missing_matches = df['CL_w'].isna().sum()
if missing_matches > 0:
    print(f"WARNING: {missing_matches} rows could not find a matching CL_w in the tail-off data.")

# 9. Save back to the main file (or a new file)
output_filename = 'fully_corrected_data.csv'
df.to_csv(output_filename, index=False)

print(f"Data successfully matched and corrections applied! Saved to {output_filename}")