import pandas as pd
import numpy as np
from scipy import stats
from scipy.interpolate import interp1d

# ─────────────────────────────────────────────────────────────────────────────
# 1. FILENAMES
# ─────────────────────────────────────────────────────────────────────────────

propoff_raw_file = 'propOff_processed_byRun.csv'
main_file        = 'blockage_corrected_data_propoff.csv'
tail_off_file    = 'processed_tailOff_beta0.csv'
output_file      = 'fully_corrected_data_propoff.csv'

# ─────────────────────────────────────────────────────────────────────────────
# 2. COMPUTE PROP-OFF CL_ALPHA FROM LINEAR FIT
#    Use AoA ≈ -4°, 0°, 8°  (same three points used in the polar analysis)
#    grouped by nominal velocity (20 or 40 m/s).
#    Only use de=0, dr=0, AoS=0 rows to stay in the clean symmetric condition.
# ─────────────────────────────────────────────────────────────────────────────

df_raw = pd.read_csv(propoff_raw_file)

# Filter: de=0, dr=0, AoS≈0
df_clean = df_raw[
    (df_raw['de'].abs() < 0.1) &
    (df_raw['dr'].abs() < 0.1) &
    (df_raw['AoS'].abs() < 0.1)
].copy()

# Assign nominal velocity bucket (20 or 40 m/s)
df_clean['V_nom'] = (df_clean['V'] / 10).round().astype(int) * 10

# Explicitly map anything around -5° or -4° to the -4° bucket for the slope fit
def nearest_aoa(aoa):
    if -6.5 <= aoa <= -2.5:   # Catches all -5 and -4 data points
        return -4
    elif -1.5 <= aoa <= 1.5:  # Catches all 0 data points
        return 0
    elif 6.5 <= aoa <= 9.5:   # Catches all 8 data points
        return 8
    return np.nan

df_clean['AoA_nom'] = df_clean['AoA'].apply(nearest_aoa)
df_clean = df_clean.dropna(subset=['AoA_nom'])

# Average CL for each (V_nom, AoA_nom) combination
df_pts = df_clean.groupby(['V_nom', 'AoA_nom'])['CL'].mean().reset_index()

# Linear regression CL vs AoA per velocity → slope = CL_alpha [1/deg]
cla_propoff = {}
print("=" * 55)
print("PROP-OFF CL_ALPHA (linear fit on AoA = -4, 0, 8°)")
print("=" * 55)
for v_nom, grp in df_pts.groupby('V_nom'):
    grp = grp.sort_values('AoA_nom')
    if len(grp) < 2:
        print(f"  V = {v_nom} m/s : not enough points — skipping")
        continue
    slope, intercept, r, p, se = stats.linregress(grp['AoA_nom'], grp['CL'])
    cla_propoff[v_nom] = slope
    print(f"  V = {v_nom:2d} m/s :  CL_alpha = {slope:.6f} /deg   "
          f"(R² = {r**2:.5f},  n = {len(grp)})")
print()

# Fall back to mean if a velocity bucket is missing
if cla_propoff:
    cla_fallback = np.mean(list(cla_propoff.values()))
else:
    cla_fallback = 0.1057   # hard-coded safety net
    print(f"WARNING: CL_alpha fit failed entirely; using fallback = {cla_fallback}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. LOAD BLOCKAGE-CORRECTED DATA & TAIL-OFF DATA
# ─────────────────────────────────────────────────────────────────────────────

df_main     = pd.read_csv(main_file)
df_tail_off = pd.read_csv(tail_off_file)

# Strip hidden whitespace from all column headers in both files
df_main.columns = df_main.columns.str.strip()
df_tail_off.columns = df_tail_off.columns.str.strip()

# ─────────────────────────────────────────────────────────────────────────────
# 4. MATCH TAIL-OFF CL_w (Using Continuous Interpolation)
# ─────────────────────────────────────────────────────────────────────────────

# Clean the tail-off data and assign nominal velocities
df_tail_off_clean = df_tail_off.dropna(subset=['CL', 'AoA', 'V']).copy()
df_tail_off_clean['V_nom'] = (df_tail_off_clean['V'] / 10).round().astype(int) * 10

# Create an interpolation function (curve fit) for each velocity bucket
clw_interpolators = {}
for v_nom, grp in df_tail_off_clean.groupby('V_nom'):
    # Sort and remove duplicate AoA points to ensure a clean interpolation curve
    grp = grp.sort_values('AoA').drop_duplicates(subset=['AoA'])
    if len(grp) > 1:
        # fill_value='extrapolate' safely estimates CL if the main dataset goes slightly past the tail-off AoA limits
        clw_interpolators[v_nom] = interp1d(grp['AoA'], grp['CL'], kind='linear', fill_value='extrapolate')

# Assign nominal velocity to the main dataset to know which curve to use
df_main['V_match'] = (df_main['V'] / 10).round().astype(int) * 10

# Safely extract the exact CL_w for the exact AoA
def get_clw(row):
    v = row['V_match']
    aoa = row['AoA'] # Uses the exact, unrounded AoA
    if v in clw_interpolators:
        return float(clw_interpolators[v](aoa))
    return np.nan

df_main['CL_w'] = df_main.apply(get_clw, axis=1)

# Transition df_main to 'df' for the rest of the script
df = df_main.copy()

# Map the fitted CL_alpha onto each row by its nominal velocity
df['CL_alpha_mapped'] = df['V_match'].map(cla_propoff).fillna(cla_fallback)

df.drop(columns=['V_match'], inplace=True)

# ─────────────────────────────────────────────────────────────────────────────
# 5. GEOMETRIC & TUNNEL CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

S      = 0.2172
St     = 0.0858
ARt    = 3.87
C      = 2.07
delta  = 0.1033
tau_2  = 0.120
tau_2t = 0.75
Cmac   = 0.165
lt     = 0.535

dCm_dalpha_t = -((St * lt) / (S * Cmac)) * ((0.1 * ARt * 0.8) / (ARt + 2))
area_ratio   = S / C

# ─────────────────────────────────────────────────────────────────────────────
# 6. CORRECTION TERMS
# ─────────────────────────────────────────────────────────────────────────────

Cla  = df['CL_alpha_mapped']

CMuw = (1/8) * tau_2 * (0.5 * Cmac) * delta * area_ratio * df['CL_w'] * Cla
CMt  = dCm_dalpha_t * delta * area_ratio * df['CL_w'] * (1 + tau_2t * lt)

df['Delta_Alpha'] = delta * area_ratio * df['CL_w'] * (1 + tau_2 * (0.5 * Cmac)) * 57.3
df['Delta_CD']    = delta * area_ratio * (df['CL_w'] ** 2)
df['Delta_Cm']    = CMuw + CMt

# ─────────────────────────────────────────────────────────────────────────────
# 7. APPLY CORRECTIONS
# ─────────────────────────────────────────────────────────────────────────────

df['AoA']   = df['AoA']   + df['Delta_Alpha']
df['CD']    = df['CD']    + df['Delta_CD']
df['CM25c'] = df['CM25c'] + df['Delta_Cm']

# ─────────────────────────────────────────────────────────────────────────────
# 8. WARNINGS
# ─────────────────────────────────────────────────────────────────────────────

missing_clw = df['CL_w'].isna().sum()
if missing_clw > 0:
    print(f"WARNING: {missing_clw} rows could not find a matching CL_w in the tail-off data.")

missing_cla = (df['CL_alpha_mapped'] == cla_fallback).sum()
if missing_cla > 0:
    print(f"NOTE: {missing_cla} rows used the fallback CL_alpha = {cla_fallback:.6f} "
          f"(no fit available for their velocity bucket).")

# ─────────────────────────────────────────────────────────────────────────────
# 8.5 FIX MISSING REYNOLDS NUMBER
# ─────────────────────────────────────────────────────────────────────────────

# Reconstruct the blank Re column using the tunnel's linear V to Re relationship
df['Re'] = df['V'] * 10965

# ─────────────────────────────────────────────────────────────────────────────
# 9. OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

desired_columns = [
    'config', 'de', 'dr', 'run', 'AoA', 'AoS', 'V', 'Re',
    'rpsM2', 'J', 'CL', 'CD', 'CM25c', 'CYaw',
    'eps_total', 'eps_sb', 'eps_wb0'
]

final_columns = [c for c in desired_columns if c in df.columns]
df_final = df[final_columns]

df_final.to_csv(output_file, index=False)
print(f"Corrections applied. Saved to '{output_file}'  ({len(df_final)} rows)")
print(f"  dCm/dalpha_t = {dCm_dalpha_t:.6f}")