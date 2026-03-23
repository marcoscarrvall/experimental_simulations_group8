import pandas as pd
import numpy as np
from scipy import stats

# ─────────────────────────────────────────────────────────────────────────────
# 1. FILENAMES
# ─────────────────────────────────────────────────────────────────────────────

propoff_raw_file = 'propOff_processed_byRun.csv'
main_file        = 'blockage_corrected_data_propoff.csv'
tail_off_file    = r'TAILOFF\TAILOFF/unc_tailOff_beta0_balance.txt'
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

# Assign nominal AoA bucket — keep only rows near -4, 0, or 8 degrees
aoa_targets = [-4, 0, 8]
def nearest_aoa(aoa):
    nearest = min(aoa_targets, key=lambda t: abs(t - aoa))
    if abs(aoa - nearest) <= 1.5:   # tolerance of ±1.5°
        return nearest
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
df_tail_off = pd.read_csv(tail_off_file, sep=r'\s+', skiprows=[1])

# ─────────────────────────────────────────────────────────────────────────────
# 4. MATCH & MERGE TAIL-OFF CL_w
# ─────────────────────────────────────────────────────────────────────────────

df_main['AoA_match']     = df_main['AoA'].round().astype(int)
df_tail_off['AoA_match'] = df_tail_off['Alpha'].round().astype(int)

df_main['V_match']     = (df_main['V'] / 10).round().astype(int) * 10
df_tail_off['V_match'] = (df_tail_off['V'] / 10).round().astype(int) * 10

df_tail_off_subset = df_tail_off[['AoA_match', 'V_match', 'CL']].copy()
df_tail_off_subset.rename(columns={'CL': 'CL_w'}, inplace=True)

df = pd.merge(df_main, df_tail_off_subset, on=['AoA_match', 'V_match'], how='left')

# Map the fitted CL_alpha onto each row by its nominal velocity
df['V_match'] = (df['V'] / 10).round().astype(int) * 10
df['CL_alpha_mapped'] = df['V_match'].map(cla_propoff).fillna(cla_fallback)

df.drop(columns=['AoA_match', 'V_match'], inplace=True)

# ─────────────────────────────────────────────────────────────────────────────
# 5. GEOMETRIC & TUNNEL CONSTANTS  (identical to original script)
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