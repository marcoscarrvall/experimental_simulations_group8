"""
Wind Tunnel Blockage Correction Script

Calculates solid blockage, wake blockage, and slipstream blockage
for each row in the processed wind tunnel CSV data.
"""

import numpy as np
import pandas as pd

def blockage_corrections(filename="ALLconfigs_processed_sortedByRun.csv"):
    
    # ─────────────────────────────────────────────
    # 1. TUNNEL & MODEL CONSTANTS
    # ─────────────────────────────────────────────

    # Tunnel geometry
    B = 1.80        # Tunnel breadth [m]
    H = 1.25        # Tunnel height [m]
    C = B * H - 4 * (0.3**2)   # Effective cross-sectional area [m²] = 1.89 m²

    # Shape / tunnel factors (from report Tables 1–2 and BRP charts)
    tau_f = 0.862;  K3_f = 0.90   # Fuselage  (K3 read from fig. 10.2, d/l=0.09375)
    tau_n = 0.862;  K3_n = 1.00   # Nacelle   (K3, d/l=0.165)
    tau_w = 0.882;  K1_w = 1.02   # Wing      (K1, NACA-64 series, t/c=0.15)
    tau_h = 0.860;  K1_h = 1.02   # H-tail    (K1, NACA-64 series, t/c=0.15)
    tau_v = 0.861;  K1_v = 0.86   # V-tail    (K1, 4-digit series, t/c=0.15)

    # Component volumes [m³]  (Table 3 of report)
    V_f = 0.0160632
    V_n = 0.0007921
    V_w = 0.0030229
    V_h = 0.0009751
    V_v = 0.0003546

    # Wing reference area (used for wake blockage)
    S_ref = 0.2172           # [m²]  (standard Fokker F27 scaled model value)

    # Propeller geometry
    D_p  = 0.2032           # Propeller diameter [m]
    S_p  = np.pi / 4 * D_p**2   # Propeller disk area [m²]

    # Air density (ISA sea-level)
    rho  = 1.225            # [kg/m³]

    # CD–CL² drag polar slope  (K in CD = CD0 + K·CL²)
    # Realistic value for a clean low-speed transport aircraft model: ~0.055
    K_polar = 0.055

    # Zero-lift drag coefficient (y-intercept of the polar, prop-off)
    CD_0 = 0.025            # Realistic baseline for this type of model

    # ─────────────────────────────────────────────
    # 2. SOLID BLOCKAGE  (constant per component)
    # ─────────────────────────────────────────────
    # Fuselage / nacelle: eps_sb = K3 * tau1 * V / C^(3/2)
    # Wing / tail:        eps_sb = K1 * tau1 * V / C^(3/2)

    C32 = C**1.5   # C^(3/2)

    eps_sb_f = K3_f * tau_f * V_f / C32
    eps_sb_n = K3_n * tau_n * V_n / C32
    eps_sb_w = K1_w * tau_w * V_w / C32
    eps_sb_h = K1_h * tau_h * V_h / C32
    eps_sb_v = K1_v * tau_v * V_v / C32

    eps_sb_total = eps_sb_f + eps_sb_n + eps_sb_w + eps_sb_h + eps_sb_v

    print("=" * 60)
    print("CONSTANT BLOCKAGE COMPONENTS")
    print("=" * 60)
    print(f"  Tunnel cross-section C          = {C:.4f} m²")
    print(f"  Solid blockage – Fuselage       = {eps_sb_f:.6f}")
    print(f"  Solid blockage – Nacelle        = {eps_sb_n:.6f}")
    print(f"  Solid blockage – Wing           = {eps_sb_w:.6f}")
    print(f"  Solid blockage – H-tail         = {eps_sb_h:.6f}")
    print(f"  Solid blockage – V-tail         = {eps_sb_v:.6f}")
    print(f"  TOTAL solid blockage eps_sb     = {eps_sb_total:.6f}")
    print()

    # ─────────────────────────────────────────────
    # 3. HELPER FUNCTIONS (row-wise)
    # ─────────────────────────────────────────────

    def thrust_coefficient_CT(J):
        """
        CT from polynomial fit (given in report):
        CT = -0.0051*J^4 + 0.0959*J^3 - 0.5888*J^2 + 1.0065*J - 0.1353
        """
        return (-0.0051*J**4 + 0.0959*J**3 - 0.5888*J**2 + 1.0065*J - 0.1353)


    def compute_blockages(row):
        V_unc  = row["V"]       # uncorrected freestream velocity [m/s]
        CL     = row["CL"]
        CD_unc = row["CD"]
        J      = row["J"]
        n      = row["rpsM2"]   # rotational speed [rev/s]  (labelled rpsM2 in CSV)

        # --- Wake blockage (attached flow) ---
        # eps_wb0 = (S_ref / (4*C)) * CD_0
        eps_wb0 = (S_ref / (4 * C)) * CD_0

        # --- Lift-induced drag ---
        CD_i = K_polar * CL**2

        # --- Separated-flow wake blockage ---
        # CD_s = CD_unc - CD_0 - CD_i  (clamped to 0 so we don't go negative)
        CD_s    = max(CD_unc - CD_0 - CD_i, 0.0)
        eps_wb_s = (5 * S_ref / (4 * C)) * CD_s

        eps_wb_total = eps_wb0 + eps_wb_s

        # --- Slipstream blockage ---
        CT = thrust_coefficient_CT(J)
        # Thrust: T = CT * rho * n^2 * D^4
        T  = CT * rho * n**2 * D_p**4
        # Specific thrust coefficient: T*_c = T / (rho * V² * S_p)
        T_star_c = T / (rho * V_unc**2 * S_p) if V_unc > 0 else 0.0

        # eps_ss = - (T*_c / 2) / sqrt(1 + 2*T*_c) * (S_p / C)
        denom    = np.sqrt(max(1 + 2 * T_star_c, 1e-12))
        eps_ss   = -(T_star_c / (2 * denom)) * (S_p / C)

        # --- Total blockage ---
        eps_total = eps_sb_total + eps_wb_total + eps_ss

        # --- Corrections ---
        V_cor   = V_unc   * (1 + eps_total)
        q_unc   = 0.5 * rho * V_unc**2
        q_cor   = q_unc   * (1 + 2 * eps_total)
        delta_cd_wb = eps_sb_total * CD_unc   # drag correction additive term

        return pd.Series({
            "eps_sb"        : eps_sb_total,
            "eps_wb0"       : eps_wb0,
            "eps_wb_s"      : eps_wb_s,
            "eps_wb"        : eps_wb_total,
            "T_star_c"      : T_star_c,
            "eps_ss"        : eps_ss,
            "eps_total"     : eps_total,
            "V_cor"         : V_cor,
            "q_unc"         : q_unc,
            "q_cor"         : q_cor,
            "delta_cd_wb"   : delta_cd_wb,
        })

    # ─────────────────────────────────────────────
    # 4. APPLY TO ALL ROWS
    # ─────────────────────────────────────────────

    df = pd.read_csv(filename)

    blockage_cols = df.apply(compute_blockages, axis=1)
    df_out = pd.concat([df, blockage_cols], axis=1)

    # ─────────────────────────────────────────────
    # 5. CORRECTED AERODYNAMIC COEFFICIENTS
    # ─────────────────────────────────────────────
    # Scale coefficients down by (1 + 2*eps) and correct drag
    df_out["CL_cor"]  = df_out["CL"]  / (1 + 2 * df_out["eps_total"])
    df_out["CD_cor"]  = (df_out["CD"] + df_out["delta_cd_wb"]) / (1 + 2 * df_out["eps_total"])
    df_out["CM_cor"]  = df_out["CM25c"] / (1 + 2 * df_out["eps_total"])

    # ─────────────────────────────────────────────
    # 6. SAVE
    # ─────────────────────────────────────────────

    output_path = "blockage_corrected_data.csv"
    df_out.to_csv(output_path, index=False)

    print(f"Output saved to: {output_path}")