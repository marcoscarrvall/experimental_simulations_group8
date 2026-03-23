import numpy as np
import pandas as pd

def blockage_corrections_propoff(
    filename="propOff_processed_byRun.csv",
    output_path="blockage_corrected_data_propoff.csv"
):
    """
    Applies wind tunnel blockage corrections to propOff data.
    Since the propeller is off:
      - Slipstream blockage (eps_ss) = 0
      - Separated-flow wake blockage (eps_wb_s) = 0  (condition J > 2 never met)
    Only solid blockage and attached-flow wake blockage are applied.
    """

    # ─────────────────────────────────────────────
    # 1. TUNNEL & MODEL CONSTANTS  (identical to main script)
    # ─────────────────────────────────────────────

    B = 1.80
    H = 1.25
    C = B * H - 4 * (0.3**2)   # Effective cross-sectional area [m²]

    tau_f = 0.862;  K3_f = 0.905
    tau_n = 0.862;  K3_n = 0.935
    tau_w = 0.882;  K1_w = 1.02
    tau_h = 0.860;  K1_h = 1.02
    tau_v = 0.861;  K1_v = 1.035
    tau_SSaft  = 0.861; K1_SSaft  = 1.02
    tau_SSwing = 0.861; K1_SSwing = 1.02

    V_f      = 0.0160632
    V_n      = 0.0007921
    V_w      = 0.0030229
    V_h      = 0.0009751
    V_v      = 0.0003546
    V_SSaft  = 0.0004491
    V_SSwing = 0.0017648

    S_ref    = 0.2172       # Wing reference area [m²]
    rho      = 1.225        # Air density [kg/m³]
    K_polar  = 0.035        # Drag polar slope (CD = CD0 + K*CL²)

    # Zero-lift drag coefficient (speed-dependent)
    CD_0 = {20: 0.0581, 40: 0.05124}

    # ─────────────────────────────────────────────
    # 2. SOLID BLOCKAGE  (constant)
    # ─────────────────────────────────────────────

    C32 = C**1.5

    eps_sb_f      = K3_f  * tau_f  * V_f      / C32
    eps_sb_n      = K3_n  * tau_n  * V_n      / C32
    eps_sb_w      = K1_w  * tau_w  * V_w      / C32
    eps_sb_h      = K1_h  * tau_h  * V_h      / C32
    eps_sb_v      = K1_v  * tau_v  * V_v      / C32
    eps_sb_SSaft  = K1_SSaft  * tau_SSaft  * V_SSaft  / C32
    eps_sb_SSwing = K1_SSwing * tau_SSwing * V_SSwing / C32

    eps_sb_total = (eps_sb_f + 2*eps_sb_n + eps_sb_w +
                    eps_sb_h + eps_sb_v +
                    eps_sb_SSaft + 2*eps_sb_SSwing)

    # ─────────────────────────────────────────────
    # 3. ROW-WISE CORRECTION  (prop-off: eps_ss = eps_wb_s = 0)
    # ─────────────────────────────────────────────

    def compute_blockages(row):
        V_unc  = row["V"]
        CL     = row["CL"]
        CD_unc = row["CD"]

        # Choose CD_0 based on velocity
        CD_0_value = CD_0[20] if V_unc <= 30 else CD_0[40]

        # Attached-flow wake blockage
        eps_wb0 = (S_ref / (4 * C)) * CD_0_value

        # No separated-flow wake blockage (prop off → no J condition met)
        eps_wb_s   = 0.0
        eps_wb_total = eps_wb0

        # No slipstream blockage
        eps_ss = 0.0

        # Total blockage
        eps_total = eps_sb_total + eps_wb_total + eps_ss

        # Corrected quantities
        V_cor       = V_unc * (1 + eps_total)
        q_unc       = 0.5 * rho * V_unc**2
        q_cor       = q_unc * (1 + eps_total)**2
        delta_cd_wb = eps_sb_total * 0.055

        return pd.Series({
            "eps_sb"      : eps_sb_total,
            "eps_wb0"     : eps_wb0,
            "eps_wb_s"    : eps_wb_s,
            "eps_wb"      : eps_wb_total,
            "eps_ss"      : eps_ss,
            "eps_total"   : eps_total,
            "V_cor"       : V_cor,
            "q_unc"       : q_unc,
            "q_cor"       : q_cor,
            "delta_cd_wb" : delta_cd_wb,
        })

    # ─────────────────────────────────────────────
    # 4. LOAD & APPLY
    # ─────────────────────────────────────────────

    df = pd.read_csv(filename)

    blockage_cols = df.apply(compute_blockages, axis=1)
    df_out = pd.concat([df, blockage_cols], axis=1)

    # ─────────────────────────────────────────────
    # 5. CORRECTED AERODYNAMIC COEFFICIENTS
    # ─────────────────────────────────────────────

    df_out["V"]  = df_out["V_cor"]
    df_out["CL"] = df_out["CL"] / (1 + df_out["eps_total"])**2
    df_out["CD"] = (df_out["CD"] + df_out["delta_cd_wb"]) / (1 + df_out["eps_total"])**2
    df_out["CM"] = df_out["CM25c"] / (1 + df_out["eps_total"])**2

    # ─────────────────────────────────────────────
    # 6. PRINT SUMMARY
    # ─────────────────────────────────────────────

    print("=" * 60)
    print("CONSTANT BLOCKAGE COMPONENTS")
    print("=" * 60)
    print(f"  Tunnel cross-section C          = {C:.4f} m²")
    print(f"  Solid blockage – Fuselage       = {eps_sb_f:.6f}")
    print(f"  Solid blockage – Nacelle (x2)   = {eps_sb_n:.6f}  each")
    print(f"  Solid blockage – Wing           = {eps_sb_w:.6f}")
    print(f"  Solid blockage – H-tail         = {eps_sb_h:.6f}")
    print(f"  Solid blockage – V-tail         = {eps_sb_v:.6f}")
    print(f"  Solid blockage – SSaft          = {eps_sb_SSaft:.6f}")
    print(f"  Solid blockage – SWing (x2)     = {eps_sb_SSwing:.6f}  each")
    print(f"  TOTAL solid blockage eps_sb     = {eps_sb_total:.6f}")
    print()
    print(f"  Slipstream blockage eps_ss      = 0.000000  (prop off)")
    print(f"  Sep.-flow wake blockage eps_wb_s= 0.000000  (prop off)")
    print()

    # ─────────────────────────────────────────────
    # 7. SAVE
    # ─────────────────────────────────────────────

    # Keep the same columns as the propOn output where possible
    base_cols = ["config", "de", "dr", "run", "AoA", "AoS",
                 "V", "Re", "CL", "CD", "CM25c", "CYaw",
                 "eps_total", "eps_sb", "eps_wb0"]

    # Only include columns that actually exist in this file
    save_cols = [c for c in base_cols if c in df_out.columns]

    df_out[save_cols].to_csv(output_path, index=False)
    print(f"Output saved to: {output_path}  ({len(df_out)} rows)")

    return df_out


if __name__ == "__main__":
    blockage_corrections_propoff()