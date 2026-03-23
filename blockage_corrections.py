import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def read_propOffData(V, dE, AoS = 0, dR = 0):
    """
    Reads propOff.xlsx, filters for V, AoS = 0°, dE, dR = 0°,
    saves a CSV with CL and CD, and returns the filtered DataFrame.
    """
    df_raw = pd.read_excel('propOff.xlsx', header=None)

    # Row 0 is the header; rows 1+ are data
    df_raw.columns = ['AoA', 'AoS', 'V', 'dE', 'dR', 'CL', 'CD', 'CYaw', 'CMroll', 'CMpitch', 'extra']
    df = df_raw.iloc[1:].apply(pd.to_numeric, errors='coerce').dropna(subset=['AoA', 'CL', 'CD'])

    # Filter: V ≈ 40 m/s (±1 m/s tolerance), AoS = 0°, dE = 0°, dR = 0°
    data = df[
        (df['V'].between(V-1, V+1)) &
        (df['AoS'].abs() < AoS+0.1)   &
        (df['dE'] == dE)            &
        (df['dR'] == dR)
    ][['AoA', 'AoS', 'V', 'CL', 'CD']].reset_index(drop=True)

    # Remove any repeats, ensure that all AoA values are unique (if not, take the mean of duplicates)
    data = data.groupby('AoA', as_index=False).mean()

    # Remove any negative CL values (if they exist) to avoid issues with the CD vs CL² plot
    data = data[data['CL'] >= 0].reset_index(drop=True)

    csv_path = f'propOff_zero_sideslip_V{int(V)}_{int(dE)}_{int(dR)}.csv'
    data.to_csv(csv_path, index=False)
    print(f"CSV saved: {csv_path}  ({len(data)} rows)")
    return data

def plot_CD_vs_CL2(V, dE, AoS = 0, dR = 0, color="k", marker = "x", linestyle = "-"):
    data = read_propOffData(V, dE, AoS, dR)
    CL2 = data['CL']**2
    CD  = data['CD']

    plt.plot(CL2, CD, color=color, zorder=5, marker = marker, linestyle = linestyle, label = f"V = {V} m/s, dE = {dE}°")

def blockage_corrections(filename="ALLconfigs_processed_sortedByRun.csv"):
    
    # ─────────────────────────────────────────────
    # 1. TUNNEL & MODEL CONSTANTS
    # ─────────────────────────────────────────────

    # Tunnel geometry
    B = 1.80        # Tunnel breadth [m]
    H = 1.25        # Tunnel height [m]
    C = B * H - 4 * (0.3**2)   # Effective cross-sectional area [m²] = 1.89 m²

    # Shape / tunnel factors (from report Tables 1–2 and BRP charts)
    tau_f = 0.862;  K3_f = 0.905   # Fuselage  (K3 read from fig. 10.2, d/l=0.09375)
    tau_n = 0.862;  K3_n = 0.935   # Nacelle   (K3, d/l=0.165)
    tau_w = 0.882;  K1_w = 1.02   # Wing      (K1, NACA-64 series, t/c=0.15)
    tau_h = 0.860;  K1_h = 1.02   # H-tail    (K1, NACA-64 series, t/c=0.15)
    tau_v = 0.861;  K1_v = 1.035   # V-tail    (K1, 4-digit series, t/c=0.15)
    tau_SSaft = 0.861; K1_SSaft = 1.02   # SSaft     (K1, 4-digit series, t/c=0.12)  
    tau_SSwing = 0.861; K1_SSwing = 1.02   # SWing     (K1, 4-digit series, t/c=0.12)


    # Component volumes [m³]  (Table 3 of report)
    V_f = 0.0160632
    V_n = 0.0007921
    V_w = 0.0030229
    V_h = 0.0009751
    V_v = 0.0003546
    V_SSaft = 0.0004491
    V_SSwing = 0.0017648

    # Wing reference area (used for wake blockage)
    S_ref = 0.2172           # [m²]  (standard Fokker F27 scaled model value)

    # Propeller geometry
    D_p  = 0.2032           # Propeller diameter [m]
    S_p  = np.pi / 4 * D_p**2   # Propeller disk area [m²]

    # Air density (ISA sea-level)
    rho  = 1.225            # [kg/m³]

    # CD–CL² drag polar slope  (K in CD = CD0 + K·CL²)
    # Realistic value for a clean low-speed transport aircraft model: ~0.055
    K_polar = 0.035

    # Zero-lift drag coefficient depends on V. For V = 20, C_D0 = 0.0581, for V = 40, C_D0 = 0.05124
    CD_0 = {20: 0.0581, 40: 0.05124}

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
    eps_sb_SSaft = K1_SSaft * tau_SSaft * V_SSaft / C32
    eps_sb_SSwing = K1_SSwing * tau_SSwing * V_SSwing / C32
    eps_sb_total = eps_sb_f + 2*eps_sb_n + eps_sb_w + eps_sb_h + eps_sb_v + eps_sb_SSaft + 2*eps_sb_SSwing

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
        if V_unc <= 30: 
            CD_0_value = CD_0[20]
        else:
            CD_0_value = CD_0[40]
        eps_wb0 = (S_ref / (4 * C)) * CD_0_value

        # --- Lift-induced drag ---
        CD_i = K_polar * CL**2

        # --- Separated-flow wake blockage ---
        # CD_s = CD_unc - CD_0 - CD_i  (clamped to 0 so we don't go negative)
        CD_s = CD_unc - CD_0_value - CD_i
        """if CD_s < 0:
            print(f"At CL = {CL} and V = {V_unc}: CD_s is negative")
        else:
            print(f"At AoA = {row['AoA']}, CL = {CL} and V = {V_unc}, J = {J}: {CD_s=:.4f} = {CD_unc:.4f} - {CD_0_value:.4f} - {CD_i:.4f}")
        """
        if V_unc > 30 and J > 2 and CD_s > 0.005:
            eps_wb_s = (5 * S_ref / (4 * C)) * CD_s
        else:
            eps_wb_s = 0.0

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
        q_cor   = q_unc   * (1 + eps_total)**2
        delta_cd_wb = eps_sb_total * 0.055   # drag correction additive term

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

    print("=" * 60)
    print("CONSTANT BLOCKAGE COMPONENTS")
    print("=" * 60)
    print(f"  Tunnel cross-section C          = {C:.4f} m²")
    print(f"  Solid blockage – Fuselage       = {eps_sb_f:.6f}")
    print(f"  Solid blockage – Nacelle        = {eps_sb_n:.6f}")
    print(f"  Solid blockage – Wing           = {eps_sb_w:.6f}")
    print(f"  Solid blockage – H-tail         = {eps_sb_h:.6f}")
    print(f"  Solid blockage – V-tail         = {eps_sb_v:.6f}")
    print(f"  Solid blockage – SSaft          = {eps_sb_SSaft:.6f}")
    print(f"  Solid blockage – SWing          = {eps_sb_SSwing:.6f}")
    print(f"  TOTAL solid blockage eps_sb     = {eps_sb_total:.6f}")
    print()

    """print("=" * 60)
    print("CORRECTION SUMMARY  (first 10 rows)")
    print("=" * 60)
    cols_show = ["run", "config", "AoA", "J", "eps_total", "V_cor", "CL", "CD", "CM25c"]
    pd.set_option("display.float_format", "{:.5f}".format)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 120)
    print(df_out[cols_show].head(10).to_string(index=False))"""

    # ─────────────────────────────────────────────
    # 5. CORRECTED AERODYNAMIC COEFFICIENTS
    # ─────────────────────────────────────────────
    # Scale coefficients down by (1 + eps)**2 and correct drag
    df_out["V"]   = df_out["V_cor"] 
    df_out["CL"]  = df_out["CL"]  / (1 + df_out["eps_total"])**2
    df_out["CD"]  = (df_out["CD"] + df_out["delta_cd_wb"]) / (1 + df_out["eps_total"])**2
    df_out["CM"]  = df_out["CM25c"] / (1 + df_out["eps_total"])**2

    # ─────────────────────────────────────────────
    # 6. SAVE
    # ─────────────────────────────────────────────

    original_columns = ["config", "de", "dr", "run", "AoA", "AoS",
                    "V", "Re", "rpsM2", "J", "CL", "CD", "CM25c", "CYaw", "eps_total", "eps_sb", "eps_wb0"]

    output_path = "blockage_corrected_data.csv"
    df_out[original_columns].to_csv(output_path, index=False)

    print(f"Output saved to: {output_path}")

    # ─────────────────────────────────────────────
    # 7. eps_wb_s vs CL
    # ─────────────────────────────────────────────
    # Use uncorrected CL/CD by re-reading the raw file
    df_wbs = df_out[["eps_wb_s", "AoA"]].reset_index(drop=True)
    # Get only non-zero wake blockage values for better plotting
    df_wbs = df_wbs[df_wbs["eps_wb_s"] > 0].reset_index(drop=True)


    wbs_path = "blockage_correction_datafiles/eps_wb_s_vs_AoA.csv"
    df_wbs.to_csv(wbs_path, index=False)
    print(f"Output saved to: {wbs_path}")

    # ─────────────────────────────────────────────
    # 8. eps_ss, V, J
    # ─────────────────────────────────────────────
    df_ss = df_out[["eps_ss", "V", "J"]].reset_index(drop=True)

    ss_path = "blockage_correction_datafiles/eps_ss_V_J.csv"
    df_ss.to_csv(ss_path, index=False)
    print(f"Output saved to: {ss_path}")

def blockage_plots():
    # ── File 1: eps_ss vs V and J ───────────────────────────────────────────────
    df1 = pd.read_csv('blockage_correction_datafiles/eps_ss_V_J.csv')
    df1['V'] = df1['V'].apply(lambda x: 20 if abs(x - 20) < abs(x - 40) else 40)
    targets = [1.6, 1.9, 2.2]
    df1['J'] = df1['J'].apply(lambda x: min(targets, key=lambda t: abs(t - x)))
    df1_clean = df1.groupby(['V', 'J'], as_index=False)['eps_ss'].mean()

    # ── File 2: eps_wb_s vs AoA ───────────────────────────────────────────────────
    df2 = pd.read_csv("blockage_correction_datafiles/eps_wb_s_vs_AoA.csv")
    df2_sorted = df2.sort_values("AoA")
    
    # ── Plotting ─────────────────────────────────────────────────────────────────    
    # — Plot 1 —
    linestyle = ["-", "--"]
    for v_val in sorted(df1_clean['V'].unique()):
        subset = df1_clean[df1_clean['V'] == v_val].sort_values('J')
        plt.plot(subset['J'], subset['eps_ss'], marker="x", linestyle=linestyle[df1_clean['V'].unique().tolist().index(v_val)], label=f'V = {v_val}', color = "k")
    
    plt.xlabel("J", fontsize=11)
    plt.ylabel("$\\epsilon_{ss}$", fontsize=11)
    plt.xlim(1.5, 2.3)
    plt.legend(title="Velocity (V)", fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig("corrections_images/epsilon_ss_plots.png", dpi=150, bbox_inches="tight")
    print("Plots saved to corrections_images/epsilon_ss_plots.png")
    plt.show()
    
    # — Plot 2 —
    plt.figure()
    # Overlay a sorted mean line (bin by rounded AoA for a trend guide)
    df2["AoA_round"] = df2["AoA"].round(2)
    agg2 = df2.groupby("AoA_round")["eps_wb_s"].mean().reset_index().sort_values("AoA_round")
    plt.plot(agg2["AoA_round"], agg2["eps_wb_s"],
            color="#000000", linewidth=2, zorder=3, label="Mean trend", marker = "x")
    
    plt.xlabel("$\\alpha$", fontsize=11)
    plt.ylabel("$\\epsilon_{wb,s}$", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.4)
    
    plt.tight_layout()
    plt.savefig("corrections_images/epsilon_wb_s_plots.png", dpi=150, bbox_inches="tight")
    print("Plots saved to corrections_images/epsilon_wb_s_plots.png")
    plt.show()

if __name__ == "__main__":
    blockage_corrections()
    blockage_plots()


    # Find K and CD0 from the raw propOff data (for wake blockage calculations)
    """V_array = [20, 40]
    deltaE_array = [-5, 0, 5]
    linestyle = ["-", "--"]
    markers = ["^", "o", "s"]
    for deltaE in deltaE_array:
        for V in V_array:
                plot_CD_vs_CL2(V, deltaE, linestyle = linestyle[V_array.index(V)], marker = markers[deltaE_array.index(deltaE)])

    plt.xlabel(r'$C_L^2$', fontsize=13)
    plt.ylabel(r'$C_D$',   fontsize=13)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(f'CD_vs_CL2_propOff.png', dpi=150)
    print(f"Plot saved: CD_vs_CL2_propOff.png")
    plt.show()"""

    
