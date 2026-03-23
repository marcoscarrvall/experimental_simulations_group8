import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from read_data import process_data

if __name__ == "__main__":
    df = process_data(csv_file='fully_corrected_data.csv')

    """# ── Lift polar ────────────────────────────────────────────────────────────
    df_V40_J16 = df[(df['de'] == 0) & (df['V'] >= 38) & (df['V'] <= 42) & (df['J'].between(1.59, 1.61))]
    df_V40_J16 = df_V40_J16.sort_values(by='AoA')

    df_V20_J19 = df[(df['de'] == 0) & (df['V'] >= 18) & (df['V'] <= 22) & (df['J'].between(1.89, 1.91))]
    df_V20_J19 = df_V20_J19.sort_values(by='AoA')

    plt.figure()
    plt.plot(df_V40_J16['AoA_corrected'], df_V40_J16['CL'], marker='x', color='k', label='V = 40 m/s, J ≈ 1.6')
    plt.plot(df_V20_J19['AoA_corrected'], df_V20_J19['CL'], marker='o', color='r', label='V = 20 m/s, J ≈ 1.9')
    plt.xlabel('Angle of Attack (°)')
    plt.ylabel('Lift Coefficient ($C_L$)')
    plt.grid()
    plt.legend()
    plt.title('Lift Polar')
    plt.show()"""

    # ── Pitching moment vs AoA ────────────────────────────────────────────────
    plt.figure()
    colors = plt.cm.Blues(np.linspace(0.5, 1, 3))
    i = 0
    coeffs_list = []
    
    for delta_e in [5, 0, -5]:
        for J_target in [1.6, 1.9, 2.2]:
            df_filtered = df[
                    (df['de'] == delta_e) &
                    (df['V'] >= 38) & (df['V'] <= 42) &
                    (df['J'].between(J_target - 0.05, J_target + 0.05)) &
                    (df['AoA_corrected'] <= 10)
                ].sort_values(by='AoA_corrected')

            if df_filtered.empty:
                print(f"  No data for de={delta_e}, J={J_target}")
                i += 1
                continue

            coeffs = np.polyfit(df_filtered['AoA_corrected'], df_filtered['CM25c_corrected'], 1)
            coeffs_list.append((delta_e, 40, J_target, coeffs))
            
            if J_target == 1.6: 
                linestyle = '--'
            elif J_target == 1.9:
                linestyle = '-.'
            elif J_target == 2.2:
                linestyle = '-'
            else:
                print(f"Unexpected J value: {J_target}, using dotted line.")
                linestyle = ':'
            plt.plot(df_filtered['AoA_corrected'], df_filtered['CM25c_corrected'],
                     marker='x', color=colors[i], label=f'$\\delta_e$={delta_e}°, J={J_target}', linestyle=linestyle)
        i += 1

    plt.xlabel('Angle of Attack (°)')
    plt.ylabel('Pitching Moment Coefficient ($C_{m_{25c}}$)')
    plt.axhline(0, color='k', linestyle='-', linewidth=0.5)
    plt.grid()
    plt.legend()
    plt.show()


    # ── Trim elevator deflection ───────────────────────────────────────────────
    if coeffs_list:
        plt.figure()
        J_values = sorted(set(item[2] for item in coeffs_list))
        colors = plt.cm.Blues(np.linspace(0.5, 1, len(J_values)))
        alphas = np.arange(-4, 9, 1)

        for j_idx, target_J in enumerate(J_values):
            J_coeffs = sorted(
                [item for item in coeffs_list if item[2] == target_J],
                key=lambda x: x[0]   # sort by delta_e
            )

            if len(J_coeffs) < 2:
                print(f"Not enough data points for J = {target_J}, skipping trim plot.")
                continue

            de_trim_for_this_J = []
            for alpha in alphas:
                de_axis = [item[0] for item in J_coeffs]
                cm_axis = [item[3][0] * alpha + item[3][1] for item in J_coeffs]

                # Linear fit if only 2 points, quadratic if 3+
                poly_order = min(2, len(de_axis) - 1)
                fit_coeffs = np.polyfit(de_axis, cm_axis, poly_order)
                roots = np.roots(fit_coeffs)
                valid_roots = [r.real for r in roots if np.isreal(r) and -20 <= r.real <= 20]

                de_trim_for_this_J.append(valid_roots[0] if valid_roots else None)

            plt.plot(alphas, de_trim_for_this_J, marker='o',
                     color=colors[j_idx], label=f'J = {target_J}')

        plt.xlabel('Angle of Attack (°)')
        plt.ylabel('Trim Elevator Deflection $\\delta_e$ (°)')
        plt.grid()
        plt.legend()
        plt.show()
