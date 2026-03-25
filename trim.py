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
                    (df['AoA'] <= 10)
                ].sort_values(by='AoA')

            if df_filtered.empty:
                print(f"  No data for de={delta_e}, J={J_target}")
                i += 1
                continue

            coeffs = np.polyfit(df_filtered['AoA'], df_filtered['CM25c'], 1)
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
            plt.plot(df_filtered['AoA'], df_filtered['CM25c'],
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

        # ── Longitudinal stability dCm/dAlpha vs J ────────────────────────────────
        plt.figure()
        stability_data = []

        for V_target, V_low, V_high in [(40, 38, 42), (20, 18, 22)]:
            for delta_e in [5, 0, -5]:
                for J_target in [1.6, 1.9, 2.2]:
                    df_filtered = df[
                        (df['de'] == delta_e) &
                        (df['V'] >= V_low) & (df['V'] <= V_high) &
                        (df['J'].between(J_target - 0.05, J_target + 0.05)) &
                        (df['AoA'] <= 10)
                    ].sort_values(by='AoA')

                    if df_filtered.empty:
                        print(f"  No data for V={V_target}, de={delta_e}, J={J_target}")
                        continue

                    coeffs = np.polyfit(df_filtered['AoA'], df_filtered['CM25c'], 1)
                    dCm_dAlpha = coeffs[0]
                    stability_data.append((V_target, delta_e, J_target, dCm_dAlpha))

        color_map  = {40: 'navy', 20: 'lightsteelblue'}
        marker_map = {40: 'o',         20: 's'}
        style_map  = {5: '--',         0: '-',   -5: '-.'}

        # One legend entry per V (color) and one per delta_e (linestyle)
        handles, labels = [], []
        for V_target, marker, color in [(40, 'o', 'navy'), (20, 's', 'lightsteelblue')]:
            for delta_e, linestyle in [(5, '--'), (0, '-'), (-5, '-.')]:
                V_de_data = [
                    (item[2], item[3]) for item in stability_data
                    if item[0] == V_target and item[1] == delta_e
                ]
                if not V_de_data:
                    print(f"No data to plot for V={V_target}, de={delta_e}")
                    continue

                J_vals, dCm_vals = zip(*sorted(V_de_data))
                line, = plt.plot(J_vals, dCm_vals,
                                marker=marker, color=color, linestyle=linestyle,
                                linewidth=1.5, markersize=7)
                handles.append(line)
                labels.append(f'V = {V_target} m/s, $\\delta_e$ = {delta_e}°')

        plt.axhline(0, color='k', linestyle='-', linewidth=0.5)
        plt.xlabel('Advance Ratio $J$')
        plt.ylabel(r'$\partial C_{m_{25c}} / \partial \alpha$ (1/°)')
        plt.ylim(-0.05, -0.03)
        plt.legend(handles, labels)
        plt.grid()
        plt.show()

        """# ── Trim angle of attack (Cm = 0) for every V, delta_e, J ────────────────
        print("\n" + "=" * 65)
        print("TRIM ANGLE OF ATTACK  (C_m = 0,  linear fit on AoA ≤ 10°)")
        print("=" * 65)
    
        trim_rows = []
    
        for V_target, V_low, V_high in [(40, 38, 42), (20, 18, 22)]:
            for delta_e in [5, 0, -5]:
                for J_target in [1.6, 1.9, 2.2]:
                    df_filtered = df[
                        (df['de'] == delta_e) &
                        (df['V'] >= V_low) & (df['V'] <= V_high) &
                        (df['J'].between(J_target - 0.05, J_target + 0.05)) &
                        (df['AoA'] <= 10)
                    ].sort_values(by='AoA')
    
                    if df_filtered.empty:
                        print(f"  No data for V={V_target}, de={delta_e}, J={J_target}  — skipped")
                        continue
    
                    # Linear fit:  Cm = m * AoA + b  →  AoA_trim = -b / m
                    m, b = np.polyfit(df_filtered['AoA'], df_filtered['CM25c'], 1)
    
                    if abs(m) < 1e-10:
                        aoa_trim = None
                        print(f"  V={V_target:2d}, de={delta_e:+3d}°, J={J_target}: "
                            f"slope ≈ 0 — no trim solution")
                    else:
                        aoa_trim = -b / m
    
                    trim_rows.append({
                        'V (m/s)'       : V_target,
                        'delta_e (°)'   : delta_e,
                        'J'             : J_target,
                        'slope dCm/dα'  : round(m, 6),
                        'AoA_trim (°)'  : round(aoa_trim, 3) if aoa_trim is not None else None,
                    })
    
        # Pretty-print as a table
        df_trim = pd.DataFrame(trim_rows)
        pd.set_option('display.float_format', '{:.3f}'.format)
        pd.set_option('display.max_rows', 50)
        print(df_trim.to_string(index=False))
    
        # Save to CSV
        trim_csv = 'trim_aoa_results.csv'
        df_trim.to_csv(trim_csv, index=False)
        print(f"\nTrim results saved to '{trim_csv}'")"""
