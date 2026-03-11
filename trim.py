import numpy as np
import matplotlib.pyplot as plt
from read_data import process_data

if __name__ == "__main__":
    df = process_data()

    # Lift polar
    df_delta_e0_V40_J16 = df[(df['delta_e'] == -5) & (df['V'] >= 38) & (df['V'] <= 42) & (df['J'] == 1.6)]
    df_delta_e_0_V40_J16 = df_delta_e0_V40_J16.sort_values(by='Alpha')
    df_delta_e0_V20_J16 = df[(df['delta_e'] == -5) & (df['V'] >= 18) & (df['V'] <= 22) & (df['J'] == 1.9)]
    df_delta_e0_V20_J16 = df_delta_e0_V20_J16.sort_values(by='Alpha')
    plt.plot(df_delta_e_0_V40_J16['Alpha'], df_delta_e_0_V40_J16['CL'], marker='x', color = "k", label = "V = 40 m/s, J = 1.6")
    plt.plot(df_delta_e0_V20_J16['Alpha'], df_delta_e0_V20_J16['CL'], marker='o', color = "r", label = "V = 20 m/s, J = 1.9")
    plt.xlabel('Angle of Attack ($\\alpha$)')
    plt.ylabel('Lift Coefficient ($c_l$)')

    plt.grid()
    plt.legend()
    plt.show()
    

    """colors = plt.cm.Blues(np.linspace(0.5, 1, 9))
    i = 0
    coeffs_list = []
    for delta_e in [5, 0, -5]:
        for V in [40]:
            for J in [1.6, 1.9, 2.2]:
                df_filtered = df[(df['delta_e'] == delta_e) & (df['V'] >= V-2) & (df['V'] <= V+2) & (df['J'] == J)]
                df_filtered = df_filtered.sort_values(by='Alpha')
                coeffs = np.polyfit(df_filtered['Alpha'], df_filtered['Cm_pitch'], 1)
                coeffs_list.append((delta_e, V, J, coeffs))
                plt.plot(df_filtered['Alpha'], df_filtered['Cm_pitch'], marker='x', color = colors[i], label = f"delta_e={delta_e}, J={J}")
                i += 1
    
    plt.xlabel('Angle of Attack ($\\alpha$)')
    plt.ylabel('Pitching Moment Coefficient ($c_m$)')
    plt.axhline(0, color='k', linestyle='-', linewidth=0.5)
    plt.grid()
    plt.legend()
    plt.show()

    plt.figure()
    J_values = sorted(list(set([item[2] for item in coeffs_list])))
    colors = plt.cm.Blues(np.linspace(0.5, 1, len(J_values)))
    alphas = np.arange(-4, 9, 1)

    for target_J in J_values:
        print(f"Processing J = {target_J}")
        J_coeffs = [item for item in coeffs_list if item[2] == target_J]
        J_coeffs.sort(key=lambda x: x[0])
        de_trim_for_this_J = []
        for alpha in alphas:
            de_axis = []
            cm_axis = []
            for delta_e, V, J, coeffs in J_coeffs:
                cm = coeffs[0] * alpha + coeffs[1]
                de_axis.append(delta_e)
                cm_axis.append(cm)
            fit_coeffs = np.polyfit(de_axis, cm_axis, 2)
            roots = np.roots(fit_coeffs)
            valid_roots = [root for root in roots if np.isreal(root) and -20 <= root <= 20]

            if valid_roots:
                de_trim_for_this_J.append((alpha, valid_roots[0].real))
            else:
                de_trim_for_this_J.append((alpha, None))
            
        plt.plot(alphas, [item[1] for item in de_trim_for_this_J], marker='o', color=colors[J_values.index(target_J)], label=f"J = {target_J}")

    plt.xlabel("Angle of Attack ($\\alpha$)")
    plt.ylabel("Elevator Deflection for Trim Condition ($\\delta_e$)")

    plt.grid()
    plt.legend()
    plt.show()"""
