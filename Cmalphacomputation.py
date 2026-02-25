import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def process_aerodynamics():
    data_folder = "BAL"
    search_path = os.path.join(data_folder, "*.txt")
    files = sorted(glob.glob(search_path))
    
    files = [f for f in files if os.path.basename(f).startswith(('corr_', 'raw_', 'unc_'))]

    if not files:
        print(f"No matching .txt files found in '{data_folder}'.")
        return

    plt.figure(figsize=(12, 7))
    print(f"{'File Name':<32} | {'Cm_alpha [/deg]':<15} | {'Cm_alpha [/rad]':<15}")
    print("-" * 70)

    for file in files:
        try:
            df = pd.read_csv(file, sep='\s+', skiprows=[1])
            df.columns = df.columns.str.strip()
            
            alpha = pd.to_numeric(df['Alpha'], errors='coerce')
            cm = pd.to_numeric(df['Cm_pitch'], errors='coerce')
            valid_mask = alpha.notna() & cm.notna()
            alpha, cm = alpha[valid_mask], cm[valid_mask]

            if alpha.empty:
                continue

            # Linear Fit
            slope_deg, intercept = np.polyfit(alpha, cm, 1)
            slope_rad = slope_deg * (180 / np.pi)

            file_label = os.path.basename(file)
            print(f"{file_label:<32} | {slope_deg:>15.5f} | {slope_rad:>15.5f}")

            # Plotting raw points
            p = plt.scatter(alpha, cm, alpha=0.4, s=15) 
            color = p.get_facecolor()[0]
            
            # Plotting regression line
            x_fit = np.array([alpha.min(), alpha.max()])
            y_fit = slope_deg * x_fit + intercept
            # Simplified label to prevent LaTeX underscore errors
            plt.plot(x_fit, y_fit, color=color, linewidth=2, 
                     label=f"{file_label} (m={slope_deg:.4f})")
            
        except Exception as e:
            print(f"Error processing {os.path.basename(file)}: {e}")

    # Formatting with Raw Strings to fix your ValueError
    plt.axhline(0, color='black', linewidth=1, alpha=0.3)
    plt.axvline(0, color='black', linewidth=1, alpha=0.3)
    plt.xlabel(r'Angle of Attack $\alpha$ [deg]')
    plt.ylabel(r'Pitching Moment $C_m$ [-]')
    plt.title(r'Longitudinal Stability Analysis: $C_m$ vs $\alpha$')
    
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='8')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    process_aerodynamics()