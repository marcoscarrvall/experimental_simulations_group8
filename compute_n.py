import numpy as np

def compute_n(J, V, D):
    n = V / (D * J)
    return n

V = np.array([20, 40])
D = 0.2032
J = np.array([1.6, 1.9, 2.2])

# For each combination of V and J, compute n

for v in V:
    for j in J:
        n = compute_n(j, v, D)
        print(f"For V={v} and J={j}, n={n:.4f}")