#%% Case study
import numpy as np
import matplotlib.pyplot as plt


# parameters
m = 100   # 1 - b = m * b
total = 0.2 * (1 + m)
a = 0.15 / total  # radius of conductor, cm
b = 0.2 / total  # conduction + insulator
L = 1   # conductor + insulator + dust layer
delta = (m + 1) ** 2 * 2e-4  # heat production rate
Q = 4  # from Joule heating

# thermal conductivity in unit W / (cm * K)
k0 = 3.85
k1 = 0.002
k2 = 0.0015

h = 0.0026  # in newton cooling

N_r = 501  # number of grids + 1
dr = L / (N_r - 1)

dt = 0.8 * (dr ** 2) / 2  # to satisfy CFL condition for explicit method
num_steps = 330000 # 1000000

snapshot_interval = 25000 # 50000
snapshots = []


# Initialization and classification
r = np.linspace(0, L, N_r)
region = np.zeros(N_r, dtype=int)

a_idx = np.argmin(np.abs(r - a))  # approximate the position of a
b_idx = np.argmin(np.abs(r - b))
# a_idx = int(a / L * (N_r - 1))
# b_idx = int(b / L * (N_r - 1))

region[:a_idx] = 0  # (0 ≤ r < a)
region[a_idx:b_idx] = 1  # (a ≤ r < b)
region[b_idx:] = 2  # (b ≤ r ≤ L)

theta = np.zeros(N_r)  # initial (dimensionless) temperature

# to plot u vs t at r = b and r = L
time_history = []
temp_b_history = []
term_1_history = []


# FDM
for n in range(num_steps):
    theta_new = np.copy(theta)

    # laplacian term for r = 0
    if np.any(r < 1e-6):
        theta_new[0] += dt * (4 * (theta[1] - theta[0]) / dr ** 2 + Q)

    # interior (exclude r = 0 and boundary points)
    d2T = (theta[2:] - 2 * theta[1:-1] + theta[:-2]) / dr ** 2
    dTdr = (theta[2:] - theta[:-2]) / (2 * dr)
    laplacian = d2T + (1 / r[1:-1]) * dTdr

    # heat diffusion term
    diff_ratio = np.array([1.0, 0.00116, 0.00076])
    diff = diff_ratio[region]

    # heat source term
    source = np.zeros_like(theta)
    source[region == 0] = Q  # conductor
    source[region == 2] = delta * np.exp(theta[region == 2])  # dust layer

    # update the temperature. now we have u for all r except r = 1
    theta_new[1:-1] += dt * (diff[1:-1] * laplacian + source[1:-1])
    # theta_new += dt * (laplacian + source)

    if np.any(~np.isfinite(theta_new)):
        nan_mask = np.isnan(theta_new)
        inf_mask = np.isinf(theta_new)
        print(f"u diverges from the time step {n} (t = {n * dt:.4f})")
        print(f"max u: {np.nanmax(theta_new):.2e}")
        break

    # (boundary) conditions
    # symmetric (r=0). used in deriving the laplacian term for r = 0

    # boundary condition (r=L). Newton cooling with u_a = 0
    theta_new[-1] = (k2 * theta_new[-2]) / (k2 + h * dr)

    # temperature continuity condition.
    # u(x-) = u(x+) will be automatically satisfied in the plotting since we are using grids
    # (r=a)
    theta_new[a_idx] = (k0 * theta[a_idx - 1] + k1 * theta[a_idx + 1]) / (k0 + k1)

    # (r=b)
    theta_new[b_idx] = (k1 * theta[b_idx - 1] + k2 * theta[b_idx + 1]) / (k1 + k2)

    # u defined on every r (grid points)
    theta = np.copy(theta_new)

    time_history.append(n * dt)
    temp_b_history.append(theta[b_idx])
    term_1_history.append(theta[-1])

    # save the snapshots
    if n % snapshot_interval == 0:
    # if n % snapshot_interval == 0 or n == num_steps - 1:
        snapshots.append((n * dt, theta.copy()))

snapshots.append((num_steps * dt, theta.copy()))

# visualization. dpi: The resolution of the figure in dots-per-inch.
plt.figure(figsize=(12, 8), dpi=100)

# temperature throughout time
for t, data in snapshots:
    plt.plot(r, data, lw=1.5,
             label=f'τ = {t:.4f}')

# surface
plt.axvline(a, color='r', ls='--', lw=1)
plt.axvline(b, color='b', ls='--', lw=1)

# indicates each part
# plt.annotate('Conductor', (a / 2, np.max(theta) * 1.05), ha='center', fontsize=12)
# plt.annotate('Insulator', ((a + b) / 2, np.max(theta) * 1.05), ha='center', fontsize=12)
plt.annotate('Dust Layer', ((b + L) / 2, np.max(theta) * 1.05), ha='center', fontsize=12)

plt.xlabel('r', fontsize=14)
plt.ylabel('u', fontsize=14)
# plt.title('Radial Temperature Distribution Evolution', fontsize=16, pad=20)

plt.grid(True, alpha=0.3)
plt.legend(loc='upper right')
plt.xlim(0, L)
plt.ylim(0, np.max(theta) * 1.1)

plt.tight_layout()
plt.show()

# visualize the dynamic boundary condition.
plt.figure(figsize=(8, 8))
plt.plot(time_history, temp_b_history, lw=1.5, color='red',label='r = b')
plt.plot(time_history, term_1_history, lw=1.5, color='blue',label='r = 1')
plt.xlabel('τ', fontsize=14)
plt.ylabel('u', fontsize=14)
# plt.title('Temperature Evolution at r = b', fontsize=16)
plt.grid(True, alpha=0.3)
plt.xlim(0, num_steps * dt)
# plt.ylim(0, max(temp_b_history) * 1.1)
plt.legend(loc='upper right', fontsize=14)
plt.tight_layout()
plt.show()

