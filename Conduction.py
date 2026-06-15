import numpy as np
import matplotlib.pyplot as plt

# material properties for steel and mineral wool
k_steel   = 50.0;    rho_steel = 7800.0;  c_steel = 480.0
k_wool    = 0.043;   rho_wool  = 10.0;    c_wool  = 837.0

alpha_steel = k_steel / (rho_steel * c_steel)
alpha_wool  = k_wool  / (rho_wool  * c_wool)

# geometry: steel | wool | steel sandwich, with hot side and convection on cold side
L_steel1  = 0.0012
L_wool    = 0.0005
L_steel2  = 0.0012
T_hot     = 500.0
T_ambient = 20.0
h         = 500.0


def build_x_array(N_steel1, N_wool, N_steel2):
    """Stitch together the three layer grids, dropping repeated endpoints."""
    x_s1 = np.linspace(0, L_steel1, N_steel1 + 1)
    x_w  = np.linspace(L_steel1, L_steel1 + L_wool, N_wool + 1)
    x_s2 = np.linspace(L_steel1 + L_wool, L_steel1 + L_wool + L_steel2, N_steel2 + 1)
    return np.concatenate([x_s1, x_w[1:], x_s2[1:]])


def build_alpha_array(N, end_steel1, end_wool):
    """Assign diffusivity layer by layer; interface nodes get the wool value."""
    alpha = np.empty(N)
    for i in range(N):
        if   i < end_steel1:   alpha[i] = alpha_steel
        elif i == end_steel1:  alpha[i] = alpha_wool
        elif i < end_wool:     alpha[i] = alpha_wool
        elif i == end_wool:    alpha[i] = alpha_steel
        else:                  alpha[i] = alpha_steel
    return alpha


def build_cn_matrices(N, end_steel1, end_wool, alpha, dx_steel1, dx_wool, dx_steel2, dt):
    """
    Crank-Nicolson system: A @ T_new = B @ T_old + b_bc

    Left face  — Dirichlet (T_hot)
    Right face — Robin convection to T_ambient
    """
    A    = np.zeros((N, N))
    B    = np.zeros((N, N))
    b_bc = np.zeros(N)

    # left boundary: just pin it to T_hot
    A[0, 0] = 1.0
    b_bc[0] = T_hot

    for i in range(1, N - 1):
        # pick the right dx; average it at layer interfaces
        if   i < end_steel1:   dx = dx_steel1
        elif i == end_steel1:  dx = 0.5 * (dx_steel1 + dx_wool)
        elif i < end_wool:     dx = dx_wool
        elif i == end_wool:    dx = 0.5 * (dx_wool + dx_steel2)
        else:                  dx = dx_steel2

        U = alpha[i] * dt / dx**2

        # standard CN stencil: half implicit, half explicit
        A[i, i-1] = -U/2;  A[i, i] = 1+U;  A[i, i+1] = -U/2
        B[i, i-1] =  U/2;  B[i, i] = 1-U;  B[i, i+1] =  U/2

    # right boundary: flux balance gives the Robin condition
    A[N-1, N-2] = -k_steel / dx_steel2
    A[N-1, N-1] =  k_steel / dx_steel2 + h
    b_bc[N-1]   =  h * T_ambient

    return A, B, b_bc


def run_simulation(N_steel1, N_wool, N_steel2, dt, t_end, store_snapshots=False, snapshot_times=None):
    """
    Runs the CN simulation over the three-layer domain.
    Returns final temperatures at the two interfaces and the cold face,
    plus the coordinate array and any saved snapshots.
    """
    dx_steel1 = L_steel1 / N_steel1
    dx_wool   = L_wool   / N_wool
    dx_steel2 = L_steel2 / N_steel2

    N          = N_steel1 + N_wool + N_steel2 + 1
    end_steel1 = N_steel1
    end_wool   = N_steel1 + N_wool

    x     = build_x_array(N_steel1, N_wool, N_steel2)
    alpha = build_alpha_array(N, end_steel1, end_wool)
    A, B, b_bc = build_cn_matrices(N, end_steel1, end_wool, alpha, dx_steel1, dx_wool, dx_steel2, dt)

    # start everything at ambient except the hot face
    T    = np.full(N, T_ambient)
    T[0] = T_hot

    snapshots = {0: T.copy()} if store_snapshots else {}
    t = 0.0

    for _ in range(int(t_end / dt)):
        T  = np.linalg.solve(A, B @ T + b_bc)
        t += dt

        if store_snapshots and snapshot_times is not None:
            for ts in snapshot_times:
                if abs(t - ts) < dt / 2:
                    snapshots[ts] = T.copy()

    return {
        'T_iface1'  : T[end_steel1],
        'T_iface2'  : T[end_wool],
        'T_cold'    : T[-1],
        'x'         : x,
        'snapshots' : snapshots,
    }


def find_converged_n(N_values, results, tol=0.01):
    """
    Walk up the refinement levels and stop when all three monitored
    temperatures change by less than tol. Falls back to the finest
    level if the tolerance is never met.
    """
    for i in range(1, len(N_values)):
        n_prev, n_curr = N_values[i-1], N_values[i]
        d1 = abs(results[n_curr][0] - results[n_prev][0])
        d2 = abs(results[n_curr][1] - results[n_prev][1])
        dc = abs(results[n_curr][2] - results[n_prev][2])
        if max(d1, d2, dc) < tol:
            return n_curr
    return N_values[-1]


def print_convergence_table(N_values, results, converged_n):
    """Print the spatial convergence table with the converged level marked."""
    print(f"  {'N/layer':>8}  {'T_iface1 (°C)':>14}  {'T_iface2 (°C)':>14}  {'T_cold (°C)':>12}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*12}")
    for n in N_values:
        T_i1, T_i2, T_c = results[n]
        marker = '  ← converged' if n == converged_n else ''
        print(f"  {n:>8}  {T_i1:>14.6f}  {T_i2:>14.6f}  {T_c:>12.6f}{marker}")

    print()
    print(f"  {'change':>8}  {'ΔT_iface1':>14}  {'ΔT_iface2':>14}  {'ΔT_cold':>12}")
    print(f"  {'-'*8}  {'-'*14}  {'-'*14}  {'-'*12}")
    for i in range(1, len(N_values)):
        n_prev, n_curr = N_values[i-1], N_values[i]
        d1 = abs(results[n_curr][0] - results[n_prev][0])
        d2 = abs(results[n_curr][1] - results[n_prev][1])
        dc = abs(results[n_curr][2] - results[n_prev][2])
        print(f"  {n_prev:>3}→{n_curr:<4}  {d1:>14.6f}  {d2:>14.6f}  {dc:>12.6f}")
    print()


def plot_spatial_convergence(ax, N_values, results, converged_n):
    """Plot temperature at key locations against nodes per layer."""
    T_i1 = [results[n][0] for n in N_values]
    T_i2 = [results[n][1] for n in N_values]
    T_c  = [results[n][2] for n in N_values]

    ax.plot(N_values, T_i1, 'o-', label='T_iface1', linewidth=2)
    ax.plot(N_values, T_i2, 's-', label='T_iface2', linewidth=2)
    ax.plot(N_values, T_c,  '^-', label='T_cold',   linewidth=2)

    ax.axvline(x=converged_n, color='black', linestyle='--', linewidth=1.5, label=f'Converged (N={converged_n})')

    ax.set_xlabel('Nodes per layer', fontsize=12)
    ax.set_ylabel('Temperature (°C)', fontsize=12)
    ax.set_title('Spatial Convergence', fontsize=13)
    ax.set_ylim([T_ambient - 10, T_hot + 10])
    ax.legend(fontsize=10)
    ax.grid(True, linestyle=':', alpha=0.6)


def plot_temperature_profiles(ax, x, snapshots, snapshot_times, converged_n):
    """Plot temperature profiles at each snapshot time with layer shading."""
    x_mm    = x * 1000
    L_total = L_steel1 + L_wool + L_steel2

    # shade each layer so it's obvious where the materials change
    ax.axvspan(0, L_steel1 * 1000, alpha=0.15, color='steelblue', label='Steel 1')
    ax.axvspan(L_steel1 * 1000, (L_steel1 + L_wool) * 1000, alpha=0.15, color='orange', label='Mineral Wool')
    ax.axvspan((L_steel1 + L_wool) * 1000, L_total * 1000, alpha=0.15, color='steelblue', label='Steel 2')

    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(snapshot_times)))
    for ts, col in zip(snapshot_times, colors):
        if ts in snapshots:
            ax.plot(x_mm, snapshots[ts], color=col, linewidth=2, label=f't = {ts} s')

    # dashed lines at the two interfaces
    ax.axvline(x=L_steel1 * 1000, color='black', linestyle='--', linewidth=1, alpha=0.5, label='Interfaces')
    ax.axvline(x=(L_steel1 + L_wool) * 1000, color='black', linestyle='--', linewidth=1, alpha=0.5)

    ax.set_xlabel('Position (mm)', fontsize=12)
    ax.set_ylabel('Temperature (°C)', fontsize=12)
    ax.set_title(f'Temperature Profiles  (N={converged_n} per layer)', fontsize=13)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xlim([0, L_total * 1000])
    ax.set_ylim([T_ambient - 10, T_hot + 10])
    ax.grid(True, linestyle=':', alpha=0.6)


# run spatial convergence study across a range of refinement levels
N_values = [5, 10, 20, 40, 80]
DT_FIXED = 0.00001   # 0.01 ms — small enough for these thin layers
T_END    = 0.05      # 50 ms — long enough to see the transient develop
CONV_TOL = 0.5       # 0.5 °C tolerance; N=5 and N=10 should fail this

print("=" * 65)
print(f"spatial convergence  (dt={DT_FIXED}, t_end={T_END})")
print("=" * 65)

spatial_results = {}
for n in N_values:
    r = run_simulation(n, n, n, dt=DT_FIXED, t_end=T_END)
    spatial_results[n] = (r['T_iface1'], r['T_iface2'], r['T_cold'])

converged_n = find_converged_n(N_values, spatial_results, tol=CONV_TOL)
print_convergence_table(N_values, spatial_results, converged_n)
print(f"  converged at N = {converged_n} per layer  (tolerance = {CONV_TOL} °C)\n")

# re-run at the converged resolution and save profiles at each snapshot time
snapshot_times = [0, 0.01, 0.02, 0.03, 0.04, 0.05]

r_converged = run_simulation(converged_n, converged_n, converged_n, dt=DT_FIXED, t_end=T_END, store_snapshots=True, snapshot_times=snapshot_times)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

plot_spatial_convergence(ax1, N_values, spatial_results, converged_n)
plot_temperature_profiles(ax2, r_converged['x'], r_converged['snapshots'], snapshot_times, converged_n)

plt.suptitle('Crank-Nicolson: Spatial Convergence & Temperature Profiles', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('convergence_and_profiles.png', dpi=150)
plt.show()