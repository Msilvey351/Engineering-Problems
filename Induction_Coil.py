import math
import numpy as np
import sympy as sp
import scipy.integrate
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go

# =============================================================
#  HBSR Induction Cooker — Circuit + Field Simulation
#
#  Circuit: half-bridge series resonant inverter, harmonic sum method
#  Field:   Biot-Savart B-field, eddy-current loss, thermal rise
#
#  I_peak and f_sw from the circuit solver feed into the field solver.
#
#  Ref: Zungor & Bodur, "Design Methodology of Series Resonant
#       Half Bridge Inverter for Induction Cooker", IEEE Access, 2023.
#
#  Units: SI (should be at least)
# =============================================================


# --- Supply --------------------------------------------------

V_MAINS_RMS = 230.0
V_0 = V_MAINS_RMS * math.sqrt(2) / 2   # half-bridge output amplitude [V]

# --- Resonant Tank -------------------------------------------
R_LOSS = 40e-3                          # coil +circuit loss [Ω]
L_COIL = 95e-6                          # coil inductance                     [H]
C_VALUES = [0.68e-6, 0.56e-6, 0.39e-6]   # resonant capacitor options           [F]
C_SELECT = 0                              # index into C_VALUES
MAX_N = 5                              # harmonic truncation order; decay ~ 1/(2n+1)^2

# --- Load ----------------------------------------------------
# R_LOAD = pan resistance referred to coil terminals.
R_LOAD = 2.0   # [Ω] typical for cast iron

# --- Coil Geometry -------------------------------------------
COIL_SHAPE = 'octa'   # 'octa' | 'circle'
COIL_DIAM = 0.022    # used for circle shape only; octa segments hardcoded [m]

# --- Pan Material --------------------------------------------
RHO_PAN = 1e-7    # electrical resistivity                    [Ω·m]
MU_R = 5000    # relative permeability                     [-]
DENSITY = 7850    # mass density                              [kg/m³]
CP = 450     # specific heat capacity                    [J/kg·K]
PAN_THICK = 0.003   # base wall thickness for thermal mass calc [m]

# --- Heatmap Grid --------------------------------------------
GRID_POINTS = 30
SPAN = 0.030   # spatial domain width, ~2.5x coil radius [m]
Z_LEVEL = 0.002   # default evaluation plane above coil     [m]
HEAT_TIME = 60      # default integration time for delta-T    [s]

# --- Derived Constants ---------------------------------------
V_RMS = V_0                           # square wave RMS = amplitude for this topology
C     = C_VALUES[C_SELECT]
AREA  = (SPAN / GRID_POINTS) ** 2    # grid cell area                            [m²]
MU_0  = 4 * np.pi * 1e-7             # permeability of free space                [H/m]

# Sympy symbols for integration of Biot-Savart
_S, _X, _Y, _Z = sp.symbols('s x y z')



#  CIRCUIT MODEL
# =============================================================

def circuit_power(R_load, f_sw, C=C, L=L_COIL, R_loss=R_LOSS,
                  V_rms=V_RMS, n_harmonics=MAX_N):
    """
    Output power delivered to R_load at switching frequency f_sw.

    Harmonic superposition over odd Fourier components of the half-bridge
    square wave (Eq. 1, Zungor & Bodur):

        P = sum_{n=0}^{N} R_load * |V_n / Z_n|^2

         Z_n = (R_load + R_loss) + j*(w_n*L - 1/w_n*C)

    Parameters
    ----------
    R_load      : float — Pan load resistance [Ω]
    f_sw        : float — Switching frequency [Hz]
    C           : float — Resonant capacitance [F]
    L           : float — Coil inductance [H]
    R_loss      : float — Series loss resistance [Ω]
    V_rms       : float — Supply RMS voltage [V]
    n_harmonics : int   — Number of odd harmonic terms to sum [-]

    Returns
    -------
    float — Output power [W]
    """
    total = 0.0
    for n in range(n_harmonics):
        k   = 2 * n + 1
        w_n = 2 * np.pi * k * f_sw
        V_n = math.sqrt(2) * V_rms / (np.pi * k)
        X   = w_n * L - 1.0 / (w_n * C)            #total reactance
        total += R_load * V_n**2 / ((R_load + R_loss)**2 + X**2)
    return total


def resonant_frequency(L=L_COIL, C_cap=C):
    """
    Natural resonant frequency of the LC tank.

        f_r = 1 / (2*pi*sqrt(L*C))

    Parameters
    L     : float — Coil inductance [H]
    C_cap : float — Resonant capacitance [F]

    Returns
    float — Resonant frequency [Hz]
    """
    return 1.0 / (2 * np.pi * math.sqrt(L * C_cap))


def peak_coil_current(f_sw, C=C, L=L_COIL, R_loss=R_LOSS, V_rms=V_RMS, R_load=R_LOAD):
    """
    Peak coil current at f_sw, fundamental component only.

    This value is passed directly to the Biot-Savart solver.

    Parameters
    ----------
    f_sw   : float — Switching frequency [Hz]
    C      : float — Resonant capacitance [F]
    L      : float — Coil inductance [H]
    R_loss : float — Series loss resistance [Ω]
    V_rms  : float — Supply RMS voltage [V]
    R_load : float — Pan load resistance [Ω]

    Returns
    float — Peak coil current [A]
    """
    w   = 2 * np.pi * f_sw
    V_1 = math.sqrt(2) * V_rms / np.pi   # fundamental amplitude [V]
    X   = w * L - 1.0 / (w * C)
    Z   = math.sqrt((R_load + R_loss)**2 + X**2)
    return V_1 / Z


def solve_circuit(f_sw=None):
    """
    Compute and print the HBSR operating point.

    ZVS condition: f_sw >= f_r (tank inductive, IGBT turn-on at zero voltage).
    Operating below f_r gives capacitive tank - inductive is needed to prevent non-zero switching and stress on system

    Parameters
    f_sw : float — Switching frequency [Hz]; defaults to f_r if None

Returns
    tuple — (f_sw [Hz], I_peak [A], P_out [W])
    """
    f_r  = resonant_frequency()
    f_op = f_sw if f_sw is not None else f_r
    I_pk = peak_coil_current(f_op)
    P    = circuit_power(R_LOAD, f_op)
    Q    = math.sqrt(L_COIL / C) / (R_LOAD + R_LOSS)
    skin = math.sqrt(RHO_PAN / (np.pi * f_op * MU_0 * MU_R)) * 1e6   # [µm]
    zvs  = 'OK (inductive)' if f_op >= f_r else 'FAIL (capacitive)'

    print('\n' + '=' * 52)
    print('  HBSR Operating Point')
    print('=' * 52)
    print(f'  f_r        = {f_r/1e3:.2f} kHz')
    print(f'  f_sw       = {f_op/1e3:.2f} kHz')
    print(f'  Q          = {Q:.2f}')
    print(f'  ZVS        = {zvs}')
    print(f'  I_peak     = {I_pk:.2f} A')
    print(f'  P_out      = {P:.0f} W')
    print(f'  delta_skin = {skin:.2f} µm')
    print(f'  C          = {C*1e6:.2f} µF')
    print(f'  R_load     = {R_LOAD} Ω')
    print('=' * 52)

    return f_op, I_pk, P


#  CIRCUIT PLOTS
# =============================================================

def plot_power_vs_freq():
    """
    P(f_sw) curves for all capacitor values, plus P(f_sw, R_load) surface.

    Left panel:  power vs frequency at fixed R_LOAD, one curve per capacitor.
    Right panel: 3D surface sweeping R_load — shows sensitivity to pan material.
    """
    freq = np.linspace(500, 50000, 500)
    f_r  = resonant_frequency()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'HBSR  —  L={L_COIL*1e6:.0f} µH, 'f'R_load={R_LOAD} Ω, R_loss={R_LOSS*1e3:.0f} mΩ', fontweight='bold')

    ax = axes[0]
    for cap in C_VALUES:
        f_r_c = resonant_frequency(C_cap=cap)
        pwr   = [circuit_power(R_LOAD, f, C=cap) for f in freq]
        ax.plot(freq / 1e3, pwr, label=f'C={cap*1e6:.2f} µF  f_r={f_r_c/1e3:.1f} kHz')
    ax.axvline(f_r / 1e3, color='crimson', linestyle='--', linewidth=1.2, label=f'f_r = {f_r/1e3:.1f} kHz (selected C)')
    ax.set_xlabel('f_sw [kHz]')
    ax.set_ylabel('P_out [W]')
    ax.set_title('P vs f_sw')
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.set_xlim(0, 50)
    ax.set_ylim(0)

    ax3d   = fig.add_subplot(1, 2, 2, projection='3d')
    freq_s = np.linspace(500, 50000, 80)
    R_s    = np.linspace(0.1, 5.0, 60)
    FF, RR = np.meshgrid(freq_s, R_s)
    PP     = np.vectorize(lambda f, r: circuit_power(r, f, C=C))(FF, RR)
    surf   = ax3d.plot_surface(FF / 1e3, RR, PP, cmap='inferno', linewidth=0, antialiased=True, alpha=0.9)
    ax3d.plot([f_r/1e3]*2, [R_s[0], R_s[-1]], [0, 0], color='cyan', linewidth=2, label=f'f_r = {f_r/1e3:.1f} kHz')
    ax3d.set_xlabel('f_sw [kHz]', labelpad=6, fontsize=9)
    ax3d.set_ylabel('R_load [Ω]', labelpad=6, fontsize=9)
    ax3d.set_zlabel('P_out [W]', labelpad=6, fontsize=9)
    ax3d.set_title(f'P surface  (C = {C*1e6:.2f} µF)', fontsize=10)
    ax3d.view_init(elev=28, azim=-50)
    ax3d.legend(fontsize=8)
    fig.colorbar(surf, ax=ax3d, shrink=0.4, aspect=10, pad=0.12, label='P [W]')

    plt.tight_layout()
    plt.show()


def plot_harmonic_breakdown(f_sw):
    """
    Bar chart of P_n vs harmonic order at f_sw.

    Confirms truncation at MAX_N=5 is valid — P_n negligible by n=3 near resonance.

    Parameters
    f_sw : float — Switching frequency [Hz]
    """
    harmonics = list(range(8))
    P_h = []
    for n in harmonics:
        k   = 2 * n + 1
        w_n = 2 * np.pi * k * f_sw
        V_n = math.sqrt(2) * V_RMS / (np.pi * k)
        X   = w_n * L_COIL - 1.0 / (w_n * C)
        Z   = math.sqrt((R_LOAD + R_LOSS)**2 + X**2)
        P_h.append(R_LOAD * (V_n / Z)**2 / 2)

    plt.figure(figsize=(9, 5))
    bars = plt.bar(
        [f'n={n}  ({2*n+1}f)' for n in harmonics], P_h, color=['darkorange'] + ['steelblue'] * 7)
    for bar, p in zip(bars, P_h):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,f'{p:.1f} W', ha='center', va='bottom', fontsize=8)
    plt.title(f'Harmonic power breakdown  f_sw = {f_sw/1e3:.1f} kHz')
    plt.xlabel('Harmonic order')
    plt.ylabel('P_n [W]')
    plt.grid(True, axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()


#  COIL GEOMETRY
# =============================================================

class Coil:
    def __init__(self, shape, diam=0):
        """
        Parametric coil geometry for Biot-Savart integration.

        Parameters
        ----------
        shape : str   — 'circle' or 'octa'
        diam  : float — Coil diameter, used for circle shape only [m]
        """
        self.shape = shape
        self.diam = diam

        if shape == 'circle':
            # l(s) = (r*cos(s), r*sin(s), 0),  s in [0, 2pi]
            r = diam / 2
            self.l = r * sp.Matrix([[-sp.sin(_S), sp.cos(_S), 0.0]])
            self.diffl = sp.diff(self.l, _S)

        elif shape == 'octa':
            # Eight linear segments approximating the octagonal winding.
            # Endpoints from original drawing in mm; converted to m here.
            MM = 1e-3
            self.l = sp.Matrix([
                [ 9*MM,                                14*MM * _S / (2*np.pi),         0.0],
                [-9*MM,                               -14*MM * _S / (2*np.pi),         0.0],
                [ 11*MM * _S / (2*np.pi),             -11*MM,                          0.0],
                [-11*MM * _S / (2*np.pi),              11*MM,                          0.0],
                [-7.25*MM + 3.5*MM * _S/(2*np.pi),  -9*MM - 4*MM * _S/(2*np.pi),     0.0],
                [ 7.25*MM + 3.5*MM * _S/(2*np.pi),  -9*MM + 4*MM * _S/(2*np.pi),     0.0],
                [-7.25*MM - 3.5*MM * _S/(2*np.pi),   9*MM - 4*MM * _S/(2*np.pi),     0.0],
                [ 7.25*MM - 3.5*MM * _S/(2*np.pi),   9*MM + 4*MM * _S/(2*np.pi),     0.0],
            ])
            self.diffl = sp.diff(self.l, _S)

        else:
            raise ValueError(f"Unknown coil shape '{shape}'. Options: 'circle', 'octa'.")


#  BIOT-SAVART FIELD SOLVER
# =============================================================

def _biot_savart_integrands(coil):
    """
    Symbolic Biot-Savart integrand for each wire segment.

    Returns (dl x r_vec) / |r_vec|^3 per segment. Prefactor
    mu_0 * I / 4*pi is applied in B_field() after numerical integration.

    Parameters
    coil : Coil — Coil geometry object

    Returns
    list of sympy.Matrix — One (3,1) integrand vector per segment
    """
    r = sp.Matrix([_X, _Y, _Z])
    integrands = []
    for i in range(coil.l.shape[0]):
        l_seg  = coil.l[i, :].T
        dl_seg = coil.diffl[i, :].T
        vec    = r - l_seg
        integrands.append(dl_seg.cross(vec) / vec.norm()**3)
    return integrands


def B_field(x, y, z, coil, I_peak):
    """
    Magnetic field at (x, y, z) via numerical Biot-Savart integration.

        B = (mu_0 * I / 4*pi) * integral[ (dl x r_vec) / |r_vec|^3 ] ds

    I_peak is supplied by peak_coil_current() — this is the link between
    the circuit and field solvers.

    Parameters
    ----------
    x       : float — x-coordinate [m]
    y       : float — y-coordinate [m]
    z       : float — z-coordinate [m]
    coil    : Coil  — Coil geometry object
    I_peak  : float — Peak coil current [A]

    Returns
    -------
    np.ndarray — B field vector [T], shape (3,)
    """
    integrands = _biot_savart_integrands(coil)
    total      = np.zeros(3)
    for integrand in integrands:
        fx    = sp.lambdify([_S, _X, _Y, _Z], integrand[0])
        fy    = sp.lambdify([_S, _X, _Y, _Z], integrand[1])
        fz    = sp.lambdify([_S, _X, _Y, _Z], integrand[2])
        field = I_peak * MU_0 / (4 * np.pi) * np.array([
            scipy.integrate.quad_vec(fx, 0, 2*np.pi, args=(x, y, z))[0],
            scipy.integrate.quad_vec(fy, 0, 2*np.pi, args=(x, y, z))[0],
            scipy.integrate.quad_vec(fz, 0, 2*np.pi, args=(x, y, z))[0],
        ])
        total += field
    return total


def _eddy_power(x, y, z, coil, I_peak, freq):
    """
    Time-averaged eddy-current dissipation per grid element.

    Skin-depth approximation — currents confined to surface layer delta.
    H evaluated on air side of pan surface (no mu_r; field not yet refracted).

        P = rho * |H_s|^2 * A / (2*delta)

    Parameters
    x      : float — x-coordinate [m]
    y      : float — y-coordinate [m]
    z      : float — z-coordinate [m]
    coil   : Coil  — Coil geometry object
    I_peak : float — Peak coil current [A]
    freq   : float — Switching frequency [Hz]

    Returns
    float — Dissipated power [W]
    """
    delta = math.sqrt(RHO_PAN / (np.pi * freq * MU_0 * MU_R))   # skin depth [m]
    B_vec = B_field(x, y, z, coil, I_peak)
    H_sq = np.sum((B_vec / MU_0)**2)                            # |H_s|^2    [A²/m²]
    return RHO_PAN * H_sq * AREA / (2 * delta)


def _delta_T(x, y, z, coil, I_peak, freq, t):
    """
    Temperature rise of a pan base column after time t.

    Heat deposited in skin layer; assumed to conduct uniformly through
    PAN_THICK on cooking timescales (>> thermal diffusion time).

        delta-T = Q / (rho_m * V * Cp)   where V = AREA * PAN_THICK

    Parameters
    ----------
    x      : float — x-coordinate [m]
    y      : float — y-coordinate [m]
    z      : float — z-coordinate [m]
    coil   : Coil  — Coil geometry object
    I_peak : float — Peak coil current [A]
    freq   : float — Switching frequency [Hz]
    t      : float — Elapsed time [s]

    Returns
    -------
    float — Temperature rise [K]
    """
    Q = _eddy_power(x, y, z, coil, I_peak, freq) * t
    return Q / (DENSITY * AREA * PAN_THICK * CP)


def _build_grid():
    return np.linspace(-SPAN / 2, SPAN / 2, GRID_POINTS)


def _draw_heatmap(values, xy_vals, title, cmap, unit):
    plt.figure(figsize=(9, 7))
    sns.heatmap(values, xticklabels=np.round(xy_vals * 1e3, 1), yticklabels=np.round(xy_vals * 1e3, 1), cmap=cmap, cbar_kws={'label': unit},)
    plt.title(title)
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.tight_layout()
    plt.show()


#  FIELD PLOTS
# =============================================================

def plot_3d_vector_field(coil, I_peak):
    """
    B-field cone plot over a 10x10x10 grid (Plotly).

    Note: O(N^3) Biot-Savart evaluations x 8 segments — expect several minutes.

    Parameters
    ----------
    coil   : Coil  — Coil geometry object
    I_peak : float — Peak coil current [A]
    """
    print('Computing 3D B field — this will take a few minutes...')
    grid       = np.linspace(-SPAN, SPAN, 10)
    xg, yg, zg = np.meshgrid(grid, grid, grid)

    Bv         = np.vectorize(
        lambda x, y, z: B_field(x, y, z, coil, I_peak),
        signature='(),(),()->(n)'
    )(xg, yg, zg)
    Bx, By, Bz = Bv[..., 0], Bv[..., 1], Bv[..., 2]

    cone = go.Cone(
        x=xg.ravel(), y=yg.ravel(), z=zg.ravel(),
        u=Bx.ravel(), v=By.ravel(), w=Bz.ravel(),
        colorscale='Inferno', colorbar=dict(title='|B| [T]'),
        sizemode='scaled', sizeref=0.5,
    )
    layout = go.Layout(
        title=f'B field  I_peak={I_peak:.1f} A',
        scene=dict(
            xaxis_title='x [m]', yaxis_title='y [m]', zaxis_title='z [m]',
            aspectratio=dict(x=1, y=1, z=1)
        )
    )
    fig = go.Figure(data=cone, layout=layout)

    phi = np.linspace(-np.pi, np.pi, 100)
    L   = coil.l
    for i in range(L.shape[0]):
        row = L[i, :]
        lx  = [float(row[0].subs(_S, v)) for v in phi]
        ly  = [float(row[1].subs(_S, v)) for v in phi]
        lz  = [float(row[2].subs(_S, v)) for v in phi]
        fig.add_scatter3d(x=lx, y=ly, z=lz, mode='lines',
                          line=dict(color='lime', width=8),
                          name='coil' if i == 0 else f'seg {i}',
                          showlegend=(i == 0))
    fig.show()


def plot_em_energy(coil, I_peak, z):
    """
    |B|^2 heatmap at height z — proportional to EM energy density [T^2].

    Parameters
    ----------
    coil   : Coil  — Coil geometry object
    I_peak : float — Peak coil current [A]
    z      : float — Evaluation height above coil [m]
    """
    print(f'Computing |B|^2 map  z={z*1e3:.1f} mm  I_peak={I_peak:.1f} A...')
    xy   = _build_grid()
    vals = np.zeros((GRID_POINTS, GRID_POINTS))
    for xi, x in enumerate(xy):
        for yi, y in enumerate(xy):
            Bv           = B_field(x, y, z, coil, I_peak)
            vals[yi, xi] = np.sum(Bv**2)
    _draw_heatmap(
        vals, xy,
        f'|B|^2 [T^2]  z={z*1e3:.1f} mm  I_peak={I_peak:.1f} A',
        'inferno', '|B|^2 [T^2]'
    )


def plot_temp(coil, I_peak, freq, z, t):
    """
    delta-T heatmap [K] at height z after t seconds.

    Parameters
    ----------
    coil   : Coil  — Coil geometry object
    I_peak : float — Peak coil current [A]
    freq   : float — Switching frequency [Hz]
    z      : float — Evaluation height above coil [m]
    t      : float — Elapsed time [s]
    """
    print(f'Computing delta-T map  z={z*1e3:.1f} mm  t={t} s  I_peak={I_peak:.1f} A...')
    xy   = _build_grid()
    vals = np.zeros((GRID_POINTS, GRID_POINTS))
    for xi, x in enumerate(xy):
        for yi, y in enumerate(xy):
            vals[yi, xi] = _delta_T(x, y, z, coil, I_peak, freq, t)
    _draw_heatmap(
        vals, xy,
        f'delta-T [K]  z={z*1e3:.1f} mm  t={t} s  '
        f'f={freq/1e3:.1f} kHz  I_peak={I_peak:.1f} A',
        'viridis', 'delta-T [K]'
    )


#  INTERFACE
# =============================================================

MENU = """
  Induction Cooker Simulation
  ---------------------------
  Circuit
    [1]  Operating point  (f_r, f_sw, Q, ZVS, I_peak, P_out, skin depth)
    [2]  P(f_sw) curves + P(f_sw, R_load) surface
    [3]  Harmonic power breakdown P_n

  Field  — I_peak and f_sw set by circuit solver
    [4]  B-field 3D vector plot
    [5]  |B|^2 energy density heatmap
    [6]  delta-T heatmap

    [0]  Quit
"""


def get_freq():
    f_r = resonant_frequency()
    raw = input(f'  f_sw [Hz]  (Enter = f_r = {f_r/1e3:.2f} kHz): ').strip()
    return float(raw.replace(',', '')) if raw else f_r


def get_z():
    raw = input(f'  z [mm]  (Enter = {Z_LEVEL*1e3:.1f} mm): ').strip()
    return float(raw) / 1e3 if raw else Z_LEVEL


def get_time():
    raw = input(f'  t [s]  (Enter = {HEAT_TIME} s): ').strip()
    return float(raw) if raw else HEAT_TIME


def main():
    coil = Coil(shape=COIL_SHAPE, diam=COIL_DIAM)

    print('\nDefault operating point (f_sw = f_r):')
    solve_circuit(resonant_frequency())

    while True:
        print(MENU)
        choice = input('Select: ').strip()

        if choice == '0':
            break
        elif choice == '1':
            solve_circuit(get_freq())
        elif choice == '2':
            plot_power_vs_freq()
        elif choice == '3':
            f_sw = get_freq()
            solve_circuit(f_sw)
            plot_harmonic_breakdown(f_sw)
        elif choice == '4':
            f_sw = get_freq()
            _, I_pk, _ = solve_circuit(f_sw)
            print('  Runtime: O(N^3) — expect several minutes.')
            plot_3d_vector_field(coil, I_pk)
        elif choice == '5':
            f_sw = get_freq()
            _, I_pk, _ = solve_circuit(f_sw)
            z = get_z()
            plot_em_energy(coil, I_pk, z)
        elif choice == '6':
            f_sw = get_freq()
            _, I_pk, _ = solve_circuit(f_sw)
            z = get_z()
            t = get_time()
            plot_temp(coil, I_pk, f_sw, z, t)
        else:
            print(f"  Unknown option '{choice}'.")


#  Note: field solver valid for single-turn coil approximation only.
#  Multi-turn coils require superposition over all turns.
# =============================================================
if __name__ == '__main__':
    main()
""" explicitly stating these here, but they are set already
    COIL_SHAPE = 'octa'    # Coil shape                              [-]
    COIL_DIAM = 0.022     # Coil diameter (circle shape only)        [m]
    R_LOAD = 2.0       # Pan load resistance                      [Ω]
    Z_LEVEL = 0.002     # Evaluation height above coil             [m]
    HEAT_TIME = 60        # Heating duration for delta-T map         [s]
"""
    #main()