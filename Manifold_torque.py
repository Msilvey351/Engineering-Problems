import math
import matplotlib.pyplot as plt
import numpy as np

# =============================================================
#  Torsional Shear Stress — 1/2" NPT Threaded Hollow Tube
#  Valid for 1/2" NPT thread only
#
#  τ_actual = Kt * T * r_root / J_root
# =============================================================

# --- 1/2" NPT Thread Constants (ASME B1.20.1) ---------------
NPT_HALF_INCH = {
    "name": '1/2" NPT',
    "TPI": 14,               # Threads per inch
    "h_mm": 0.05714 * 25.4,    # Thread depth = 1.451 mm (radial)
    "taper_per_mm": (0.0625 / 25.4),  # Diametral taper per mm of length
}


def calculate_npt_thread_depth(thread_length_mm: float) -> dict:
    """
    Calculate the maximum effective thread depth for a 1/2" NPT thread
    over a given engagement length. The worst case (deepest) is always
    at the small end of the taper.

    Note: Only valid for 1/2" NPT threads.

    Parameters
    ----------
    thread_length_mm : float — Thread engagement length [mm]

    Returns
    -------
    dict with thread geometry values
    """
    h     = NPT_HALF_INCH["h_mm"]
    taper = NPT_HALF_INCH["taper_per_mm"]

    taper_contribution_radial = (taper * thread_length_mm) / 2
    max_radial_depth          = h + taper_contribution_radial

    return {
        "Thread Standard": NPT_HALF_INCH["name"],
        "TPI": NPT_HALF_INCH["TPI"],
        "Thread Engagement Length (mm)": thread_length_mm,
        "Constant Thread Depth h (mm)": h,
        "Taper Contribution over length (mm)": taper_contribution_radial,
        "Max Radial Depth (mm)": max_radial_depth,
    }


def polar_moment_of_inertia(d_outer_mm: float, d_inner_mm: float) -> float:
    """
    Polar moment of inertia (J) for a hollow circular section.

    Parameters
    ----------
    d_outer_mm : float — Outer diameter [mm]
    d_inner_mm : float — Inner diameter [mm]

    Returns
    -------
    float — J [mm⁴]
    """
    if d_inner_mm >= d_outer_mm:
        raise ValueError("Inner diameter must be less than outer diameter.") # catch for invalid geometry
    return (math.pi / 32) * (d_outer_mm**4 - d_inner_mm**4)


def calculate_thread_shear_stress(
    T_Nm             : float,
    d_outer_mm       : float,
    d_inner_mm       : float,
    thread_length_mm : float,
    Kt               : float,
) -> dict:
    """
    Calculate shear stress at the first thread root (worst case — small end)
    for a 1/2" NPT threaded hollow tube.

    Also reports nominal wall thickness and minimum remaining wall thickness
    at the thread root (the thinnest, structurally critical cross-section).

    Note: Only valid for 1/2" NPT threads.

    Parameters
    ----------
    T_Nm: float — Applied torque [N·m]
    d_outer_mm: float — Outer diameter [mm]
    d_inner_mm: float — Inner diameter [mm]
    thread_length_mm: float — Thread engagement length [mm]
    Kt: float — Stress concentration factor [-]

    Returns
    -------
    dict with all computed values
    """
    thread_info = calculate_npt_thread_depth(thread_length_mm)
    max_radial_depth = thread_info["Max Radial Depth (mm)"]

    d_root_mm = d_outer_mm - (2 * max_radial_depth)
    r_root_mm = d_root_mm / 2

    # ── Wall thickness values ───────────────────────────────
    nominal_wall_mm = (d_outer_mm - d_inner_mm) / 2
    min_wall_mm = (d_root_mm  - d_inner_mm) / 2

    if r_root_mm <= d_inner_mm / 2:
        raise ValueError(
            f"Thread depth ({max_radial_depth:.2f} mm) exceeds wall thickness "
            f"({(d_outer_mm - d_inner_mm) / 2:.2f} mm) — invalid geometry."
        )

    r_root_m = r_root_mm / 1000
    J_root = polar_moment_of_inertia(d_root_mm, d_inner_mm) / 1e12  # mm⁴ → m⁴

    tau_Pa = Kt * (T_Nm * r_root_m) / J_root
    tau_MPa = tau_Pa / 1e6

    return {
        **thread_info,
        "--- Tube Geometry ---"                : "---",
        "Outer Diameter (mm)"                  : d_outer_mm,
        "Inner Diameter (mm)"                  : d_inner_mm,
        "Nominal Wall Thickness (mm)"          : nominal_wall_mm,
        "Thread Root Diameter (mm)"            : d_root_mm,
        "Thread Root Radius (mm)"              : r_root_mm,
        "Min Wall at Thread Root (mm)"         : min_wall_mm,
        "--- Load ---"                         : "---",
        "Applied Torque T (N·m)"               : T_Nm,
        "--- Stress ---"                       : "---",
        "Stress Concentration Kt"              : Kt,
        "J_root (m⁴)"                          : J_root,
        "τ_thread (Pa)"                        : tau_Pa,
        "τ_thread (MPa)"                       : tau_MPa,
    }


def print_results(results: dict) -> None:
    """Pretty-print the results."""
    print("\n" + "=" * 64)
    print('  1/2" NPT Threaded Tube — Torsional Shear Stress (Worst Case)')
    print("=" * 64)
    for label, value in results.items():
        if isinstance(value, float):
            print(f"  {label:<46} {value:.6g}")
        else:
            print(f"  {label:<46} {value}")
    print("=" * 64 + "\n")


# =============================================================
#  GRAPH — Shear Stress vs Torque
# =============================================================

def plot_thread_shear_stress(
    d_outer_mm: float,
    d_inner_mm: float,
    thread_length_mm: float,
    Kt: float,
    UTS_MPa: float,
    T_max_Nm: float,
    T_applied_Nm: float = None,
) -> None:
    """
    Plot shear stress at the 1/2" NPT first thread root vs applied torque,
    with a dashed line at the shear strength limit (0.6 × UTS) and an
    annotation showing nominal and minimum wall thickness.

    Note: Only valid for 1/2" NPT threads.

    Parameters
    ----------
    d_outer_mm: float — Outer diameter [mm]
    d_inner_mm: float — Inner diameter [mm]
    thread_length_mm: float — Thread engagement length [mm]
    Kt: float — Stress concentration factor [-]
    UTS_MPa: float — Ultimate Tensile Strength [MPa]
    T_max_Nm: float — Maximum torque on x-axis [N·m]
    T_applied_Nm: float — (Optional) highlight a specific applied torque [N·m]
    """

    thread_info = calculate_npt_thread_depth(thread_length_mm)
    max_radial_depth = thread_info["Max Radial Depth (mm)"]
    d_root_mm = d_outer_mm - (2 * max_radial_depth)
    r_root_m = (d_root_mm / 2) / 1000
    J_root = polar_moment_of_inertia(d_root_mm, d_inner_mm) / 1e12

    nominal_wall_mm = (d_outer_mm - d_inner_mm) / 2
    min_wall_mm = (d_root_mm  - d_inner_mm) / 2

    shear_strength_MPa = 0.6 * UTS_MPa

    T_values = np.linspace(0, T_max_Nm, 500)
    tau_MPa = Kt * (T_values * r_root_m) / J_root / 1e6

    T_failure = (shear_strength_MPa * 1e6 * J_root) / (Kt * r_root_m)

    # ── Plot ──────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(T_values, tau_MPa,
            color="darkorange", linewidth=2.5,
            label=rf"Thread Root Stress  $\tau = K_t \cdot T \cdot r_{{root}}/J_{{root}}$"
                  rf"  ($K_t={Kt}$)")

    ax.axhline(shear_strength_MPa, color="crimson", linewidth=2,
               linestyle="--",
               label=rf"Shear Strength $\approx 0.6 \times$ UTS = {shear_strength_MPa:.1f} MPa")

    ax.axvline(T_failure, color="crimson", linewidth=1.5,
               linestyle=":",
               label=f"Failure Torque = {T_failure:.1f} N·m")

    ax.fill_betweenx([0, max(tau_MPa)], T_failure, T_max_Nm,
                     color="crimson", alpha=0.15, label="Failure Zone")

    if T_applied_Nm is not None:
        tau_app = Kt * (T_applied_Nm * r_root_m) / J_root / 1e6
        ax.scatter(T_applied_Nm, tau_app, color="darkorange", zorder=5, s=70)
        ax.annotate(
            f"T = {T_applied_Nm} N·m\n"
            f"τ = {tau_app:.2f} MPa",
            xy=(T_applied_Nm, tau_app),
            xytext=(T_applied_Nm + T_max_Nm * 0.05, tau_app + shear_strength_MPa * 0.08),
            arrowprops=dict(arrowstyle="->", color="darkorange"),
            fontsize=9, color="darkorange",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="darkorange", alpha=0.8)
        )

    wall_text = (
        f"Nominal wall thickness : {nominal_wall_mm:.3f} mm\n"
        f"Min wall at thread root: {min_wall_mm:.3f} mm"
    )
    ax.text(
        0.98, 0.05, wall_text,
        transform=ax.transAxes,
        fontsize=8.5, verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.4", fc="white",
                  ec="steelblue", alpha=0.85),
        color="steelblue",
    )

    ax.set_xlabel("Applied Torque  T  [N·m]", fontsize=12)
    ax.set_ylabel("Shear Stress  τ  [MPa]",   fontsize=12)
    ax.set_title(
        '1/2" NPT — Worst-Case Shear Stress vs Applied Torque\n'
        rf"($D_o$={d_outer_mm:.1f} mm, $D_i$={d_inner_mm:.1f} mm, "
        rf"Thread Length={thread_length_mm:.1f} mm, $K_t$={Kt})",
        fontsize=11
    )

    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(0, T_max_Nm)
    ax.set_ylim(0)
    ax.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig("npt_half_inch_shear_stress.png", dpi=150)
    plt.show()


# =============================================================
#  note: Only valid for 1/2" NPT threads
# =============================================================
if __name__ == "__main__":

    T_applied_Nm     = 70.0          # Applied torque                   [N·m]
    d_outer_mm       = 22            # Tube outer diameter              [mm]
    d_inner_mm       = 18            # Tube inner diameter              [mm]
    thread_length_mm = 17            # Thread engagement length         [mm]
    Kt               = 3.02          # Stress concentration factor      [-]
    UTS_MPa          = 400.0         # UTS of material                  [MPa]
    T_max_Nm         = 200.0         # Max torque on x-axis             [N·m]


    results = calculate_thread_shear_stress(
        T_applied_Nm, d_outer_mm, d_inner_mm,
        thread_length_mm, Kt
    )
    print_results(results)

    plot_thread_shear_stress(
        d_outer_mm       = d_outer_mm,
        d_inner_mm       = d_inner_mm,
        thread_length_mm = thread_length_mm,
        Kt               = Kt,
        UTS_MPa          = UTS_MPa,
        T_max_Nm         = T_max_Nm,
        #T_applied_Nm     = T_applied_Nm,
    )