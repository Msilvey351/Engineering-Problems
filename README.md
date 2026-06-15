````markdown
# Engineering-Problems
Engineering problems I have faced and the code used to solve them

---

## Manifold_Torque

Calculates worst-case torsional shear stress at the first thread root of a 1/2" NPT threaded hollow tube. Thread geometry is taken from ASME B1.20.1.

Valid for 1/2" NPT threads only.

### What it does

Given a tube geometry and applied torque, it computes:

- Thread root diameter and remaining wall thickness accounting for NPT taper
- Polar moment of inertia at the reduced (thread root) cross-section
- Peak shear stress with stress concentration applied: tau = Kt * T * r_root / J_root
- A plot of shear stress vs torque with a shear strength limit

### Usage

Just adjust parameters and run.

### Inputs

Adjust parameters based on manifold.

Kt = 3.02 is calculated based on a radial notch with angle 60 degrees and radius 0.01mm. This has been found to be accurate from testing.

### Thread Geometry

The 1/2" NPT constants are hardcoded from ASME B1.20.1. See standard for specific thread dimensions. The taper means the effective thread depth increases along the engagement length. The worst case is always at the small end, so that is what gets reported. Change if the thread snaps at an earlier point.

### Limitations

- Hardcoded for 1/2" NPT only. Other thread sizes would need their own constants.
- Shear strength estimate of 0.6 * UTS is approximate.

---

## Induction_Coil

A Python simulation of a Half-Bridge Series Resonant (HBSR) induction cooker. Models the power electronics circuit and the resulting magnetic field and thermal behaviour in a pan.

Based on: Zungor & Bodur, "Design Methodology of Series Resonant Half Bridge Inverter for Induction Cooker", IEEE Access, 2023.

### Features

- Circuit: operating point (f_r, Q, ZVS status, I_peak, P_out, skin depth)
- Field: magnetic field B at any point via numerical Biot-Savart integration
- Thermal: per-element eddy-current loss and temperature rise using the skin-depth approximation

### Usage

On startup it prints the default operating point at f_r, then presents a menu:

```
[1]  Operating point          (f_r, f_sw, Q, ZVS, I_peak, P_out, skin depth)
[2]  P(f_sw) curves           (power vs frequency for all capacitor values + 3D surface)
[3]  Harmonic breakdown       (bar chart of power per harmonic order)
[4]  B-field 3D vector plot   (Plotly cone plot over a 3D grid)
[5]  |B|^2 energy density map (2D heatmap at a chosen height z)
[6]  delta-T heatmap          (temperature rise map after t seconds)
```

Options 4 through 6 can take several minutes due to repeated Biot-Savart integration over the coil segments.

### Key Parameters

Edit the constants at the top of the file to match your hardware.

### Structure

The circuit solver computes the resonant frequency, Q factor, and peak coil current. That current is then passed into the field solver, which runs the Biot-Savart integration over the coil geometry. The thermal functions sit on top of the field solver and accumulate eddy-current losses over time.

I_peak is passed to the field solver.

### Limitations

- Single-turn coil approximation only. Multi-turn coils would need superposition over all turns.
- The thermal model assumes the heat distributes uniformly through the full pan thickness. No convection or radiation is included.
- Does not account for cookwares' effect on the magnetic field.

---

## Damped_Spring

A Python simulation of a mass-spring-damper system using SciPy's ODE solver. Computes and plots the displacement and velocity of a mass over time given a set of physical parameters and initial conditions.

### Physics Background

The system follows the second-order equation of motion:

    m·ẍ + c·ẋ + k·x = 0

Where:
- `x` — displacement (m)
- `ẋ` — velocity (m/s)
- `ẍ` — acceleration (m/s²)
- `m` — mass (kg)
- `c` — damping coefficient (N·s/m)
- `k` — spring constant (N/m)

This is solved as a system of two first-order ODEs by treating displacement and velocity as separate state variables.

### Damping Behaviour

The system's behaviour is determined by the damping ratio `ζ = c / (2√(mk))`:

| Damping Ratio | Behaviour |
|---|---|
| ζ < 1 | Underdamped — oscillatory decay |
| ζ = 1 | Critically damped — fastest return to equilibrium |
| ζ > 1 | Overdamped — slow return, no oscillation |

### Configuration

All system parameters and initial conditions are defined at the top of the file and can be freely adjusted:

```
m = 1.0    # Mass (kg)
k = 2.0    # Spring constant (N/m)
c = 0.5    # Damping coefficient (N·s/m)

x0 = 1.0   # Initial displacement (m)
v0 = 0.0   # Initial velocity (m/s)
```

### Output

Running the script prints the system's physical properties to the console and produces a plot of velocity and displacement.

---

## Crank_Nicolson

Simulates transient heat conduction through a three-layer sandwich (steel / mineral wool / steel) using the finite difference method (Crank-Nicolson discretisation of the conduction equation). The script runs a spatial convergence study to confirm convergence then plots the temperature profile at several time steps.

### Problem Setup

A thin three-layer wall is exposed to a fixed hot temperature (Dirichlet) on one side and convective cooling (Robin) on the other:

| Layer | Material | Thickness |
|---|---|---|
| 1 | Steel | 1.2 mm |
| 2 | Mineral Wool | 0.5 mm |
| 3 | Steel | 1.2 mm |

**Boundary conditions:**
- Left face: Dirichlet — fixed at 500 °C
- Right face: Robin — convection to 20 °C ambient with h = 500 W/m²K

### Method

The 1D transient conduction equation is:

    ∂T/∂t = α ∂²T/∂x²

The domain is discretised into three separate uniform grids (one per layer) that are stitched together. At each interface the grid spacing is averaged and the thermal diffusivity takes the value of the softer material (wool).

**Spatial discretisation**

The second derivative is approximated with a standard central difference:

    ∂²T/∂x² ≈ (T[i-1] - 2T[i] + T[i+1]) / dx²

**Crank-Nicolson time stepping**

Rather than evaluating the right-hand side at the old time level (explicit) or the new time level (fully implicit), Crank-Nicolson averages the two:

    (T[i]^(n+1) - T[i]^n) / dt = (α / dx²) * (1/2)(δ²T^n + δ²T^(n+1))

where δ²T denotes the central difference stencil. This gives second-order accuracy in time and is unconditionally stable, which matters here because α_steel >> α_wool — an explicit scheme would need an extremely small timestep to stay stable in the steel layers.

**Rearranging into matrix form**

Collecting the n+1 terms on the left and the n terms on the right, and defining:

    U = α[i] * dt / dx²

the update equation at each interior node i becomes:

    -U/2 * T[i-1]^(n+1)  +  (1+U) * T[i]^(n+1)  -  U/2 * T[i+1]^(n+1)
        = U/2 * T[i-1]^n  +  (1-U) * T[i]^n  +  U/2 * T[i+1]^n

Writing this for every node simultaneously gives:

    A @ T_new = B @ T_old + b_bc

where A carries the implicit (n+1) coefficients, B carries the explicit (n) coefficients, and b_bc is a vector that holds the boundary contributions.

**Boundary conditions**

The left face is Dirichlet. Row 0 of A is set to the identity and b_bc[0] = T_hot, which simply pins the first node every timestep.

The right face is a Robin (convective) condition. An energy balance on the half-cell at the last node gives:

    k * (T[N-2] - T[N-1]) / dx = h * (T[N-1] - T_ambient)

Rearranging:

    -k/dx * T[N-2]  +  (k/dx + h) * T[N-1]  =  h * T_ambient

This fills the last row of A directly, with h * T_ambient going into b_bc[N-1]. B is left as zero for that row since the boundary condition replaces the stencil entirely. The system is then solved with a standard direct solver at each timestep.

### Spatial Convergence Study

The script tests five mesh refinement levels — 5, 10, 20, 40, and 80 nodes per layer — and monitors the temperature at three locations:

- T_iface1 — steel/wool interface
- T_iface2 — wool/steel interface
- T_cold — cold face (right boundary)

The converged mesh is the first level where all three temperatures change by less than 0.5 °C relative to the previous refinement. Results are printed as a table with the converged level marked.

### Output

Running the script produces one figure saved as `convergence_and_profiles.png` with two panels:

- Left — temperatures at the three monitored locations plotted against nodes per layer, with a vertical line at the converged mesh
- Right — full temperature profiles across the wall at t = 0, 0.01, 0.02, 0.03, 0.04, and 0.05 seconds, with layer shading

### Parameters

These are all set near the top of the script and are easy to change:

| Variable | Value | Description |
|---|---|---|
| DT_FIXED | 1e-5 s | Timestep |
| T_END | 0.05 s | Simulation end time |
| CONV_TOL | 0.5 °C | Convergence tolerance |
| T_hot | 500 °C | Hot face temperature |
| T_ambient | 20 °C | Ambient temperature |
| h | 500 W/m²K | Convection coefficient |

### Disclaimer

This solution was developed with assistance from a large language model (Claude Sonnet 4.6). The LLM was used specifically for generating the matplotlib plotting code and for general code cleanup and formatting. The numerical method, boundary conditions, material properties, and convergence logic were written and verified independently.
````