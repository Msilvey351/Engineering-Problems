# Engineering-Problems
Engineering Problems I have faced and the code used to solve them


**Manifold_Torque**

Calculates worst-case torsional shear stress at the first thread root of a 1/2" NPT threaded hollow tube. Thread geometry is taken from ASME B1.20.1. 

Valid for 1/2" NPT threads only.

## What it does

Given a tube geometry and applied torque, it computes:

- Thread root diameter and remaining wall thickness accounting for NPT taper
- Polar moment of inertia at the reduced (thread root) cross-section
- Peak shear stress with stress concentration applied: tau = Kt * T * r_root / J_root
- A plot of shear stress vs torque with a shear strength limit


## Usage
Just adjust parameters and run

## Inputs
Adjust parameters based on manifold

Kt = 3.02 is calculated based on a radial notch with angle 60 degrees and radius 0.01mm. This has been found to be accurate from testing


## Thread Geometry
The 1/2" NPT constants are hardcoded from ASME B1.20.1. See standard for specific thread dimensions. 
The taper means the effective thread depth increases along the engagement length. The worst case is always at the small end, so that is what gets reported. Change if the thread snaps at an earlier point

## Limitations
- Hardcoded for 1/2" NPT only. Other thread sizes would need their own constants.
- Shear strength estimate of 0.6 * UTS is approximate




**Induction_Coil**

A Python simulation of a Half-Bridge Series Resonant (HBSR) induction cooker. Models the power electronics circuit and the resulting magnetic field and thermal behaviour in a pan.

Based on: Zungor & Bodur, "Design Methodology of Series Resonant Half Bridge Inverter for Induction Cooker", IEEE Access, 2023.

## Features

- Circuit: operating point (f_r, Q, ZVS status, I_peak, P_out, skin depth) 
- Field: magnetic field B at any point via numerical Biot-Savart integration
- Thermal: per-element eddy-current loss and temperature rise using the skin-depth approximation

## Usage
On startup it prints the default operating point at f_r, then presents a menu:

```
[1]  Operating point         (f_r, f_sw, Q, ZVS, I_peak, P_out, skin depth)
[2]  P(f_sw) curves          (power vs frequency for all capacitor values + 3D surface)
[3]  Harmonic breakdown       (bar chart of power per harmonic order)
[4]  B-field 3D vector plot   (Plotly cone plot over a 3D grid)
[5]  |B|^2 energy density map (2D heatmap at a chosen height z)
[6]  delta-T heatmap          (temperature rise map after t seconds)
```

Options 4 through 6 can take several minutes due to repeated Biot-Savart integration over the coil segments.

## Key Parameters
Edit the constants at the top of the file to match your hardware.

## Structure
The circuit solver computes the resonant frequency, Q factor, and peak coil current. That current is then passed into the field solver, which runs the Biot-Savart integration over the coil geometry. The thermal functions sit on top of the field solver and accumulate eddy-current losses over time.

I_peak is passed to the field solver. 

## Limitations
- Single-turn coil approximation only. Multi-turn coils would need superposition over all turns.
- The thermal model assumes the heat distributes uniformly through the full pan thickness. No convection or radiation is included.
- Does not account for cookwares' effect on the magnetic field




**Damped Spring**

A Python simulation of a mass-spring-damper system
using SciPy's ODE solver. Computes and plots the displacement and velocity of a mass over time given a set of physical parameters and initial conditions.

---

## Physics Background

The system follows the second-order equation of motion:

    m·ẍ + c·ẋ + k·x = 0

Where:
- `x` — displacement (m)
- `ẋ` — velocity (m/s)
- `ẍ` — acceleration (m/s²)
- `m` — mass (kg)
- `c` — damping coefficient (N·s/m)
- `k` — spring constant (N/m)

This is solved as a system of two first-order ODEs by treating displacement
and velocity as separate state variables.

### Damping Behaviour

The system's behaviour is determined by the **damping ratio** `ζ = c / (2√(mk))`:

| Damping Ratio | Behaviour |
|---|---|
| `ζ < 1` | Underdamped — oscillatory decay |
| `ζ = 1` | Critically damped — fastest return to equilibrium |
| `ζ > 1` | Overdamped — slow return, no oscillation |


## Configuration

All system parameters and initial conditions are defined at the top of the file
and can be freely adjusted:

    # System Parameters
    m = 1.0    # Mass (kg)
    k = 2.0    # Spring constant (N/m)
    c = 0.5    # Damping coefficient (N·s/m)

    # Initial Conditions
    x0 = 1.0   # Initial displacement (m)
    v0 = 0.0   # Initil velocity (m/s)

---

## Output

Running the scrit prints the system's physical properties to the console and produces a plot of velocity and displacement

