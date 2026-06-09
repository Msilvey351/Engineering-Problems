# Engineering-Problems
Engineering Problems I have faced and the code used to solve them


**Manifold_torque**
Uses calculates the torque needed to deform or break the thread on a cylindrical tube with a 1/2" NPT thread. 
It takes inputs of material ultimate tensile strength, inner and outer (nominal) tube diameters, length of thread engaged (ie, length of thread to the break point), and a stress concentration factor Kt (calculation for Kt done seperately). 
It plots the torque vs the stress on the break point, highlighting the region where the stress is greater than shear strength and material will deform. 


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



