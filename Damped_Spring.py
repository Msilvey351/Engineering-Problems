from scipy.integrate import solve_ivp
from numpy import linspace, sqrt
import matplotlib.pyplot as plt

# System Parameters
m = 1.0    # Mass (kg)
k = 2.0    # Spring constant (N/m)
c = 0.5    # Damping coefficient (N·s/m)

# Initial Conditions
x0 = 1.0   # Initial displacement (m)
v0 = 0.0   # Initial velocity (m/s)

# Natural frequency (rad/s): the frequency at which the system oscillates without damping
omega_n = sqrt(k / m)

# Damping ratio (dimensionless)
zeta = c / (2 * sqrt(m * k))

print(f"Natural frequency: {omega_n:.3f} rad/s")
print(f"Damping ratio: {zeta:.3f}")

if zeta < 1:
    print("System is underdamped (oscillatory decay)")
elif zeta == 1:
    print("System is critically damped")
else:
    print("System is overdamped (no oscillation)")


def derivatives(t, state, m, k, c):
    """
    Defines the system of first-order ODEs for a damped mass-spring system.

    The second-order equation of motion:
        m·ẍ + c·ẋ + k·x = 0
    is decomposed into two first-order equations by treating displacement
    and velocity as separate state variables.

    Parameters
    t : float
        Current time (s). Not used explicitly since the system is
        autonomous (time-independent), but required by the solver interface.
    state : array-like of shape (2,)
        Current state vector where:
            state[0] = x  : displacement (m)
            state[1] = ẋ  : velocity (m/s)
    m : float
        Mass (kg)
    k : float
        Spring constant (N/m)
    c : float
        Damping coefficient (N·s/m)

    Returns
    list of float
        Derivatives of the state vector:
            [velocity, acceleration] → [ẋ, ẍ]
    """

    velocity = state[1]

    # Compute acceleration from the equation of motion:
    # ẍ = -(k·x + c·ẋ) / m
    acceleration = -(k * state[0] + c * state[1]) / m

    return [velocity, acceleration]


t = linspace(0, 10, 1000)

# Solve the ODE System
soln = solve_ivp(derivatives, (0, 10), (x0, v0), t_eval=t, args=(m, k, c))

# Check Solver Success
if not soln.success:
    raise RuntimeError(f"ODE solver failed: {soln.message}")

# Plot Results
fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

# Displacement plot
axes[0].plot(t, soln.y[0], color="steelblue", label="Displacement")
axes[0].set_ylabel("Displacement (m)")
axes[0].legend()
axes[0].grid()

# Velocity plot
axes[1].plot(t, soln.y[1], color="tomato", label="Velocity")
axes[1].set_ylabel("Velocity (m/s)")
axes[1].set_xlabel("Time (s)")
axes[1].legend()
axes[1].grid()

fig.suptitle(
    f"Damped Mass-Spring System  |  ζ={zeta:.2f}, ωₙ={omega_n:.2f} rad/s",
    fontsize=13
)
plt.tight_layout()
plt.savefig("damped_spring.png", dpi=150, bbox_inches="tight")
print("Plot saved to damped_spring.png")