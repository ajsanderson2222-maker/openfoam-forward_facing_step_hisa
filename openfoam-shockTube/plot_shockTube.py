"""
Plot pressure, density, velocity, and temperature profiles along the shock tube
at all written time steps. Produces individual PNGs per variable plus a
combined 4-panel figure.
"""
import pyvista as pv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

pv.OFF_SCREEN = True

case_dir = "/home/andrew/openfoam_jobs/openfoam-shockTube/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-shockTube"

reader = pv.OpenFOAMReader(os.path.join(case_dir, "case.foam"))
times  = sorted(reader.time_values)
print(f"Times available: {times}")

# Collect 1D profiles by sampling along x at y=z=0
# Tube runs along x; collapse to 1D by sorting cell centres
profiles = {}

for t in times:
    reader.set_active_time_value(t)
    mesh = reader.read()["internalMesh"]

    centres = mesh.cell_centers().points          # (N,3)
    x       = centres[:, 0]
    order   = np.argsort(x)
    x       = x[order]

    p   = mesh.cell_data["p"][order]
    T   = mesh.cell_data["T"][order]
    rho = mesh.cell_data["rho"][order] if "rho" in mesh.cell_data.keys() else (mesh.cell_data["p"][order] / (287.0 * mesh.cell_data["T"][order]))
    U   = mesh.cell_data["U"][order]
    ux  = U[:, 0]

    profiles[t] = dict(x=x, p=p, T=T, rho=rho, ux=ux)

# Color map over time
cmap   = plt.cm.viridis
colors = [cmap(i / max(len(times) - 1, 1)) for i in range(len(times))]

fields = [
    ("p",  "Pressure (Pa)",        "pressure.png"),
    ("rho","Density (kg/m³)",       "density.png"),
    ("ux", "Velocity x (m/s)",      "velocity.png"),
    ("T",  "Temperature (K)",       "temperature.png"),
]

# Individual plots
for key, ylabel, fname in fields:
    fig, ax = plt.subplots(figsize=(10, 4))
    for i, t in enumerate(times):
        d = profiles[t]
        ax.plot(d["x"], d[key], color=colors[i], label=f"t={t:.3f}s", linewidth=1.2)
    ax.set_xlabel("x (m)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"Shock Tube — {ylabel}")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, fname), dpi=150)
    plt.close(fig)
    print(f"Saved {fname}")

# Combined 4-panel
fig, axes = plt.subplots(2, 2, figsize=(14, 8))
axes = axes.flatten()
for ax, (key, ylabel, _) in zip(axes, fields):
    for i, t in enumerate(times):
        d = profiles[t]
        ax.plot(d["x"], d[key], color=colors[i], label=f"t={t:.3f}s", linewidth=1.2)
    ax.set_xlabel("x (m)")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel)
    ax.legend(fontsize=7, loc="best")
    ax.grid(True, alpha=0.3)

fig.suptitle("Shock Tube — hisa / OpenFOAM 13", fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "shock_tube_profiles.png"), dpi=150)
plt.close(fig)
print("Saved shock_tube_profiles.png")
