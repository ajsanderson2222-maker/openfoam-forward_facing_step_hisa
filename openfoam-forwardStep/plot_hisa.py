"""
Generate mesh PNG and Mach contour PNG from hisa forwardStep solution using pyvista.
No display required — runs fully offscreen.
"""
import pyvista as pv
import numpy as np
import os

pv.OFF_SCREEN = True
pv.global_theme.background = "white"

case_dir = "/home/andrew/openfoam_jobs/openfoam-forwardStep/simulation"
out_dir = "/home/andrew/openfoam_jobs/openfoam-forwardStep"

# Find latest time directory (reconstructed)
time_dirs = []
for d in os.listdir(case_dir):
    try:
        time_dirs.append(float(d))
    except ValueError:
        pass
latest_time = max(time_dirs)
print(f"Latest time: {latest_time}")

# Load mesh via OpenFOAM reader
foam_file = os.path.join(case_dir, "case.foam")
reader = pv.OpenFOAMReader(foam_file)
reader.set_active_time_value(latest_time)
mesh = reader.read()

# Internal mesh block
internal = mesh["internalMesh"]
print(f"Mesh: {internal.n_cells} cells, bounds: {internal.bounds}")

# ── Mesh plot ─────────────────────────────────────────────────────────────────
b = internal.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
cx = (b[0] + b[1]) / 2
cy = (b[2] + b[3]) / 2
# window matches domain aspect ratio: x_span / y_span
x_span = b[1] - b[0]
y_span = b[3] - b[2]
scale = 3200
win_w = scale
win_h = max(int(scale * y_span / x_span), 1)

pl = pv.Plotter(off_screen=True, window_size=[win_w, win_h])
pl.add_mesh(internal, style="wireframe", color="black", line_width=0.4)
pl.background_color = "white"
pl.view_xy()
# fit camera exactly to domain extents with small padding
pl.camera.parallel_projection = True
pl.camera.position = (cx, cy, 1.0)
pl.camera.focal_point = (cx, cy, 0.0)
pl.camera.up = (0, 1, 0)
pl.camera.parallel_scale = y_span / 2 * 1.02  # 2% padding
pl.screenshot(os.path.join(out_dir, "mesh.png"))
print("Saved mesh.png")

# ── Mach contours ─────────────────────────────────────────────────────────────
# Convert cell data to point data for smoother contours
internal_pt = internal.cell_data_to_point_data()

U = internal_pt["U"]
T = internal_pt["T"]
gamma = 1.4
# hisa uses non-dimensionalised fields: a = sqrt(gamma * T), no R factor
speed_of_sound = np.sqrt(gamma * T)
mach = np.linalg.norm(U, axis=1) / speed_of_sound
internal_pt["Mach"] = mach

print(f"Mach range: {mach.min():.3f} – {mach.max():.3f}")

pl2 = pv.Plotter(off_screen=True, window_size=[1600, 800])
pl2.add_mesh(
    internal_pt,
    scalars="Mach",
    cmap="coolwarm",
    clim=[0.0, min(3.5, mach.max())],
    scalar_bar_args={"title": "Mach", "vertical": True},
)
pl2.view_xy()
pl2.background_color = "white"
pl2.screenshot(os.path.join(out_dir, "mach_contours.png"))
print("Saved mach_contours.png")
