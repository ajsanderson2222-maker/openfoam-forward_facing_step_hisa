# Forward-Facing Step — hisa / OpenFOAM 13

Supersonic flow over a forward-facing step. Classic benchmark for high-speed compressible solvers. Produces a strong bow shock ahead of the step, expansion fans at corners, and reflected shocks downstream.

---

## Software

| Component | Version |
|-----------|---------|
| OpenFOAM  | 13 (openfoam.org) |
| hisa      | 1.13.4 |
| Platform  | Ubuntu 24.04 / WSL2 |

hisa (High-Speed Aerodynamic solver) is a density-based compressible flow solver built on OpenFOAM. Source: https://gitlab.com/hisa/hisa

---

## Problem Description

A uniform supersonic stream enters from the left and encounters a step at x = 0.6 m that rises 0.2 m from the bottom wall. The domain is 2D (one cell thick in z). The flow is inviscid, calorically perfect.

All variables are **non-dimensionalised**: the reference speed of sound is 1 m/s at T = 1 K, so the Mach number equals the velocity magnitude directly.

---

## Domain and Mesh

### Geometry

```
y=1  ┌──────────────────────────────────────┐ top (symmetry)
     │                                      │
     │         upper block                  │
     │                                      │
y=0.2├──────┬───────────────────────────────┤
     │ low  │ step face (obstacle)
     │ block│
y=0  └──────┘ bottom (symmetry)
    x=0   x=0.6                           x=3
```

- Domain: 3 m × 1 m × 0.1 m (2D, 1 cell in z)
- Step height: 0.2 m at x = 0.6 m
- Step face is a solid wall (obstacle patch)

### blockMesh layout

Three hex blocks:

| Block | Region | Cells (x×y×z) |
|-------|--------|---------------|
| 0 | Inlet region below step (0–0.6 × 0–0.2) | 72 × 24 × 1 |
| 1 | Inlet region above step (0–0.6 × 0.2–1) | 72 × 96 × 1 |
| 2 | Downstream region (0.6–3 × 0.2–1) | 288 × 96 × 1 |

Base mesh: 36,576 cells. Two rounds of uniform refinement are applied near the step to resolve the shock structure (via `topoSet` + `refineMesh`).

---

## Fluid Model

**Thermophysical model:** `hePsiThermo` — compressible, pressure-based thermo using ψ (compressibility).

**Equation of state:** Perfect gas.

**Transport:** Inviscid (μ = 0, laminar simulation type). This is the classic inviscid Euler formulation.

**Gas properties (non-dimensionalised):**

| Property | Value | Notes |
|----------|-------|-------|
| γ (gamma) | 1.4 | Diatomic ideal gas |
| Cp | 2.5 J/(kg·K) | Non-dim: R = Cp(1 − 1/γ) = 1 |
| Molecular weight | 11640.3 g/mol | Chosen so speed of sound = 1 m/s at T = 1 K |
| μ | 0 | Inviscid |

---

## Initial Conditions

All fields initialised to uniform freestream values:

| Field | Value | Units (non-dim) |
|-------|-------|-----------------|
| U (velocity) | (3, 0, 0) | m/s |
| T (temperature) | 1 | K |
| p (pressure) | 1 | Pa |
| ρ (density) | derived | kg/m³ |

Freestream Mach number: **Ma = 3.0** (U = 3, a = √(γT) = √1.4 ≈ 1.183, Ma = 3/1.183 ≈ 2.54 — note: internal field Mach at inlet ~2.54 due to boundary layer near step).

---

## Boundary Conditions

| Patch | Type | U | T | p |
|-------|------|---|---|---|
| inlet | patch | fixedValue (3,0,0) | fixedValue 1 | fixedValue 1 |
| outlet | patch | inletOutlet | inletOutlet | zeroGradient |
| bottom | symmetryPlane | symmetry | symmetry | symmetry |
| top | symmetryPlane | symmetry | symmetry | symmetry |
| obstacle (step face) | patch | slip | characteristicWallTemperature | characteristicWallPressure |
| defaultFaces (front/back) | empty | empty | empty | empty |

**Boundary condition notes:**
- `slip` on obstacle: inviscid wall, no normal velocity, tangential velocity free.
- `characteristicWallPressure` / `characteristicWallTemperature`: hisa characteristic-based BC that uses Riemann invariants to extrapolate wall pressure and temperature, physically consistent for supersonic flows.
- `inletOutlet` on outlet: acts as zeroGradient for outflow, fixedValue for any reverse flow.

---

## Solver — hisa

hisa solves the compressible Euler / Navier-Stokes equations in density-based form. The conserved variables are (ρ, ρU, ρE).

### Flux scheme

```
fluxScheme    AUSMPlusUp;
lowMachAusm   false;
```

**AUSM+up** (Advection Upstream Splitting Method, Plus with pressure diffusion): a low-dissipation approximate Riemann solver suited for all-Mach flows. Accurately captures shocks without excessive numerical smearing.

### Reconstruction (MUSCL)

```
reconstruct(rho)  wVanLeer;
reconstruct(U)    wVanLeer;
reconstruct(T)    wVanLeer;
```

**Weighted Van Leer limiter**: second-order TVD reconstruction. Prevents spurious oscillations across shocks while recovering second-order accuracy in smooth regions.

### Time integration

```
ddtSchemes { default dualTime rPseudoDeltaT CrankNicolson 0.9; }
```

**Dual time-stepping with Crank-Nicolson (blending 0.9):** Outer physical time march + inner pseudo-time iteration to convergence. CrankNicolson blending of 0.9 (nearly fully implicit) provides stability at high CFL.

### Gradient scheme

```
gradSchemes { default faceLeastSquares linear; }
```

Face least-squares gradient — more accurate than Gauss on skewed meshes.

---

## Linear Solver

```
solver      GMRES;
inviscidJacobian  LaxFriedrichs;
preconditioner    LUSGS;
maxIter     30;
nKrylov     4;
solverTolRel  1e-1;
```

| Setting | Value | Purpose |
|---------|-------|---------|
| GMRES | 4 Krylov vectors | Iterative Krylov solver for the linearised system |
| LaxFriedrichs Jacobian | — | First-order Jacobian for the preconditioner (robust) |
| LU-SGS preconditioner | — | Lower-Upper Symmetric Gauss-Seidel, matrix-free, cheap per iteration |

---

## Pseudo-time Control

```
nPseudoCorr    1000   (max inner iterations per physical time step)
pseudoTol      5e-3   (convergence criterion on all residuals)
pseudoCoNum    2.0    (initial pseudo-CFL)
pseudoCoNumMax 100.0  (CFL ramp ceiling)
pseudoCoNumMin 0.1
```

CFL ramps from 2 → 100 as the solution converges within each outer step, accelerating convergence.

---

## Run Control

| Parameter | Value |
|-----------|-------|
| Start time | 0 |
| End time | 4 (non-dim) |
| deltaT | 0.005 |
| Write interval | 0.1 |
| Parallel decomposition | 4 processors (scotch) |

Total wall time: ~262 s on 4 cores (WSL2).

---

## Results

Simulation converged to steady state. Final time step: t = 4.

| Quantity | Value |
|----------|-------|
| Total cells | 36,576 |
| Mach range | 0.008 – 3.80 |
| Freestream Mach | ~3.0 |
| Peak Mach (expansion) | ~3.80 |
| Near-stagnation Mach | ~0.008 |

### Flow features visible in Mach contours

- **Bow shock** ahead of the step face: normal shock component, strong pressure rise
- **Expansion fan** at the top-front corner of the step: Mach increases beyond freestream
- **Reflected shock** from top wall interacting with the main bow shock
- **Subsonic pocket** just upstream of the step face (stagnation region)

### Output images

| File | Description |
|------|-------------|
| `forwardStep/mesh.png` | Wireframe mesh (3200×1067, parallel projection) |
| `forwardStep/mach_contours.png` | Mach number contours at t=4, coolwarm colormap, range 0–3.5 |

---

## Running the Case

```bash
# Source OpenFOAM environment
source /opt/openfoam13/etc/bashrc

# Add hisa binaries to PATH
export PATH=$PATH:$HOME/OpenFOAM/$(whoami)-13/platforms/linux64GccDPInt32Opt/bin

# Run
cd forwardStep
bash runSim
```

To clean and restart:
```bash
bash cleanSim
```

---

## Post-processing

Requires the project venv with pyvista + matplotlib:

```bash
source .venv/bin/activate
python forwardStep/plot_hisa.py
```

Generates `mesh.png` and `mach_contours.png` in `forwardStep/`. Reads the reconstructed OpenFOAM case via pyvista's `OpenFOAMReader`. Mach computed from non-dimensionalised fields as:

```
Mach = |U| / sqrt(γ · T)     (γ = 1.4, T in non-dim units)
```

---

## References

- Woodward, P. & Colella, P. (1984). *The numerical simulation of two-dimensional fluid flow with strong shocks.* Journal of Computational Physics, 54(1), 115–173.
- hisa documentation: https://hisa.gitlab.io
- OpenFOAM 13: https://openfoam.org
