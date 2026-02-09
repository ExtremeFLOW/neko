# Initial Condition Mismatch Investigation: Advecting Drop

## Problem Statement

When comparing the analytical initial condition formula (from `advecting_drop.f90`) with the stored initial field (`field0.f00000`), we observe significant differences with maximum absolute error = 1.0.

## Analytical Formula (from Fortran code)

```fortran
rad = sqrt(x² + y² + z²)
phi = 0.5 * (1 + tanh((rad - 0.15) / (2*eps)))
```

**Key Parameters:**
- Droplet center: **(0, 0, 0)** - at the origin
- Droplet radius: **0.15**
- Interface width: **2 × epsilon**

## Root Cause: Domain Geometry Mismatch

### Domain Configuration

**From `genmeshbox.log`:**
```
xmin, xmax: 0.0 → 2.0
ymin, ymax: 0.0 → 2.0
zmin, zmax: -0.1 → 0.0
```

**Actual GLL point ranges** (from loaded field):
```
x: [0.013, 1.987]
y: [0.013, 1.987]
z: [-0.087, -0.013]
```

**Domain center:** (1.0, 1.0, -0.05)

### The Issue
The analytical formula places the droplet at **(0, 0, 0)**, which is the **bottom-left corner** of the domain, not the center. This explains the massive discrepancy.

## Evidence from Diagnostics

### 1. Error Magnitude
- **Maximum absolute error**: 1.0 (complete field inversion at some locations)
- **L2 error**: 0.51
- **Mean absolute error**: 0.28

### 2. Spatial Error Distribution
Points near origin (0, 0, 0):
```
Location: (0.013, 0.013, 0.0) → phi_file=1.0, phi_anal=0.0, error=1.0
Location: (0.013, 0.013, -0.013) → phi_file=1.0, phi_anal=0.0, error=1.0
```

Points far from origin:
```
Location: (1.967, 1.987, -0.013) → phi_file=0.0, phi_anal=1.0, error=-1.0
Location: (1.987, 1.967, -0.013) → phi_file=0.0, phi_anal=1.0, error=-1.0
```

### 3. Error by Radius from Origin
| Radius Range | Mean Error | Points |
|--------------|------------|--------|
| r < 0.10 | 0.02 | 5,346 |
| 0.10 ≤ r < 0.15 | 0.47 | 4,050 |
| 0.15 ≤ r < 0.20 | 0.50 | 6,318 |
| 1.00 ≤ r < 1.50 | 0.10 | 148,914 |
| r ≥ 1.50 | 0.49 | 112,320 |

The interface region (r ≈ 0.15) shows mean error ≈ 0.47, confirming the droplet is *not* centered at the origin in the stored field.

## The Mystery: Formula vs. Reality

**Confirmed Facts:**
1. The case file specifies `"initial_condition": {"type": "user"}` → uses Fortran formula
2. The Fortran formula is `rad = sqrt(x² + y² + z²)` → mathematically centers at origin (0,0,0)
3. The domain is [0, 2] × [0, 2] × [-0.1, 0] → origin is at corner
4. **BUT**: The stored field AND simulation show droplet at domain CENTER (1.0, 1.0)

**The Problem:**
- Mathematical analysis predicts droplet at corner (0,0,0)
- Actual simulation has droplet at center (1.0, 1.0)
- Both use the same code with `rad = sqrt(x² + y² + z²)`

**Possible Explanations to Investigate:**

1. **Hidden coordinate transformation:**
   - Perhaps mesh coordinates are transformed before `initial_conditions()` is called?
   - Or there's a coordinate shift in the mesh generation/loading?

2. **L2 projection effects:**
   - Spectral element L2 projection from DOF points to GLL nodes might introduce shifts?
   - But this seems unlikely to shift by exactly 1.0 units

3. **Misunderstanding of coordinate system:**
   - Perhaps `s%dof%x/y/z` doesn't reference absolute coordinates but local coordinates?

4. **Formula is actually correct for this coordinate system:**
   - Maybe the mesh/coordinate system is set up differently than `genmeshbox.log` suggests?

## Recommendations

1. **Verify the actual initial condition generation:**
   - Check if there's coordinate transformation in the field initialization
   - Review the mesh generation and coordinate system setup
   - Confirm the intended droplet center location

2. **Update analytical formula** to match the stored field:
   ```python
   # Correct centering at domain center
   center_x, center_y, center_z = 1.0, 1.0, 0.0  # or -0.05?
   rad = sqrt((x - center_x)² + (y - center_y)² + (z - center_z)²)
   phi = 0.5 * (1 + tanh((rad - 0.15) / (2*eps)))
   ```

3. **For analysis scripts:** Load `phi_initial` directly from `field0.f00000` rather than recomputing from the formula.

---

**Generated:** 2026-02-07
**Script:** `investigate_initial_condition_mismatch.py`
**Case:** `gamma_0.01_epsilon_0.02`
