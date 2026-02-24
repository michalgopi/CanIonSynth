# Volume Calculation Methods

Accurate volume estimation is crucial for validating the generated phantoms. The system employs two complementary approaches.

## 1. Numerical Calculation
This method estimates volume by counting voxels in the generated binary masks.

*   **Formula**:
    $$ V_{\text{numerical}} = \sum(\text{mask}) \times (spacing_x \times spacing_y \times spacing_z) $$
*   **Pros**: Handles arbitrary complex shapes, intersections, and Boolean operations perfectly.
*   **Cons**: Subject to discretization errors (partial volume effects) at boundaries. Accuracy depends on resolution.

## 2. Analytical Calculation
This method calculates volume using mathematical formulas for the geometric primitives used.

*   **Formulas**:
    *   **Cylinder**: $V = \pi r^2 h$
    *   **Ellipsoid**: $V = \frac{4}{3}\pi r_x r_y r_z$
    *   **Torus**: $V = 2\pi^2 R r_{xy} r_z$
    *   **Spherical Cap**: $V = \frac{\pi h^2(3r - h)}{3}$
*   **Pros**: Mathematically precise for ideal shapes. Fast computation.
*   **Cons**: Difficult to apply to complex intersections (e.g., a tilted fluid plane cutting through a meniscus and a pipe).

## Fluid Volume Calculation Algorithm
The system uses a sophisticated algorithm (`compute_accurate_fluid_volume_fixed`) to estimate the analytical volume of fluid in a can, accounting for:

1.  **Base Volume**: Cylinder volume up to the fluid height.
2.  **Tilt Correction**: Adjusts the effective height based on the angle of the fluid surface cut.
3.  **Meniscus**: Adds the volume of the curved meniscus rim.
4.  **Bottom Geometry**:
    *   **Flat**: Adds/subtracts ellipsoid contributions.
    *   **Rounded**: Accounts for the torus and spherical dome volumes.
5.  **Subtractions**: Subtracts volumes of internal objects (pipes, balls) that displace fluid.

## Validation
The system compares Numerical and Analytical volumes during generation.
*   **Metric**: Relative Difference = $|V_{analytical} - V_{numerical}| / V_{numerical}$
*   **Target**: Typically < 5%. Larger differences may indicate resolution limits or mask generation issues.
