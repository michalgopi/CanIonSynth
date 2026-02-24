# Phantom Types

The system simulates two main classes of objects commonly scanned in industrial CT: Cans and Ionic Chambers.

## 1. Cans

Cans are simulated as cylindrical containers holding fluids and objects.

### Steel Cans
*   **Geometry**: Cylindrical body with a **flat bottom**.
*   **Characteristics**:
    *   Thinner walls (0.019–0.025 cm).
    *   Higher material density.
    *   Bottom constructed from two ellipsoids creating a flat profile.

### Aluminum Cans
*   **Geometry**: Cylindrical body with a **rounded bottom** (domed).
*   **Characteristics**:
    *   Thicker walls (0.028–0.044 cm).
    *   Lower material density.
    *   **Torus Geometry**: The bottom edge is modeled using a torus to create a smooth, curved transition to the walls. The center of the bottom features a spherical dome ("curvature").

### Internal Features
*   **Fluids**:
    *   **Single or Dual Phase**: Can contain one or two distinct fluid layers with different densities.
    *   **Meniscus**: Simulation of surface tension effects (curved fluid surface) at the walls.
    *   **Tilting**: Fluid surface can be tilted (simulating a non-level scan) using x/y cut angles.
*   **Objects**:
    *   **Balls**: Small spheres placed at the bottom (single or pairs).
    *   **Pipe**: A cylindrical pipe structure inserted into the fluid.
    *   **Dispenser**: A mechanism at the top of the can.

## 2. Ionic Chambers

Ionic chambers are precision instruments used for radiation detection. The simulation models their complex internal multi-material structure.

### Common Structure
All chambers share a concentric layer design:
1.  **Central Copper Electrode**: Innermost conductive element.
2.  **Inner Insulator**: Polyethylene layer.
3.  **Inner Aluminum Layer**: Shielding.
4.  **Outer Insulator**: Polyethylene layer.
5.  **Outer Aluminum Layer**: Secondary shielding.
6.  **Graphite Housing**: Outermost structural layer.
7.  **Air Cavity**: The sensitive volume between electrodes.

### Shapes
1.  **Ball-shaped**: Spherical top. Used for isotropic response.
2.  **Rounded-top Cylinder**: Cylindrical base with a half-ellipsoid cap.
3.  **Flat-top Cylinder** (Square-top): Cylindrical base with a flat top.
4.  **Lollipop-shaped**: Thin stem connected to a flat, wider cylindrical head.

### Standardized Sizes (`new_flat_sizes`)
For systematic testing, standardized volumes can be generated:
*   **Version 1**: ~60,000 mm³ (Radius 20mm, Height 48mm).
*   **Version 2**: ~30,000 mm³ (Radius 15mm, Height 42.9mm).
*   **Version 3**: ~10,000 mm³ (Radius 10mm, Height 32.5mm).
