# Volume Calculation Methods

Accurate volume estimation is one of the main validation signals for this repository. The code compares voxel-count-based measurements against analytical geometric estimates wherever possible.

## Where The Logic Lives

- `in_docker_organized/volume_integration.jl`: can phantom analytical and numerical volume comparison
- `in_docker_organized/main_create_phantom_can.jl`: logs and stores the fluid-volume results
- `in_docker_organized/main_create_phantom_ionic_chamber.jl`: reports analytical and numerical air-volume comparisons for chamber runs

## Numerical Volume Calculation

The numerical method counts voxels belonging to a binary mask and multiplies by the physical voxel volume.

Formula:

`V_numerical = sum(mask) * (spacing_x * spacing_y * spacing_z)`

Strengths:

- works for arbitrary voxelized geometry
- naturally handles intersections and mask logic
- directly matches the discretized phantom that is written to disk

Limitations:

- sensitive to voxel spacing and discretization
- boundary voxels introduce partial-volume error
- lower resolutions generally increase the mismatch relative to analytical geometry

## Analytical Volume Calculation

The analytical method derives volume from the underlying geometric primitives and correction terms.

Representative formulas used by the project include:

- cylinder: `V = pi * r^2 * h`
- ellipsoid: `V = (4/3) * pi * r_x * r_y * r_z`
- torus: `V = 2 * pi^2 * R * r_xy * r_z`
- spherical cap: `V = pi * h^2 * (3r - h) / 3`

Strengths:

- provides a geometry-based reference value
- highlights discretization or mask-generation problems
- is especially useful for validating parametric phantom construction

Limitations:

- harder to maintain for complex composite shapes
- requires explicit correction terms for tilted fluids, meniscus geometry, rounded bottoms, and occluding internal objects

## Can Phantom Calculation Strategy

The can workflow uses `compute_accurate_fluid_volume_fixed(...)` for the most detailed comparison path.

That logic accounts for:

1. base fluid volume inside the main can body
2. tilted fluid-surface corrections
3. meniscus geometry near the walls
4. flat-bottom or rounded-bottom lower geometry
5. displaced volume from internal objects such as pipes and balls

The final values are printed during execution and stored alongside the generated outputs.

## Ionic Chamber Calculation Strategy

The ionic chamber workflow reports analytical and numerical air-volume estimates for the internal air cavity. This is a useful sanity check for the chamber geometry, especially when changing layer thicknesses or standardized size variants.

## Interpreting The Difference

The main comparison metric is:

`abs(V_analytical - V_numerical) / V_numerical`

Use the result as a validation signal rather than a hard physics guarantee.

- small differences usually indicate a healthy voxelization of the intended geometry
- moderate differences may be expected for low-resolution runs or highly curved boundaries
- large differences usually justify checking spacing, masks, tilted cuts, and rounded-bottom geometry

## Practical Guidance

- use higher resolutions when you need tighter numerical and analytical agreement
- keep JSON-driven regression cases in `tests/configs/` for stable comparisons
- inspect the generated masks whenever a volume difference changes unexpectedly
