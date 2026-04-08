# Phantom Types

The repository generates two phantom families: cans and ionic chambers.

## Can Phantoms

Can phantoms model filled industrial containers with internal structures and fluid interfaces.

### Shared Can Features

- cylindrical outer wall with configurable thickness
- configurable top and bottom curvature
- one or two fluid phases with different densities
- optional meniscus geometry and tilted fluid surface
- optional internal objects such as balls, pipe structures, and dispenser components

### Flat-Bottom Can

Use `rounded_bottom=false` to create a can with a flatter base profile. This variant is useful when you want a sharper wall-to-bottom transition.

### Rounded-Bottom Can

Use `rounded_bottom=true` to create a can with a curved base and smoother lower transition. The checked-in `tests/configs/can_phantom.json` follows this family.

### Important Can Parameters

- `bigger_cyl_size`: main can dimensions
- `cylinder_wall_thickness`: wall thickness
- `density_inside` and `density_inside_b`: fluid densities
- `dual_phase_percentage`: relative fluid split when using two phases
- `x_cut_angle` and `y_cut_angle`: tilted fluid surface controls
- `menisc_radius` and related meniscus parameters: fluid curvature controls
- `add_pipe`, `first_ball`, `second_ball`: optional internal features

## Ionic Chamber Phantoms

Ionic chamber phantoms model layered radiation-measurement devices with multiple materials and head shapes.

### Shared Chamber Structure

Most chamber variants are built from concentric or overlapping cylindrical elements representing:

- graphite housing
- aluminum shielding layers
- insulation layers
- graphite and copper electrode components
- internal air cavity

### Supported Chamber Variants

- `rounded_top=true`: rounded cylindrical chamber
- `square_top=true`: flat-top chamber
- `ball_like=true`: ball-like head shape
- `lolipop_like=true`: lollipop-style head with narrow neck

The repository already contains reproducible selectors for several of these shapes:

- `tests/configs/ionic_chamber_square.json`
- `tests/configs/ionic_chamber_ball.json`
- `tests/configs/ionic_chamber_lolipop.json`

### Standardized Sizes

When `new_flat_sizes=true`, the generator can produce standardized chamber dimensions controlled by `rand_ver`.

- `rand_ver=1`: largest standard size family
- `rand_ver=2`: medium standard size family
- `rand_ver=3`: smallest standard size family

### Important Chamber Parameters

- `total_len`, `base_len`, `main_radius`: main geometry dimensions
- `air_thickness`: air cavity thickness
- `graphite_density`, `copper_density`, `aluminium_density`, `insulation_density`: material densities
- `graphite_electrode_radius`, `copper_radius`: electrode sizes
- `inner_insluation_thickness`, `outer_insluation_thickness`: insulation thicknesses
- `aluminium_inner_thicness`, `aluminium_outer_thicness`: shield thicknesses
- `add_graphite_in_copper`, `elongate_copper`, `add_spiral`: optional structure variations

## Choosing A Phantom Type

Use can phantoms when you need fluid interfaces, meniscus behavior, or internal inclusions inside a container. Use ionic chamber phantoms when you need layered multi-material geometry with well-defined internal air cavities and controlled electrode structure.
