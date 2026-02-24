#!/bin/bash
# Generate an Aluminum Can (Rounded bottom)

echo "Generating Aluminum Can Phantom..."
cd in_docker_organized

# Using a config file is the best way to ensure specific geometry.
# We will create a temporary config for this example.

cat <<EOF > temp_alum_config.json
{
    "rounded_bottom": true,
    "cylinder_wall_thickness": 0.03,
    "two_fluids": true
}
EOF

julia --project=.. main_create_phantom_can.jl 64x64x64 false false "alum_example" false false 0.0 temp_alum_config.json

rm temp_alum_config.json
echo "Done. Check in_docker_organized/ for output."
