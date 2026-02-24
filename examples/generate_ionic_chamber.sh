#!/bin/bash
# Generate an Ionic Chamber (Ball shaped)

echo "Generating Ionic Chamber Phantom..."
cd in_docker_organized

# Create config for ball shape
cat <<EOF > temp_ionic_config.json
{
    "ball_like": true,
    "lollipop_like": false,
    "new_flat_sizes": false
}
EOF

julia --project=.. main_create_phantom_ionic_chamber.jl 64x64x64 false false "ionic_example" false false 0.0 temp_ionic_config.json

rm temp_ionic_config.json
echo "Done. Check in_docker_organized/ for output."
