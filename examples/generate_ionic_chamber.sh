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

# Run Julia
OUTPUT_LOG=$(julia --project=.. main_create_phantom_ionic_chamber.jl 64x64x64 false false "ionic_example" false false 0.0 temp_ionic_config.json)
echo "$OUTPUT_LOG"

# Extract output directory
OUTPUT_DIR=$(echo "$OUTPUT_LOG" | grep "Output stored in:" | awk '{print $4}')

rm temp_ionic_config.json

if [ -n "$OUTPUT_DIR" ]; then
    echo ""
    echo "Generation complete."
    echo "Output directory: $OUTPUT_DIR"
    echo "To visualize:"
    echo "  python3 ../scripts/visualize_nifti.py $OUTPUT_DIR/ionic_chamber.nii.gz"
else
    echo "Could not determine output directory."
fi
