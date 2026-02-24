#!/bin/bash
# Generate an Aluminum Can (Rounded bottom)

echo "Generating Aluminum Can Phantom..."
cd in_docker_organized

# Create config
cat <<EOF > temp_alum_config.json
{
    "rounded_bottom": true,
    "cylinder_wall_thickness": 0.03,
    "two_fluids": true
}
EOF

# Run Julia
OUTPUT_LOG=$(julia --project=.. main_create_phantom_can.jl 64x64x64 false false "alum_example" false false 0.0 temp_alum_config.json)
echo "$OUTPUT_LOG"

# Extract output directory
OUTPUT_DIR=$(echo "$OUTPUT_LOG" | grep "Output stored in:" | awk '{print $4}')

rm temp_alum_config.json

if [ -n "$OUTPUT_DIR" ]; then
    echo ""
    echo "Generation complete."
    echo "Output directory: $OUTPUT_DIR"
    echo "To visualize:"
    echo "  python3 ../scripts/visualize_nifti.py $OUTPUT_DIR/example_can.nii.gz"
else
    echo "Could not determine output directory."
fi
