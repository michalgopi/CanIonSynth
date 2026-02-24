#!/bin/bash
# Generate a Steel Can (Flat bottom)

echo "Generating Steel Can Phantom..."
cd in_docker_organized

# Note: We simulate a steel can by configuring the JSON parameters.
echo '{"rounded_bottom": false, "cylinder_wall_thickness": 0.02}' > temp_steel.json

# Run Julia and capture output to find the directory
# We use 'tee' to show output to user while capturing
OUTPUT_LOG=$(julia --project=.. main_create_phantom_can.jl 64x64x64 false false "steel_example" false false 0.0 temp_steel.json)
echo "$OUTPUT_LOG"

# Extract output directory
OUTPUT_DIR=$(echo "$OUTPUT_LOG" | grep "Output stored in:" | awk '{print $4}')

rm temp_steel.json

if [ -n "$OUTPUT_DIR" ]; then
    echo ""
    echo "generation complete."
    echo "Output directory: $OUTPUT_DIR"
    echo "To visualize:"
    echo "  python3 ../scripts/visualize_nifti.py $OUTPUT_DIR/example_can.nii.gz"
else
    echo "Could not determine output directory."
fi
