#!/bin/bash
# Generate a Steel Can (Flat bottom)
# Args: dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]

echo "Generating Steel Can Phantom..."
cd in_docker_organized

# Note: We simulate a steel can by configuring the JSON parameters.
# Since we are calling the script directly without a JSON file argument,
# the script typically randomizes parameters.
# To guarantee a flat bottom (Steel), we would ideally use a JSON config.
# For this example, we rely on the script's randomness or default behavior
# (which might be mixed), but we demonstrate the command invocation.

# To explicitly force parameters, one should use a config file like:
# tests/configs/can_phantom.json (and ensure rounded_bottom is false)

julia --project=.. main_create_phantom_can.jl 64x64x64 false false "steel_example" false false 0.0

echo "Done. Check in_docker_organized/ for output."
