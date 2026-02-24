import os
import uuid
import subprocess
import json

# Example of batch generation orchestration
# This script mimics what coordinate_phantom_create.py does

NUM_PHANTOMS = 2
OUTPUT_DIR = "batch_output"

def generate_batch():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Generating {NUM_PHANTOMS} phantoms...")

    for i in range(NUM_PHANTOMS):
        run_uuid = str(uuid.uuid4())
        print(f"  Generating phantom {i+1}/{NUM_PHANTOMS} (UUID: {run_uuid})")

        # Command to run Julia script
        # We assume we are running from the repo root
        cmd = [
            "julia",
            "--project=.",
            "in_docker_organized/main_create_phantom_can.jl",
            "64x64x64",  # Dims
            "false",     # Radon
            "false",     # Variable spacing
            run_uuid,    # UUID
            "true",      # Randomize parameters
            "false",     # Smooth
            "0.0"        # Noise
        ]

        # Set environment to skip upload for local test
        env = os.environ.copy()
        env["SKIP_UPLOAD"] = "true"
        env["SKIP_WANDB"] = "true"

        try:
            subprocess.run(cmd, env=env, check=True, capture_output=True)
            print(f"  Success: {run_uuid}")
        except subprocess.CalledProcessError as e:
            print(f"  Error generating {run_uuid}: {e}")
            print(e.stderr.decode())

if __name__ == "__main__":
    generate_batch()
