import os
import uuid
import subprocess
import json
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

# Example of batch generation orchestration
# This script mimics what coordinate_phantom_create.py does

NUM_PHANTOMS = 2
OUTPUT_DIR = "batch_output"

def generate_batch():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    log.info(f"Generating {NUM_PHANTOMS} phantoms...")

    for i in range(NUM_PHANTOMS):
        run_uuid = str(uuid.uuid4())
        log.info(f"Generating phantom {i+1}/{NUM_PHANTOMS} (UUID: {run_uuid})")

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

        try:
            subprocess.run(cmd, check=True, capture_output=True)
            log.info(f"Success: {run_uuid}")
        except subprocess.CalledProcessError as e:
            log.error(f"Error generating {run_uuid}: {e}")
            log.error(e.stderr.decode())

if __name__ == "__main__":
    generate_batch()
