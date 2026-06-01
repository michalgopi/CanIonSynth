#!/usr/bin/env python3
import os
import subprocess
import shutil
import json
import re
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

# Configuration
RESOLUTION = "128x128x128"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTPUT_BASE_DIR = os.path.join(PROJECT_ROOT, "manual_test_outputs")
IN_DOCKER_DIR = os.path.join(PROJECT_ROOT, "in_docker_organized")
VISUALIZE_SCRIPT = os.path.join(PROJECT_ROOT, "scripts", "visualize_nifti.py")

OS_ENV = os.environ.copy()

def run_command(cmd, cwd=None):
    log.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, env=OS_ENV, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"Error running command: {result.stderr}")
        return None
    return result.stdout

def extract_output_dir(output_log):
    match = re.search(r"Output stored in:\s*(\S+)", output_log)
    if match:
        return match.group(1)
    return None

def process_test_case(name, script, args, config_json=None):
    log.info(f"Processing Test Case: {name}")
    
    config_path = None
    if config_json:
        config_path = os.path.join(IN_DOCKER_DIR, f"temp_{name}.json")
        with open(config_path, "w") as f:
            json.dump(config_json, f)
    
    # Prepare command
    cmd = ["julia", "--project=..", script, RESOLUTION] + args
    if config_path:
        cmd.append(os.path.basename(config_path))
    
    output_log = run_command(cmd, cwd=IN_DOCKER_DIR)
    
    if config_path and os.path.exists(config_path):
        os.remove(config_path)
    
    if not output_log:
        log.error(f"Failed to generate {name}")
        return

    tmp_dir = extract_output_dir(output_log)
    if not tmp_dir:
        log.error(f"Could not find output directory for {name}")
        log.info(f"Full log output: {output_log}")
        return

    log.info(f"Found temp output directory: {tmp_dir}")

    # Identify NIfTI files to copy
    nifti_files = [f for f in os.listdir(tmp_dir) if f.endswith(".nii.gz")]
    if not nifti_files:
        log.warning(f"No NIfTI files found in {tmp_dir}")
        return

    for nf in nifti_files:
        src = os.path.join(tmp_dir, nf)
        dest_name = f"{name}_{nf}"
        dest = os.path.join(OUTPUT_BASE_DIR, dest_name)
        shutil.copy2(src, dest)
        log.info(f"Copied {nf} to {dest}")
        
        # Run visualization
        vis_cmd = ["python3", VISUALIZE_SCRIPT, dest]
        run_command(vis_cmd)

def main():
    if not os.path.exists(OUTPUT_BASE_DIR):
        os.makedirs(OUTPUT_BASE_DIR)
        log.info(f"Created output directory: {OUTPUT_BASE_DIR}")
    else:
        log.info(f"Using existing output directory: {OUTPUT_BASE_DIR}")

    test_cases = [
        {
            "name": "steel_can",
            "script": "main_create_phantom_can.jl",
            "args": ["false", "false", "test_steel", "false", "false", "0.0"],
            "config": {"rounded_bottom": False, "cylinder_wall_thickness": 0.02}
        },
        {
            "name": "alum_can",
            "script": "main_create_phantom_can.jl",
            "args": ["false", "false", "test_alum", "false", "false", "0.0"],
            "config": {"rounded_bottom": True, "cylinder_wall_thickness": 0.035, "two_fluids": True}
        },
        {
            "name": "radon_can",
            "script": "main_create_phantom_can.jl",
            "args": ["true", "false", "test_radon", "false", "false", "0.0"],
            "config": None
        },
        {
            "name": "ionic_ball",
            "script": "main_create_phantom_ionic_chamber.jl",
            "args": ["false", "false", "test_ball", "false", "false", "0.0"],
            "config": {"ball_like": True, "lollipop_like": False}
        },
        {
            "name": "ionic_lollipop",
            "script": "main_create_phantom_ionic_chamber.jl",
            "args": ["false", "false", "test_lollipop", "false", "false", "0.0"],
            "config": {"lollipop_like": True}
        },
        {
            "name": "ionic_standard",
            "script": "main_create_phantom_ionic_chamber.jl",
            "args": ["false", "false", "test_std", "false", "false", "0.0"],
            "config": {"new_flat_sizes": True, "rand_ver": 1}
        }
    ]

    for tc in test_cases:
        process_test_case(tc["name"], tc["script"], tc["args"], tc["config"])

    log.info("All tests completed")
    log.info(f"Results are available in: {OUTPUT_BASE_DIR}")

if __name__ == "__main__":
    main()
