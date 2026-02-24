#!/usr/bin/env python3
import os
import subprocess
import shutil
import json
import re

# Configuration
RESOLUTION = "128x128x128"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUTPUT_BASE_DIR = os.path.join(PROJECT_ROOT, "manual_test_outputs")
IN_DOCKER_DIR = os.path.join(PROJECT_ROOT, "in_docker_organized")
VISUALIZE_SCRIPT = os.path.join(PROJECT_ROOT, "scripts", "visualize_nifti.py")

# Environment variables
OS_ENV = os.environ.copy()
OS_ENV["SKIP_UPLOAD"] = "true"
OS_ENV["SKIP_WANDB"] = "true"

def run_command(cmd, cwd=None):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, env=OS_ENV, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        return None
    return result.stdout

def extract_output_dir(output_log):
    match = re.search(r"Output stored in:\s*(\S+)", output_log)
    if match:
        return match.group(1)
    return None

def process_test_case(name, script, args, config_json=None):
    print(f"\n--- Processing Test Case: {name} ---")
    
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
        print(f"Failed to generate {name}")
        return

    tmp_dir = extract_output_dir(output_log)
    if not tmp_dir:
        print(f"Could not find output directory for {name}")
        # Try to find it in the log if "Output stored in" is not there
        print("Full log output:")
        print(output_log)
        return

    print(f"Found temp output directory: {tmp_dir}")
    
    # Identify NIfTI files to copy
    nifti_files = [f for f in os.listdir(tmp_dir) if f.endswith(".nii.gz")]
    if not nifti_files:
        print(f"No NIfTI files found in {tmp_dir}")
        return

    for nf in nifti_files:
        src = os.path.join(tmp_dir, nf)
        dest_name = f"{name}_{nf}"
        dest = os.path.join(OUTPUT_BASE_DIR, dest_name)
        shutil.copy2(src, dest)
        print(f"Copied {nf} to {dest}")
        
        # Run visualization
        vis_cmd = ["python3", VISUALIZE_SCRIPT, dest]
        run_command(vis_cmd)

def main():
    if not os.path.exists(OUTPUT_BASE_DIR):
        os.makedirs(OUTPUT_BASE_DIR)
        print(f"Created output directory: {OUTPUT_BASE_DIR}")
    else:
        print(f"Using existing output directory: {OUTPUT_BASE_DIR}")

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

    print("\n--- All tests completed ---")
    print(f"Results are available in: {OUTPUT_BASE_DIR}")

if __name__ == "__main__":
    main()
