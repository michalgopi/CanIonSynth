import subprocess
import json
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
log = logging.getLogger(__name__)

# Read the configuration JSON file.
# Set CONFIG_JSON_PATH env var to override the default location.
config_json_path = os.environ.get(
    "CONFIG_JSON_PATH",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "control_json_can_128.json")
)

log.info(f"Reading configuration from {config_json_path}...")
with open(config_json_path, 'r') as f:
    config = json.load(f)

log.info(f"Config: {config}")
dims = config.get("dims", "128x128x128")
add_radon = str(config.get("add_radon", False)).lower()
add_smooth = str(config.get("add_smooth", False)).lower()
additive_noise = str(config.get("additive_noise", 0.10))
randomize = str(config.get("randomize", False)).lower()
variable_spacing = str(config.get("variable_spacing", False)).lower()
num_exec = int(config.get("num_exec", 1))
json_folders_path = config.get("json_folders_path", "").strip()

# Load filenames from local directory if a path is provided
json_files = []

if json_folders_path:
    try:
        log.info(f"Listing files in {json_folders_path}...")
        json_files = sorted([f for f in os.listdir(json_folders_path) if f.endswith(".json")])
        log.info(f"Found {len(json_files)} files in the path")
    except OSError as e:
        log.error(f"Failed to list files in path: {e}")

    json_files = [file for file in json_files if file]
    num_exec = len(json_files)
    log.info(f"num_exec {num_exec}")

log.info(f"json_folders_path: {json_folders_path}, json_files: {json_files}")

script_name = "in_docker_organized/main_create_phantom_ionic_chamber.jl"
if config.get("is_can", False):
    script_name = "in_docker_organized/main_create_phantom_can.jl"

for i in range(num_exec):
    try:
        if json_folders_path:
            file = json_files[i]
            unique_id = dims + '_' + file[:-5]
            json_file_path = os.path.join(json_folders_path, file)

            with open(json_file_path, 'r') as json_file:
                json_data = json.load(json_file)
            log.info(f"json_data: {json_data}")
            log.info(f"json_file_path: {json_file_path}")
            log.info(f"julia {script_name} {dims} {add_radon} {variable_spacing} {unique_id} {randomize} {add_smooth} {additive_noise} {json_file_path}")

            result = subprocess.run(
                ['julia', script_name, dims, add_radon, variable_spacing, unique_id,
                 randomize, add_smooth, additive_noise, json_file_path],
                check=True, capture_output=True, text=True
            )
        else:
            unique_id = f"{dims}_run_{i}"
            log.info(f"julia {script_name} {dims} {add_radon} {variable_spacing} {unique_id} {randomize} {add_smooth} {additive_noise}")

            result = subprocess.run(
                ['julia', script_name, dims, add_radon, variable_spacing, unique_id,
                 randomize, add_smooth, additive_noise],
                check=True, capture_output=True, text=True
            )

    except subprocess.CalledProcessError as e:
        log.error(f"Command execution failed: {e}")
        log.error(f"Error output: {e.stderr}")
    except json.JSONDecodeError as e:
        log.error(f"Failed to parse JSON: {e}")
    except Exception as e:
        log.error(f"An unexpected error occurred: {e}")
    finally:
        log.info(f"Run {i} complete")
