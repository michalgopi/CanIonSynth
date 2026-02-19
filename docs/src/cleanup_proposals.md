# Proposed Repository Cleanup

Based on an analysis of the repository structure and file contents, the following cleanup actions are recommended to improve organization and remove redundant or obsolete files.

**Note:** This is a proposal only. No files have been deleted.

## Files Recommended for Deletion

### 1. Root Directory Redundancies
The following files in the root directory appear to be duplicates or older versions of files now better organized in `in_docker_organized/`.

*   **`get_geometry_main.jl`**
    *   **Reason:** This file is a duplicate of `in_docker_organized/get_geometry_main.jl`. The version in `in_docker_organized/` is the one referenced in the documentation and seems to be the active one.
    *   **Action:** Delete.

*   **`generate-random_can_low_res.jl`**
    *   **Reason:** This script appears to be an older or alternative version of `in_docker_organized/main_create_phantom_can.jl`. It contains hardcoded API keys and less structured logic compared to the organized version.
    *   **Action:** Delete (ensure any unique logic is ported to `main_create_phantom_can.jl` if needed, though a quick check suggests `main_create_phantom_can.jl` is superior).

### 2. Temporary and Backup Files
These files appear to be artifacts of development or backup copies that are no longer needed.

*   **`Dockerfile copy`**
    *   **Reason:** A backup copy of a Dockerfile.
    *   **Action:** Delete.

*   **`oldDockerfileOld.txt`**
    *   **Reason:** clearly an old backup file ("old...Old.txt").
    *   **Action:** Delete.

*   **`test.py`**
    *   **Reason:** Contains only `import numpy as np`. It serves no functional purpose.
    *   **Action:** Delete.

*   **`download_data.py`**
    *   **Reason:** Scripts for downloading "Task09_Spleen" dataset (MONAI), which seems unrelated to the synthetic geometric phantom generation which is the core purpose of this repo.
    *   **Action:** Delete (unless there is a hidden dependency on this external data, which seems unlikely for *synthetic* generation).

## Proposed Restructuring

To further clean up the repository, the following structural changes are proposed:

1.  **Move `CtFanArc_params.jl` to `in_docker_organized/`**
    *   **Reason:** This file is closely related to the projection logic used by the scripts in `in_docker_organized/`. Keeping it in the root while other dependencies are in the subdirectory is inconsistent.
    *   **Action:** Move `CtFanArc_params.jl` to `in_docker_organized/` and update `include` paths in `main_create_phantom_can.jl` and `main_create_phantom_ionic_chamber.jl`.

2.  **Rename `in_docker_organized` to `src`**
    *   **Reason:** `in_docker_organized` is a descriptive but non-standard name. `src` is the conventional directory for source code.
    *   **Action:** Rename directory.

## Summary of Actions

| File/Directory | Recommendation | Reason |
| :--- | :--- | :--- |
| `get_geometry_main.jl` | Delete | Duplicate of `in_docker_organized/get_geometry_main.jl` |
| `generate-random_can_low_res.jl` | Delete | Obsolete/Duplicate script |
| `Dockerfile copy` | Delete | Backup file |
| `oldDockerfileOld.txt` | Delete | Backup file |
| `test.py` | Delete | Empty test file |
| `download_data.py` | Delete | Unrelated to core functionality |
| `CtFanArc_params.jl` | Move | Move to source directory for consistency |
| `in_docker_organized/` | Rename | Rename to `src/` for convention |
