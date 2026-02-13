using Pkg
using UUIDs
using Test

# Ensure setup environment script is run
println("Skipping setup_env.jl (already setup manually)...")
# include("setup_env.jl")

const START_DIR = pwd()
const SCRIPT_DIR = "in_docker_organized"

function run_script_and_verify(script_name, args, check_files, description)
    println("\n---------------------------------------------------")
    println("Testing $description")
    println("Script: $script_name")
    println("Args: $args")

    # We need to run from the directory where the script is to ensure imports work
    cd(SCRIPT_DIR)

    # Construct command
    # We use a UUID for the run
    run_uuid = string(UUIDs.uuid4())
    # Replace placeholder <uuid> in args
    final_args = replace.(args, "<uuid>" => run_uuid)
    # Ensure JSON paths are absolute or relative to SCRIPT_DIR
    # If json path starts with ../tests/configs, it should be correct relative to in_docker_organized/

    cmd = `julia $script_name $final_args`

    # Set environment variables to skip upload/wandb
    env = copy(ENV)
    env["SKIP_UPLOAD"] = "true"
    env["SKIP_WANDB"] = "true"

    # Capture output to find the output directory
    output_buffer = IOBuffer()

    try
        # Run the command and capture output
        # We use pipeline to capture stdout
        process = run(pipeline(setenv(cmd, env), stdout=output_buffer, stderr=output_buffer))

        output_str = String(take!(output_buffer))
        # println("Script Output:\n", output_str)

        # Parse output to find the output directory
        # Look for "Output stored in: ..."
        m = match(r"Output stored in: (.+)", output_str)
        if m !== nothing
            output_dir = strip(m.captures[1])
            println("Found output directory: '$output_dir'")
            flush(stdout)

            # Verify files exist
            all_files_exist = true
            for file in check_files
                file_path = joinpath(output_dir, file)
                if isfile(file_path)
                    println("[PASS] File exists: $file")
                else
                    println("[FAIL] File missing: $file")
                    println("Checked path: '$file_path'")
                    all_files_exist = false
                end
                flush(stdout)
            end

            if all_files_exist
                println("SUCCESS: All expected files found.")
                flush(stdout)

                # Cleanup
                println("Cleaning up output directory...")
                rm(output_dir, recursive=true, force=true)

            else
                println("Contents of output directory $output_dir:")
                try
                    foreach(println, readdir(output_dir))
                catch e
                    println("Could not read directory: $e")
                end
                flush(stdout)
                error("FAILURE: Some expected files were missing.")
            end

        else
            println("WARNING: Could not find 'Output stored in:' message in output.")
            println("Cannot verify output files location.")
            println("Full output for debugging:")
            println(output_str)
            error("Verification failed: Could not locate output directory.")
        end

    catch e
        println("Error running script: $e")
        if position(output_buffer) > 0
             println("Partial Output:\n", String(take!(output_buffer)))
        end
        rethrow(e)
    finally
        cd(START_DIR)
    end
end

@testset "Synthetic Data Generation Tests" begin

    # --- CAN PHANTOM TESTS ---

    # 1. Defaults (Random)
    # Args: dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]
    run_script_and_verify(
        "main_create_phantom_can.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0"],
        ["example_can.nii.gz", "fluid_mask.nii.gz", "argss.json"],
        "Can Phantom - Default (Random)"
    )

    # 2. With Radon Transform
    run_script_and_verify(
        "main_create_phantom_can.jl",
        ["32x32x32", "true", "false", "<uuid>", "false", "false", "0.0"],
        ["example_can.nii.gz", "fluid_mask.nii.gz", "after_radon.nii.gz", "after_radon_plus_before.nii.gz"],
        "Can Phantom - With Radon Transform"
    )

    # 3. With Smoothing and Noise
    run_script_and_verify(
        "main_create_phantom_can.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "true", "0.1"],
        ["example_can.nii.gz"],
        "Can Phantom - With Smoothing and Noise"
    )

    # 4. From JSON Configuration
    # Note: path relative to in_docker_organized/
    run_script_and_verify(
        "main_create_phantom_can.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0", "../tests/configs/can_phantom.json"],
        ["example_can.nii.gz", "argss.json"],
        "Can Phantom - From JSON Config"
    )


    # --- IONIC CHAMBER TESTS ---

    # 1. Defaults (Random)
    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0"],
        ["ionic_chamber.nii.gz", "ionic_chamber_params.json"],
        "Ionic Chamber - Default (Random)"
    )

    # 2. With Radon Transform
    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "true", "false", "<uuid>", "false", "false", "0.0"],
        ["ionic_chamber.nii.gz", "after_radon.nii.gz"],
        "Ionic Chamber - With Radon Transform"
    )

    # 3. Variable Spacing
    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "false", "true", "<uuid>", "false", "false", "0.0"],
        ["ionic_chamber.nii.gz"],
        "Ionic Chamber - Variable Spacing"
    )

    # 4. From JSON Configuration
    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0", "../tests/configs/ionic_chamber.json"],
        ["ionic_chamber.nii.gz", "ionic_chamber_params.json"],
        "Ionic Chamber - From JSON Config"
    )

end
