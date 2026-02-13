using Pkg
using UUIDs
using Test

# Ensure setup environment script is run
println("Skipping setup_env.jl (already setup manually)...")
# include("setup_env.jl")

const START_DIR = pwd()
const SCRIPT_DIR = "in_docker_organized"

function run_script_and_verify(script_name, args, check_files)
    println("\n---------------------------------------------------")
    println("Testing $script_name with args: $args")

    # We need to run from the directory where the script is to ensure imports work
    cd(SCRIPT_DIR)

    # Construct command
    # We use a UUID for the run
    run_uuid = string(UUIDs.uuid4())
    # Replace placeholder <uuid> in args
    final_args = replace.(args, "<uuid>" => run_uuid)

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
        println("Script Output:\n", output_str)

        # Parse output to find the output directory
        # Look for "Output stored in: ..."
        m = match(r"Output stored in: (.+)", output_str)
        if m !== nothing
            output_dir = strip(m.captures[1])
            open("debug_test_run.log", "a") do io
                println(io, "Found output directory: '$output_dir'")
            end
            println("Found output directory: '$output_dir'")
            flush(stdout)

            # Verify files exist
            all_files_exist = true
            for file in check_files
                file_path = joinpath(output_dir, file)
                if isfile(file_path)
                    open("debug_test_run.log", "a") do io
                        println(io, "[PASS] File exists: $file")
                    end
                    println("[PASS] File exists: $file")
                else
                    open("debug_test_run.log", "a") do io
                        println(io, "[FAIL] File missing: $file")
                        println(io, "Checked path: '$file_path'")
                    end
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
                open("debug_test_run.log", "a") do io
                    println(io, "Contents of output directory $output_dir:")
                    try
                        foreach(line -> println(io, line), readdir(output_dir))
                    catch e
                        println(io, "Could not read directory: $e")
                    end
                end

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
    # Test Can Phantom
    # Args: dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]
    # Example: 32x32x32 false false <uuid> false false 0.0
    run_script_and_verify(
        "main_create_phantom_can.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0"],
        ["example_can.nii.gz", "fluid_mask.nii.gz", "argss.json"]
    )

    # Test Ionic Chamber Phantom
    # Args: dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]
    # Example: 32x32x32 false false <uuid> false false 0.0
    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0"],
        ["ionic_chamber.nii.gz", "ionic_chamber_params.json"]
    )
end
