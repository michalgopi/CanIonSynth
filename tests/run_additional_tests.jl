using Pkg
using UUIDs
using Test

const START_DIR = pwd()
const SCRIPT_DIR = "in_docker_organized"

function run_script_and_verify(script_name, args, description)
    println("\n---------------------------------------------------")
    println("Testing $description")

    cd(SCRIPT_DIR)

    run_uuid = string(UUIDs.uuid4())
    # Replace placeholder <uuid> in args
    final_args = replace.(args, "<uuid>" => run_uuid)

    cmd = `julia --project=$START_DIR $script_name $final_args`

    env = copy(ENV)
    env["SKIP_UPLOAD"] = "true"
    env["SKIP_WANDB"] = "true"

    try
        run(setenv(cmd, env))
        println("SUCCESS: $description")
    catch e
        println("FAILURE: $description - $e")
        rethrow(e)
    finally
        # Cleanup
        if isdir(run_uuid)
            rm(run_uuid, recursive=true, force=true)
        end
        cd(START_DIR)
    end
end

@testset "Extended Phantom Tests" begin

    # Test standardized sizes for Ionic Chamber (new_flat_sizes=true)
    # We need to construct a JSON config to pass this parameter because it might not be a positional arg
    # or we can rely on randomness if we can't set it via args.
    # Looking at main_create_phantom_ionic_chamber.jl args:
    # dims add_radon variable_spacing uuid randomize add_smooth additive_noise [json_path]

    # We will create a temp json file
    config_path = joinpath(SCRIPT_DIR, "temp_test_config.json")
    open(config_path, "w") do f
        println(f, "{\"new_flat_sizes\": true, \"rand_ver\": 1}")
    end

    run_script_and_verify(
        "main_create_phantom_ionic_chamber.jl",
        ["32x32x32", "false", "false", "<uuid>", "false", "false", "0.0", "temp_test_config.json"],
        "Ionic Chamber - Standardized Size (Ver 1)"
    )

    rm(config_path; force=true)

end
