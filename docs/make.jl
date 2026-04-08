using Documenter

# We need to make sure the package is available (even if it's not a proper package structure yet, we can just document the scripts)
# Ideally, we would have a module. For now, we will just generate pages from markdown.

makedocs(
    sitename = "Synthetic Tomography",
    source = @__DIR__,
    build = joinpath(@__DIR__, "..", "docs-build"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        inventory_version = "dev"
    ),
    pages = [
        "Home" => "index.md",
        "User Guide" => "user_guide.md",
        "Workflow Summary" => "workflow_summary.md",
        "System Architecture" => "system_architecture.md",
        "Phantom Types" => "phantom_types.md",
        "Functions Reference" => "functions_reference.md",
        "Volume Calculation" => "volume_calculation.md",
        "Manual Testing" => "manual_testing_visualization.md",
        "Manual Verification" => "manual_verification.md"
    ]
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "github.com/jakubMitura14/synth_data_generation.git",
    )
end
