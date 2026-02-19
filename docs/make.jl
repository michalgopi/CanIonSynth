using Documenter

# We need to make sure the package is available (even if it's not a proper package structure yet, we can just document the scripts)
# Ideally, we would have a module. For now, we will just generate pages from markdown.

makedocs(
    sitename = "Synthetic Tomography",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "user_guide.md",
        "Workflow Summary" => "workflow_summary.md",
        "Functions Reference" => "functions_reference.md",
        "Cleanup Proposals" => "cleanup_proposals.md"
    ]
)

deploydocs(
    repo = "github.com/jakubMitura14/SyntheticTomo.git", # This should be updated by the user to their repo
)
