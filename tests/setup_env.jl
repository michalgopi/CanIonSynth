using Pkg

println("Instantiating project environment...")
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Ensure PyCall uses the intended Python interpreter.
if haskey(ENV, "PYTHON") && !isempty(ENV["PYTHON"]) && isfile(ENV["PYTHON"])
    println("Using explicit PYTHON for PyCall: $(ENV["PYTHON"])")
    Pkg.build("PyCall")
elseif haskey(ENV, "VIRTUAL_ENV")
    venv_python = Sys.iswindows() ? joinpath(ENV["VIRTUAL_ENV"], "Scripts", "python.exe") : joinpath(ENV["VIRTUAL_ENV"], "bin", "python")
    if isfile(venv_python)
        println("Found active virtual environment. Configuring PyCall to use it: $venv_python")
        ENV["PYTHON"] = venv_python
        Pkg.build("PyCall")
    else
        println("VIRTUAL_ENV is set but no Python executable was found there. PyCall will keep its current configuration.")
    end
else
    println("No explicit PYTHON or VIRTUAL_ENV found. PyCall will use its current Python configuration.")
end

println("Environment instantiated.")
