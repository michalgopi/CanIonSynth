using Pkg

println("Instantiating project environment...")
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Ensure PyCall uses the active Virtual Environment Python
if haskey(ENV, "VIRTUAL_ENV")
    venv_python = Sys.iswindows() ? joinpath(ENV["VIRTUAL_ENV"], "Scripts", "python.exe") : joinpath(ENV["VIRTUAL_ENV"], "bin", "python")
    if isfile(venv_python)
        println("Found active virtual environment. Configuring PyCall to use it: $venv_python")
        ENV["PYTHON"] = venv_python
        Pkg.build("PyCall")
    end
else
    println("No VIRTUAL_ENV found. PyCall will use its default Python environment.")
end

println("Environment instantiated.")
