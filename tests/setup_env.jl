using Pkg

println("Instantiating project environment...")
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()
println("Environment instantiated.")
