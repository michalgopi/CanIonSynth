using Pkg

println("Installing dependencies...")

# Add standard packages
packages = [
    "Sinograms",
    "ImageGeoms",
    "MIRTjim",
    "Unitful",
    "HDF5",
    "PyCall",
    "LazyGrids",
    "ImageFiltering",
    "Accessors",
    "FFTW",
    "Plots",
    "Revise",
    "Statistics",
    "UUIDs",
    "JSON",
    "Dates",
    "Logging"
]

for pkg in packages
    try
        Pkg.add(pkg)
        println("Installed $pkg")
    catch e
        println("Error installing $pkg: $e")
    end
end

# Add ImagePhantoms from URL
println("Installing ImagePhantoms from GitHub...")
try
    Pkg.add(url="https://github.com/jakubMitura14/ImagePhantoms.jl.git")
    println("Installed ImagePhantoms")
catch e
    println("Error installing ImagePhantoms: $e")
end

println("Dependencies installation complete.")
