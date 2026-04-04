using PyCall
try
    sitk = pyimport("SimpleITK")
    println("SimpleITK version: ", sitk.__version__)
catch e
    println("Failed to import SimpleITK: ", e)
end
