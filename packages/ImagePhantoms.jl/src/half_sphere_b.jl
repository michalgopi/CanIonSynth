# src/half_sphere.jl
#=
half_sphere.jl
=#
export HalfSphere_b_x, half_sphere_x_b,HalfSphere_b_y, half_sphere_y_b,HalfSphere_b_z, half_sphere_z_b

"""
    HalfSphere_b <: AbstractShape{3}
"""
struct HalfSphere_b_x <: AbstractShape{3}end
struct HalfSphere_b_y <: AbstractShape{3}end
struct HalfSphere_b_z <: AbstractShape{3}end
# constructor
"""
    half_sphere(cx, cy, cz, r, Φ=0, Θ=0, value::Number = 1)
    half_sphere(center::NTuple{3,RealU}, radius::RealU, angle::NTuple{3,RealU}, v)
Construct `Object{HalfSphere_b}` from parameters.
"""
half_sphere_x_b(args... ; kwargs...) = Object(HalfSphere_b_x(), args...; kwargs...)
half_sphere_y_b(args... ; kwargs...) = Object(HalfSphere_b_y(), args...; kwargs...)
half_sphere_z_b(args... ; kwargs...) = Object(HalfSphere_b_z(), args...; kwargs...)

# methods
volume1(::HalfSphere_b_x) = 2/3 * π # volume of unit half-sphere
volume1(::HalfSphere_b_y) = 2/3 * π # volume of unit half-sphere
volume1(::HalfSphere_b_z) = 2/3 * π # volume of unit half-sphere

ℓmax1(::HalfSphere_b_x) = 2 # maximum chord through a unit half-sphere
ℓmax1(::HalfSphere_b_y) = 2 # maximum chord through a unit half-sphere
ℓmax1(::HalfSphere_b_z) = 2 # maximum chord through a unit half-sphere

"""
evaluates weather a point is inside the unit half-sphere
first part is identical as in unit sphere second part is added to consider only the upper half-sphere
"""
phantom1(ob::Object3d{HalfSphere_b_x}, xyz::NTuple{3,Real}) = (sum(abs2, xyz) ≤ 1) && (xyz[1] <= 0)
phantom1(ob::Object3d{HalfSphere_b_y}, xyz::NTuple{3,Real}) = (sum(abs2, xyz) ≤ 1) && (xyz[2] <= 0)
phantom1(ob::Object3d{HalfSphere_b_z}, xyz::NTuple{3,Real}) = (sum(abs2, xyz) ≤ 1) && (xyz[3] <= 0)

# x-ray transform (line integral) of unit half-sphere
# `u,v` should be unitless
function xray1_main_b(
    u::Ru,
    v::Rv,
    ϕ::RealU,
    θ::RealU, #
    axis::Int
) where {Ru <: Real, Rv <: Real}
    T = promote_type(Ru, Rv, Float32)
    r2 = u^2 + v^2

    (sϕ, cϕ) = sincos(ϕ)
    (sθ, cθ) = sincos(θ)

    p1 = u * cϕ + v * sϕ * sθ
    p2 = u * sϕ - v * cϕ * sθ
    p3 = v * cθ

    e1 = -sϕ * cθ # x = p1 + ℓ * e1
    e2 = cϕ * cθ  # y = p2 + ℓ * e2
    e3 = sθ       # z = p3 + ℓ * e3

    ℓxmin, ℓxmax = cube_bounds(p1, e1)
    ℓymin, ℓymax = cube_bounds(p2, e2)
    ℓzmin, ℓzmax = cube_bounds(p3, e3)

    minℓ = max(ℓxmin, ℓymin, ℓzmin)
    maxℓ = min(ℓxmax, ℓymax, ℓzmax)
    ℓ = max(maxℓ - minℓ, zero(T))
    x = p1 + ℓ * e1
    y = p2 + ℓ * e2
    z = p3 + ℓ * e3

    if(axis == 1)
        if(x>0)
            return zero(T)
        end
    elseif(axis == 2)
        if(y>0)
            return zero(T)
        end
    elseif(axis == 3)
        if(z>0)
            return zero(T)
        end
    end
    if r2 < 1
        return  sqrt(one(T) - r2)
    end

    return zero(T)
end

# volume1(::HalfSphere_b_x) = 2/3 * π # volume of unit half-sphere
# volume1(::HalfSphere_b_y) = 2/3 * π # volume of unit half-sphere
# volume1(::HalfSphere_b_z) = 2/3 * π # volume of unit half-sphere

function xray1(
    hs::HalfSphere_b_x,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
return xray1_main_b( u,v,ϕ,θ,1)
end

function xray1(
    hs::HalfSphere_b_y,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
return xray1_main_b( u,v,ϕ,θ,2)
end

function xray1(
    hs::HalfSphere_b_z,
    u::Ru,
    v::Rv,
    ϕ::RealU, # irrelevant
    θ::RealU, # irrelevant
) where {Ru <: Real, Rv <: Real}
return xray1_main_b( u,v,ϕ,θ,3)
end
