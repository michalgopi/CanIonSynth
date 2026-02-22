#=
cylinder.jl
=#


export Cylinder_min_half_elipsoid, cylinder_min_half_elipsoid


"""
    Cylinder <: AbstractShape{3}
"""
struct Cylinder_min_half_elipsoid <: AbstractShape{3} end


# constructor


"""
    cylinder(cx, cy, cz, wx, wy, wz, Φ, Θ, value::Number)
    cylinder(center::NTuple{3,RealU}, width::NTuple{3,RealU}, angle::NTuple{3,RealU}, v)
Construct `Object{Cylinder}` from parameters;
here `width` is the *radius* in x,y and the *height* in z.
"""
cylinder_min_half_elipsoid(args... ; kwargs...) = Object(Cylinder_min_half_elipsoid(), args...; kwargs...)


# methods
volume1(::Cylinder_min_half_elipsoid) = π-(2/3 * π)
ℓmax1(::Cylinder_min_half_elipsoid) = √5 # max line integral through unit cylinder
ℓmax(ob::Object3d{Cylinder_min_half_elipsoid}) = sqrt(sum(abs2, (2maximum(ob.width[1:2]), ob.width[3])))

"""
    phantom1(ob::Object3d{Cylinder_min_half_elipsoid}, (x,y,z))
Evaluate unit cylinder at `(x,y,z)`,
for unitless coordinates.
"""
phantom1(ob::Object3d{Cylinder_min_half_elipsoid}, xyz::NTuple{3,Real}) =
    (((abs(xyz[3]) ≤ 0.5) && (sum(abs2, xyz[1:2]) ≤ 1)) && (!((sum(abs2, xyz) ≤ 1) && (xyz[1] >= 0))) )

# line integral through rectangle of width (wx,wy) at (r,ϕ)

function get_ray_half_sphere(u,v,ϕ,θ,axis)
    T = Float32
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

    # if(axis == 1)
    #     if(x<0)
    #         return zero(T)
    #     end
    # elseif(axis == 2)
    #     if(y<0)
    #         return zero(T)
    #     end
    # elseif(axis == 3)
    #     if(z<0)
    #         return zero(T)
    #     end
    # end

    if r2 < 1
        return  sqrt(one(T) - r2)
    end

    return zero(T)
end

function get_xyz(u,v,ϕ,θ,axis)
    T = Float32
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

    return x,y,z
end


# x-ray transform (line integral) of unit cylinder
# `u,v` should be unitless
function xray1(
    ::Cylinder_min_half_elipsoid,
    u::Real,
    v::Real,
    ϕ::RealU, # azim (irrelevant)
    θ::RealU, # polar
)
    T = promote_type(typeof(u), typeof(v), Float32)

    r = abs(u)
    if r > 1
        return zero(T)
    end
    # rectangle in plane of distance `r` from origin
    wz = 1 # from unit-height of cylinder
    wy = 2 * sqrt(1 - r^2)

    half_sphere_axis=2
    res= T(_rect_proj(wz, wy, v, θ))
    # res=res-(get_ray_half_sphere(u,v,ϕ,θ,half_sphere_axis)*1.4)
    x,y,z=get_xyz(u,v,ϕ,θ,half_sphere_axis)


    if(z<0)
        return zero(T)
    end

    if(res<0)
        return zero(T)
    end

    return res
end
