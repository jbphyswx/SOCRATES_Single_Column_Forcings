using ForwardDiff
using PCHIPInterpolation
using Dierckx

function pyinterp(
    x,
    xp,
    fp;
    bc::String = "error",
    k::Int = 1,
    method::Symbol = :Spline1D,
    return_spl::Bool = false,
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing,
    )
    if method == :pchip
        spl = PCHIPInterpolation.Interpolator(xp, fp) # Has a continuous (but not smooth) derivative (piecewise cubic hermite interpolating polynomial, is unique so no k)
        if bc == "error"
            ;
        elseif bc == "extrapolate"
            spl = pchip_extrapolate(spl) # allow extrapolating beyond the edges
        end

    elseif method ∈ [:pchip_smooth_derivative, :pchip_smooth]
        spl = pchip_smooth_derivative(xp, fp; bc=bc, f_enhancement_factor=f_enhancement_factor, f_p_enhancement_factor=f_p_enhancement_factor) # how our current implementation works
    elseif method ∈ [:Spline1D, :Dierckx]
        spl = Dierckx.Spline1D(xp, fp; k = k, bc = bc) # k=1 has disjoint derivative, k=2,3 is smooth in derivative but not monotonicity preserving... this is matlab pchip
    else
        error("method not recognized")
    end
    if return_spl
        return spl
    else
        return spl.(vec(x)) # have to broadcast bc pchip requires it...
    end
end


function pchip_extrapolate(spl::PCHIPInterpolation.Interpolator)
    # extrapolating beyond the bounds of the data, we'll just continue the slope at the edge of the data for continuous derivative (though not necessarily smooth)
    xmin = spl.xs[1]
    xmax = spl.xs[end]
    ymin = spl.ys[1]
    ymax = spl.ys[end]

    return x -> begin
        if x < xmin
            dydx_xmin = PCHIPInterpolation._derivative(spl, Val(:begin), 1) # slope at beginning of interval 1
            return ymin - (xmin - x) * dydx_xmin
        elseif x > xmax
            dydx_xmax = PCHIPInterpolation._derivative(spl, Val(:end), length(spl.xs)-1) # slope at end of interval N-1 (there are only N-1 intervals on N points) (N = length(spl.xs))
            return ymax + (x - xmax) * dydx_xmax
        else
            return spl(x)
        end
    end
end

"""
Create a function with a smooth second derivative by integrating a pchip approximation as the data's first derivative
Ideally you'd fit the spline to f(x) but f(x) may not be smooth. So we fit the spline to f'(x) and then integrate it to get f(x)
This comes at the cost of not necessarily hitting the points in f(x) exactly and rounding corners but it will be smooth...

To reduce the rounding, we can increase the resolution of the data by calculating the derivative spline at more points before creating the spline representation of the derivative
We do this via `enhancement_factor`, which allows us to add `enhancement_factor-1` points between each point in `xp` to constrain the spline of the spline derivative to more closely match the spline derivative and thus f(x) when integrated

We can also reduce rounding by boosting linear interpolation between the given data points, though this is somewhat of a guess as you don't actually know what the function looks like

Essentially:
    f_enhancement_factor: How closely the outcome matches linear interpolation between the data points. How rounded can the curve be?
    f_p_enhancement_factor: How closely the outcome matches the spline fit of f_p (possibly enhanced) -- i.e. does it cut corners or actually go to the points like the spline does?
# To Do: Support bc = "nearest" and bc = "NaN" for nearest neighbor and NaN respectively outside the bounds of xp
"""
function pchip_smooth_derivative(
    xp,
    fp;
    bc::String = "error",
    f_enhancement_factor::Union{Int, Nothing} = nothing,
    f_p_enhancement_factor::Union{Int, Nothing} = nothing)

    if !isnothing(f_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
        # add enhancement_factor-1 points between each point (note this can lead to inexact errors if the interpolated points can't be cast to exterior type (like range needing float but x being in int))
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_enhancement_factor-1))
        for i in 1:length(xp)-1
            xp_new[((i-1) * f_enhancement_factor + 1):(i * f_enhancement_factor)] .= range(xp[i], stop=xp[i+1], length=f_enhancement_factor+1)[1:end-1]
        end
        xp_new[end] = xp[end]
        xp, fp = xp_new, pyinterp(xp_new, xp, fp; method=:Dierckx, return_spl=false)
    end

    # create a pchip interpolator to xp
    f_pchip_spl = PCHIPInterpolation.Interpolator(xp, fp) # should this be the extrapolatory one
    # differentiate it
    dfdx =  ForwardDiff.derivative.(Ref(f_pchip_spl), xp)

    if !isnothing(f_p_enhancement_factor) # increase the resolution of xp, yp by linear interpolation to get more points to constrain the smooth fcn
        # add enhancement_factor-1 points between each point (note this can lead to inexact errors if the interpolated points can't be cast to exterior type (like range needing float but x being in int))
        xp_new = Array{eltype(xp)}(undef, length(xp) + (length(xp) - 1) * (f_p_enhancement_factor-1))
        for i in 1:length(xp)-1
            xp_new[((i-1) * f_p_enhancement_factor + 1):(i * f_p_enhancement_factor)] .= range(xp[i], stop=xp[i+1], length=f_p_enhancement_factor+1)[1:end-1]
        end
        xp_new[end] = xp[end]
        xp, dfdx = xp_new, ForwardDiff.derivative.(Ref(f_pchip_spl), xp_new)
    end
    
    # create a pchip interpolator to that dxdp
    spl_dfdx = PCHIPInterpolation.Interpolator(xp, dfdx)
    dfp_dx_xmin = PCHIPInterpolation._derivative(spl_dfdx, Val(:begin), 1) # slope at beginning of interval 1
    dfp_dx_xmax = PCHIPInterpolation._derivative(spl_dfdx, Val(:end), length(xp)-1) # slope at end of interval N-1 (there are only N-1 intervals on N points)
    dfdx_min = spl_dfdx.ys[1]
    dfdx_max = spl_dfdx.ys[end]    
    ymin = fp[1]
    ymax = fp[end]

    xmin = xp[1]
    xmax = xp[end]
    
    # integrate it (only takes definite bounds so we can integrate between adjacent points and then cumsum?)
    return x -> begin 
        if x < xp[1] # integrate from x to x_0
            if bc == "error"
                error("Requested x is below the minimum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate")
            elseif bc == "extrapolate"
                x_0 = xp[1] # aka xmin
                Δx = x_0 - x
                # assume derivative from x to x_0 is a line  f'(x) = x-> dfdx_min - dfp_dx_xmin * Δx, aka the second derivative is continous, the derivative is smooth | integrate to get f(x), but go from x_0 to x (so negative of x to x_0) 
                return ymin - ( dfdx_min * Δx  - dfp_dx_xmin * ( + x_0 * Δx - (x_0^2 - x^2) / 2) ) #+ (ymin+xmin)/dfdx_min + ymin  # -(...) bc we integrate from x_0 to x then take the negative 
            else
                error("Unsupported bc option $bc")
            end
        elseif x > xp[end] # integrate from x to x_N
            if bc == "error"
                error("Requested x is above the maximum x of the spline but bc is set to error, use bc=\"extrapolate\" to extrapolate")
            elseif bc == "extrapolate"
                x_0 = xp[end] # aka xmax
                Δx = x - x_0
                ymax = ymin + PCHIPInterpolation.integrate(spl_dfdx, xp[1], xp[end]) # more accurate bc integration means you're slightly off on the right side, not sure if it's just numerical error or what
                return ymax + ( dfdx_max * Δx  + dfp_dx_xmax * ( (x^2 - x_0^2) / 2 - x_0 * Δx ) ) # -(...) bc we integrate from x_0 to x then take the negative 
            else
                error("Unsupported bc option $bc")
            end
        else
            return ymin + PCHIPInterpolation.integrate(spl_dfdx, xmin, x) # start at xmin,ymin and integrate to x
        end
    end

end