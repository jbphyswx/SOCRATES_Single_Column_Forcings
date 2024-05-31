
"""
This is a function to get the pressure and density profiles for a given flight number/forcing to use as the reference profiles in TC.jl
    We will:
        Read p, t
        Calculate ρ
        return a named tuple with p_c, p_f, ρ_c, ρ_f

SAM seems to hold ρ constant with time, so we can probably just rely on that
"""
function get_LES_reference_profiles(
    flight_number::Int;
    obs_or_ERA5 = "Obs"::Union{String, Symbol},
    new_zc::Union{Nothing, AbstractArray} = nothing,
    new_zf::Union{Nothing, AbstractArray} = nothing,
    # thermo_params,
)

    # initial conditions
    data = open_atlas_les_input(flight_number)

    return_values = (:p_c, :p_f, :ρ_c, :ρ_f)

    # if nothing is given, read in z from grid data to serve as new zf (plus 0 at sfc)
    # if only zf is given, calculate zc as midpoints between zf
    # if both zf and zc are given, use them as is (trust lol)
    # if only zc is given, throw an error since we don't want to do math to ensure we have repeated dzs to make that work


    # Setup new zc and zf based on input
    if isnothing(new_zc) && isnothing(new_zf)
        new_zf = [0; data[:grid_data]] # by default, the grid we read in becomes zf the way TC.jl is set up
        new_zc = (new_zf[1:end-1] .+ new_zf[2:end]) ./ 2
        # named tuple repeated for each in return values
        new_z = NamedTuple{return_values}((new_zc, new_zf, new_zc, new_zf))
    elseif isnothing(new_zc) && isa(new_zf, AbstractArray)
        new_zc = (new_zf[1:end-1] .+ new_zf[2:end]) ./ 2
        new_z = NamedTuple{return_values}((new_zc, new_zf, new_zc, new_zf))
    elseif isa(new_zc, AbstractArray) && isa(new_zf, AbstractArray)
        new_z = NamedTuple{return_values}((new_zc, new_zf, new_zc, new_zf))
    else
        error("You must provide either new_zc or new_zf, or both, but only providing new_zc is not allowed")
    end

    # resolve forcing keyword
    if obs_or_ERA5 ∈ ["Obs", :obs_data]
        forcing = :obs_data
        # data = data[(:obs_data, :ERA5_data)]
        data = data[(:obs_data,)] # drop ERA if we're doing obs cause we don't need it (11 has no obs either i think)        
    elseif obs_or_ERA5 ∈ ["ERA5", :ERA5_data]
        forcing = :ERA5_data
        data = data[(:ERA5_data,)] # drop obs if we're doing era5 cause we don't need it (11 has no obs either i think)        
    else
        if obs_or_ERA5 ∉ [:obs_data, :ERA5_data, "Obs", "ERA5"]
            error("obs_or_ERA5 must be either \"Obs\"/[:obs_data] or \"ERA5\"/[:ERA5_data]")
        else
            forcing = obs_or_ERA5
        end
    end


    LES_data = open_atlas_les_output(flight_number)[forcing]
    z = LES_data["z"][:]

    p = LES_data["p"][:] .* 100. # Pressure variations in SAM are under a milibar, so we can use the 1D p variable rather than the 2D PRES variable
    ps = LES_data["Ps"][1] * 100. # surface pressure
    ρ = LES_data["RHO"][:,1] # use t=0 as our reference
    # extrapolate to get ρs since that's not given (not calculating from first principles probably is safer too w/ uncertainty in q)
    ρs = pyinterp([ps], reverse(p), reverse(ρ); method=:Spline1D, bc="extrapolate")[1] # switch to increasing for interpolation

    p = [ps; p]
    ρ = [ρs; ρ]
    z = [0; z]

    p_c = pyinterp(new_zc, z, p; method=:Spline1D, bc="extrapolate")
    p_f = pyinterp(new_zf, z, p; method=:Spline1D, bc="extrapolate")
    ρ_c = pyinterp(new_zc, z, ρ; method=:Spline1D, bc="extrapolate")
    ρ_f = pyinterp(new_zf, z, ρ; method=:Spline1D, bc="extrapolate")

    return NamedTuple{return_values}((p_c, p_f, ρ_c, ρ_f))

end